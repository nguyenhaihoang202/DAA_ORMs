# ============================================================
# Frequentist ORM: HMP_2012 Gingival V35 SUBSET (Gamboa-Tuz)
# ============================================================
library(tidyverse)
library(TreeSummarizedExperiment)

# -----------------------------
# Function to run ORM with error handling
# -----------------------------
run_orm <- function(abundance, metadata, formula){
  mm <- model.matrix(formula, metadata) |> 
    cbind(abundance) |> 
    tibble::as_tibble() |> 
    dplyr::select(-"(Intercept)")
  
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  res <- data.frame(variable = vars, estimate = NA, se = NA, p_value = NA)
  
  fit_1 <- try(rms::orm(abundance ~ ., data = mm, maxiter = 100), silent = TRUE) 
  if (inherits(fit_1, "try-error")) {
    return(res |> tibble::rownames_to_column("variable"))
  }
  
  score_1 <- fit_1$stats["Score"]
  res <- data.frame(estimate = fit_1$coefficients[vars],
                    se = sqrt(diag(vcov(fit_1))[vars]),
                    p_value = NA)
  
  if (length(inds) > 1) {
    for (i in inds) {
      fit_0 <- try(rms::orm(abundance ~ ., data = mm[, -i], maxiter = 100), silent = TRUE)
      if (inherits(fit_0, "try-error")) {
        res$p_value[i] <- NA
      } else {
        score_0 <- fit_0$stats["Score"]
        res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
      }
    }
  } else {
    res$p_value <- as.numeric(1 - pchisq(score_1 - 0, df = 1))
  }
  
  return(res |> tibble::rownames_to_column("variable"))
}

# -----------------------------
# Load gingival SUBSET dataset
# -----------------------------
load("mb_datasets_gamboa_tuz.rds")

tse <- data_mbd_raw$HMP_2012_16S_gingival_V35_subset

# Keep only the two relevant body subsites
tse <- tse[, tse$body_subsite %in% c("subgingival_plaque", "supragingival_plaque")]

# Filter taxa by prevalence (adjust threshold as needed)
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)

# Extract counts and normalize to relative abundances
counts <- assay(tse) |> t()
rel_abs <- as.data.frame(counts / rowSums(counts))

# Metadata with group coding
meta <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(group = factor(body_subsite,
                        levels = c("subgingival_plaque", "supragingival_plaque"),
                        labels = c("subgingival", "supragingival")))

# -----------------------------
# Ground truth from annotations
# -----------------------------
taxonomy_df <- rowData(tse) |> as.data.frame() |> 
  tibble::rownames_to_column("taxon")

# Define ground truth (aerobic -> supragingival, anaerobic -> subgingival)
bvs <- taxonomy_df %>%
  filter(taxon_annotation %in% c("aerobic", "anaerobic")) %>%
  mutate(
    ground_truth = case_when(
      taxon_annotation == "aerobic"   ~ "supragingival",
      taxon_annotation == "anaerobic" ~ "subgingival",
      TRUE                            ~ "none"
    )
  ) %>%
  select(taxon, ground_truth)

# -----------------------------
# Run ORM per taxon
# -----------------------------
res <- rel_abs |> 
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = "taxon")

# -----------------------------
# Compare to ground truth
# -----------------------------
sl <- 0.1  # BH-FDR significance level

res2 <- res |> 
  left_join(bvs, by = "taxon") |> 
  drop_na() |> 
  mutate(
    q = p.adjust(p_value, method = "BH"),
    res = case_when(
      q < sl & estimate > 0 ~ "supragingival",
      q < sl & estimate < 0 ~ "subgingival",
      TRUE                  ~ "ns"
    ),
    correct   = res == ground_truth & res != "ns",
    incorrect = res != ground_truth & res != "ns"
  )

# -----------------------------
# Precision & Recall evaluation
# -----------------------------
eval_summary <- res2 %>%
  filter(!is.na(ground_truth)) %>%
  summarise(
    TP = sum(correct),
    FP = sum(incorrect),
    FN = sum(res == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) %>%
  mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN)
  )

print(eval_summary)

# ============================================================
# Bayesian hierarchical zero-inflated ordered logistic model
# ============================================================
library(cmdstanr)
library(posterior)

# 1) Encode abundances into ordinal categories with a structural zero bin
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x)
  z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L
  nx <- x[!z]
  if (length(nx) == 0) return(out)
  qs <- quantile(nx, probs = seq(0, 1, length.out = K_nonzero + 1), na.rm = TRUE, type = 8)
  eps <- 1e-9 * (max(nx, na.rm = TRUE) - min(nx, na.rm = TRUE) + 1)
  for (j in 2:length(qs)) if (qs[j] <= qs[j-1]) qs[j] <- qs[j-1] + eps
  out[!z] <- cut(nx, breaks = qs, include.lowest = TRUE, labels = FALSE) + 1L
  out
}

K_nonzero <- 3
# rel_abs is samples x taxa
rel_cats <- apply(rel_abs, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)  # samples x taxa
y_mat <- t(rel_cats)                                                           # taxa x samples
M <- nrow(y_mat); N <- ncol(y_mat); K <- max(y_mat)
taxa_names <- colnames(rel_abs)

# 2) Vectorize for Stan (group: 0 = subgingival, 1 = supragingival)
group_num <- ifelse(meta$group == "supragingival", 1L, 0L)
y_vec     <- as.vector(y_mat)                 # sample 1: all taxa, sample 2: all taxa, ...
group_vec <- rep(group_num, each = M)
taxon_idx <- rep(seq_len(M), times = N)
stan_data <- list(MN = M * N, M = M, K = K, y = y_vec, group = group_vec, taxon_idx = taxon_idx)

# 3) Compile and sample the hierarchical model
stan_file <- "hierarchical_model_test.stan"   
mod_hier  <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))
stan_init <- function() list(c_base = qlogis((1:(K-1))/K), delta_c = rep(0, M))

fit_hier <- mod_hier$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 500, iter_sampling = 500,
  seed = 1, init = stan_init
)

# 4) Summarize taxon-level group effects beta[m]
draws_hier <- fit_hier$draws("beta") |> as_draws_df()
res_bayes <- draws_hier |>
  dplyr::select(dplyr::starts_with("beta[")) |>
  tidyr::pivot_longer(everything(), names_to = "taxon_ix", values_to = "beta") |>
  dplyr::group_by(taxon_ix) |>
  dplyr::summarise(
    est   = median(beta),
    lwr95 = quantile(beta, 0.025),
    upr95 = quantile(beta, 0.975),
    p_gt0 = mean(beta > 0),
    p_lt0 = mean(beta < 0),
    .groups = "drop"
  ) |>
  dplyr::mutate(ix = as.integer(gsub("beta\\[|\\]", "", taxon_ix)),
                taxon = taxa_names[ix]) |>
  dplyr::select(taxon, est, lwr95, upr95, p_gt0, p_lt0)

# 5) Posterior decision rule and join with ground truth
pp_thresh <- 0.995  
res_bayes2 <- res_bayes |>
  dplyr::mutate(
    call = dplyr::case_when(
      p_gt0 >= pp_thresh ~ "supragingival",     # group = 1
      p_lt0 >= pp_thresh ~ "subgingival",       # group = 0
      TRUE                  ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  tidyr::drop_na()

# 6) Precision and recall for Bayesian model
eval_bayes <- res_bayes2 |>
  dplyr::summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "bayesian"
  )

# 7) Precision and recall for ORM (reuse res2)
eval_orm <- res2 |>
  dplyr::summarise(
    TP = sum(correct),
    FP = sum(incorrect),
    FN = sum(res == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "orm"
  )

# 8) Side by side summary
bind_rows(eval_orm, eval_bayes) |> print()

# The 0.995 threshold is too strict for this dataset -> no calls for the bayesian model.
----------------------------------------------------------

# 1) How concentrated are posterior sign probabilities?
# If most p_gt0 and p_lt0 are below 0.995, a 0.995 threshold will call nothing.
summary(res_bayes$p_gt0)
summary(res_bayes$p_lt0)

# Fraction of taxa exceeding common cutoffs on either side
frac_ge <- tibble(
  cutoff = c(0.95, 0.975, 0.99, 0.995),
  frac_any = c(
    mean(res_bayes$p_gt0 >= 0.95  | res_bayes$p_lt0 >= 0.95),
    mean(res_bayes$p_gt0 >= 0.975 | res_bayes$p_lt0 >= 0.975),
    mean(res_bayes$p_gt0 >= 0.99  | res_bayes$p_lt0 >= 0.99),
    mean(res_bayes$p_gt0 >= 0.995 | res_bayes$p_lt0 >= 0.995)
  )
)
print(frac_ge)
# frac_any at 0.995 is 0, the threshold is simply too strict for this dataset.

# 2) Sanity check the direction:
# Aerobic should lean supragingival (positive beta), anaerobic should lean subgingival (negative beta).
dir_check <- res_bayes %>%
  dplyr::left_join(bvs, by = "taxon") %>%
  dplyr::filter(ground_truth %in% c("supragingival","subgingival")) %>%
  dplyr::group_by(ground_truth) %>%
  dplyr::summarise(
    mean_est   = mean(est, na.rm = TRUE),
    median_est = median(est, na.rm = TRUE),
    mean_p_gt0 = mean(p_gt0, na.rm = TRUE),
    mean_p_lt0 = mean(p_lt0, na.rm = TRUE),
    n = dplyr::n()
  )
print(dir_check)
# Hierarchical shrinkage in the model is pulling many supragingival betas toward the global mean.


# 4) Threshold sweep to see precision and recall at more practical Bayesian cutoffs
sweep_pp <- function(th) {
  res_bayes %>%
    dplyr::mutate(
      call = dplyr::case_when(
        p_gt0 >= th ~ "supragingival",
        p_lt0 >= th ~ "subgingival",
        TRUE        ~ "ns"
      )
    ) %>%
    dplyr::left_join(bvs, by = "taxon") %>%
    tidyr::drop_na() %>%
    dplyr::summarise(
      TP = sum(call == ground_truth & call != "ns"),
      FP = sum(call != ground_truth & call != "ns"),
      FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
    ) %>%
    dplyr::mutate(
      precision = TP / (TP + FP),
      recall    = TP / (TP + FN),
      pp_thresh = th
    )
}

pp_grid <- c(0.95, 0.975, 0.99, 0.995)
bayes_sweep <- purrr::map_dfr(pp_grid, sweep_pp)
print(bayes_sweep)
# Comment: this shows how recall grows as we relax pp_thresh from 0.995 down to 0.975 or 0.95.
# Should pick a threshold that balances recall against false positives for this dataset.

# ------------------------------------------------------------
# REPLACE with a more realistic threshold:
# ------------------------------------------------------------
pp_thresh <- 0.975  # was 0.995; 0.975 is often a good tradeoff on weaker-signal datasets

res_bayes2 <- res_bayes %>%
  dplyr::mutate(
    call = dplyr::case_when(
      p_gt0 >= pp_thresh ~ "supragingival",
      p_lt0 >= pp_thresh ~ "subgingival",
      TRUE               ~ "ns"
    )
  ) %>%
  dplyr::left_join(bvs, by = "taxon") %>%
  tidyr::drop_na()

# Recompute precision and recall with the new threshold
eval_bayes <- res_bayes2 %>%
  dplyr::summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) %>%
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "bayesian"
  )
print(eval_bayes)

# This seems to be more realistic, can you help me to review?