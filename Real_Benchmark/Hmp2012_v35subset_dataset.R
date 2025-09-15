# ============================================================
# Frequentist ORM: HMP_2012 Gingival V35 SUBSET (Gamboa-Tuz). Here I run 1000 iter warmup and 1000 iter sampling; take 4 hours for the ordinal model.
# Should run 500 + 500 for testing
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

# Filter taxa by prevalence
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
eval_orm <- res2 %>%
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

print(eval_orm)

# ============================================================
# Bayesian model
# ============================================================
library(cmdstanr)
library(posterior)

# ---- Ordinal binning (4 categories; zeros fall into the lowest bin) ----
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x)
  z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L
  nx <- x[!z]
  if (length(nx) == 0) return(out)
  qs <- quantile(nx, probs = seq(0, 1, length.out = K_nonzero + 1), na.rm = TRUE, type = 8)
  eps <- 1e-9 * (max(nx, na.rm=TRUE) - min(nx, na.rm=TRUE) + 1)
  for (j in 2:length(qs)) if (qs[j] <= qs[j-1]) qs[j] <- qs[j-1] + eps
  out[!z] <- cut(nx, breaks = qs, include.lowest = TRUE, labels = FALSE) + 1L
  out
}

# ---- Build ordinal data for the new model (K = 4 total categories) ----
K_nonzero <- 3
rel_cats  <- apply(rel_abs, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)  # samples x taxa
y_mat     <- t(rel_cats)                        # taxa x samples
M         <- nrow(y_mat)                        # taxa
N         <- ncol(y_mat)                        # samples
K_ord     <- max(y_mat)                         # ordinal categories (should be 4)
taxa_names <- colnames(rel_abs)

# Group coding: 0 = subgingival, 1 = supragingival (model internally centers by -0.5)
group_num <- ifelse(meta$group == "supragingival", 1L, 0L)
y_vec     <- as.vector(y_mat)                   # sample1: all taxa, sample2: all taxa, ...
group_vec <- rep(group_num, each = M)
taxon_idx <- rep(seq_len(M), times = N)

stan_data_ord <- list(
  MN = M * N, M = M, K = K_ord,
  y = y_vec, group = group_vec, taxon_idx = taxon_idx
)

# ---- Fit the newest model with cmdstanr ----
mod_bayes <- cmdstan_model("bayes_model_210825.stan", cpp_options = list(stan_threads = TRUE))

fit_bayes <- mod_bayes$sample(
  data = stan_data_ord,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 1000, iter_sampling = 1000,
  seed = 1
)

# ---- Summarize per-taxon beta and convert to calls ----
draws_bayes <- fit_bayes$draws("beta") |> as_draws_df()
res_bayes   <- draws_bayes |>
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

pp_thresh_bayes <- 0.975   # a realistic posterior sign threshold for this dataset
res_bayes_calls <- res_bayes |>
  dplyr::mutate(
    call = dplyr::case_when(
      p_gt0 >= pp_thresh_bayes ~ "supragingival",
      p_lt0 >= pp_thresh_bayes ~ "subgingival",
      TRUE                     ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  tidyr::drop_na()

eval_bayes210825 <- res_bayes_calls |>
  dplyr::summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "bayes_210825"
  )

# ============================================================
# BayDiPAn model
# ============================================================
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- Build binary presence/absence + covariates for prevalence model ----
# counts: samples x taxa (already available above)
y_bin_KN <- t( (counts > 0) * 1L )                # K x N (taxa x samples)
x_vec    <- as.numeric(meta$group == "supragingival")  # N-length vector in {0,1}
reads_vec <- as.numeric(scale(log10(rowSums(counts)), center = TRUE, scale = FALSE)) # centered

stan_data_sup <- list(
  N = nrow(counts),
  K = ncol(counts),
  y = y_bin_KN,
  x = x_vec,
  reads = reads_vec
)

sm_sup   <- rstan::stan_model(file = "stan_blogrrall_120625.stan")
fit_sup  <- rstan::sampling(sm_sup, data = stan_data_sup,
                            chains = 4, iter = 2000, warmup = 1000)

# ---- Summarize supervisor beta and convert to calls ----
post_beta <- rstan::extract(fit_sup, pars = "beta")$beta  # iterations x K
taxa_sup  <- colnames(counts)                             # taxa names align with K
res_sup   <- tibble::tibble(
  taxon = taxa_sup,
  est   = apply(post_beta, 2, stats::median),
  lwr95 = apply(post_beta, 2, stats::quantile, 0.025),
  upr95 = apply(post_beta, 2, stats::quantile, 0.975),
  p_gt0 = colMeans(post_beta > 0),
  p_lt0 = colMeans(post_beta < 0)
)

pp_thresh_sup <- 0.975
res_sup_calls <- res_sup |>
  dplyr::mutate(
    call = dplyr::case_when(
      p_gt0 >= pp_thresh_sup ~ "supragingival",
      p_lt0 >= pp_thresh_sup ~ "subgingival",
      TRUE                   ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  tidyr::drop_na()

eval_supervisor <- res_sup_calls |>
  dplyr::summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "supervisor_120625"
  )

# ============================================================
# Side-by-side summary with ORM
# ============================================================
side_by_side <- dplyr::bind_rows(
  eval_orm,
  eval_bayes210825,
  eval_supervisor
)
print(side_by_side)

# It seems like frequentist ORM perform the best in this case? But both bayesian models seems feasible.