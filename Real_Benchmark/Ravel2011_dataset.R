library(tidyverse)
library(TreeSummarizedExperiment)

# Function to run ORM
run_orm <- function(abundance, metadata, formula){
  
  # Create the design matrix of the model
  mm <- model.matrix(formula, metadata) |> 
    cbind(abundance) |> 
    tibble::as_tibble() |> 
    dplyr::select(-"(Intercept)")
  
  # Get the indices and names of the variables in the model matrix
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  # Fit ordinal regression model.
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  
  # Extract the score statistic that is used to calculate p-values
  score_1 <- fit_1$stats["Score"]
  
  # Extract the log odds estimate and its standard error for each variable
  res <- data.frame(estimate = fit_1$coefficients[vars],
                    se = sqrt(diag(vcov(fit_1))[vars]),
                    p_value = NA)
  
  # Calculate the p-value based on score test for each variable
  if(length(inds) > 1){
    
    for(i in inds){
      
      #Fit the "null model" (the model without the variable of interest)
      fit_0 <- rms::orm(abundance ~ ., data = mm[, -i])
      score_0 <- fit_0$stats["Score"]
      
      # p-value is based on the difference of score statistics. Under the null
      # hypothesis it follows the chi squared distribution with 1 degree of
      # freedom. 
      res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
    }
    
  }else{
    
    # If there is only one variable in the model, the "null model" would include
    # only intercepts. The score statistic for such null model would be zero.
    res$p_value <- as.numeric(1 - pchisq(score_1 - 0, df = 1))
  }
  
  return(res |> tibble::rownames_to_column("variable"))
}


# Load data and prepare data----------------------------------------------------
load('mb_datasets_gamboa_tuz.rds')

# Extract the Ravel 2011 BV dataset
tse <- data_mbd_raw$Ravel_2011_16S_BV
tse$study_condition
# Filter the object to include only healthy and BV samples
tse <- tse[, tse$study_condition %in% c('healthy', 'bacterial_vaginosis')]

# Filter out taxa with prevalence < 0.01
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)


# Extract counts from the TSE object
counts <- tse |> assay() |> t()

# Relative abundances
rel_abs <- as.data.frame(counts / rowSums(counts))

# Metadata
meta <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(group = factor(study_condition,
                        levels = c('healthy', 'bacterial_vaginosis'),
                        labels = c('healthy', 'bv')))


# BV signatures. These define the "ground truth".
# BV associated bacteria should be more abundant in women with bacterial vaginosis
# HV associated bacteria should be more abundant in healthy women
bvs <- rowData(tse) |> 
  as.data.frame() |> 
  rownames_to_column('taxon') |> 
  mutate(ground_truth = case_when(taxon_annotation == 'bv-associated' ~ 'bv',
                                  taxon_annotation == 'hv-associated' ~ 'healthy',
                        T ~ 'none')) |>
  select(taxon, ground_truth)


# Run DAA-----------------------------------------------------------------------

res <- rel_abs |> 
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')


sl <- .10 # Significance level

# Compare the results against the ground truth
# Positive estimate indicates that the taxon is more abundant in women with BV   
res2 <- res |> 
  left_join(bvs, by = 'taxon') |> 
  mutate(q = p.adjust(p_value, method = 'BH'),
         res = case_when(q < sl & estimate > 0 ~ 'bv',
                         q < sl & estimate < 0 ~ 'healthy',
                         T ~ 'ns'),
         correct = ifelse(res == ground_truth, T, F),
         incorrect = ifelse(res == 'bv' & ground_truth == 'healthy' |
                      res == 'healthy' & ground_truth == 'bv', T, F))

# Summarize the results against the ground truth
res2 |> summarize(correct = sum(correct),
                  incorrect = sum(incorrect),
                  n = sum(ground_truth != 'none'))

# 20 correctly identified taxa among 28 taxa with known ground truth
# No incorrect findings

# =========================
# Bayesian models
# =========================
library(cmdstanr)
library(posterior)
library(rstan)

# ---- Ordinal binning for bayes_model_210825 (4 categories; zeros -> bin 1) ----
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

# ---- Build ordinal data (samples x taxa -> taxa x samples) ----
K_nonzero <- 3
rel_cats  <- apply(rel_abs, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)  # samples x taxa
y_mat     <- t(rel_cats)                                                        # taxa x samples
M         <- nrow(y_mat)
N         <- ncol(y_mat)
K_ord     <- max(y_mat)
taxa_names <- colnames(rel_abs)

# Group coding for Bayesian model: 0 = healthy, 1 = bv  (model centers internally at -0.5/+0.5)
group_num <- ifelse(meta$group == "bv", 1L, 0L)
y_vec     <- as.vector(y_mat)                     # sample1: all taxa, sample2: all taxa, ...
group_vec <- rep(group_num, each = M)
taxon_idx <- rep(seq_len(M), times = N)

stan_data_ord <- list(
  MN = M * N, M = M, K = K_ord,
  y = y_vec, group = group_vec, taxon_idx = taxon_idx
)

# ---- Fit the new model ----
mod_bayes <- cmdstan_model("bayes_model_210825.stan", cpp_options = list(stan_threads = TRUE))
fit_bayes <- mod_bayes$sample(
  data = stan_data_ord,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 500, iter_sampling = 500,
  seed = 1
)

# ---- Summarize β per taxon and call BV/healthy with posterior sign threshold ----
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

pp_thresh_bayes <- 0.975
res_bayes_calls <- res_bayes |>
  dplyr::mutate(
    call = dplyr::case_when(
      p_gt0 >= pp_thresh_bayes ~ "bv",
      p_lt0 >= pp_thresh_bayes ~ "healthy",
      TRUE                     ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  dplyr::mutate(
    correct   = (call == ground_truth) & (ground_truth %in% c("bv","healthy")),
    incorrect = (call != "ns") & (ground_truth %in% c("bv","healthy")) & (call != ground_truth)
  )

eval_bayes210825 <- res_bayes_calls |>
  dplyr::filter(ground_truth %in% c("bv","healthy")) |>
  dplyr::summarise(
    TP = sum(correct, na.rm = TRUE),
    FP = sum(incorrect, na.rm = TRUE),
    FN = sum(call == "ns", na.rm = TRUE)
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "bayes_210825"
  )

# ============================================================
# BayDIPAN model
# ============================================================
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Presence/absence and covariates (supervisor model works on prevalence)
y_bin_KN <- t( (counts > 0) * 1L )                                   # K x N (taxa x samples)
x_vec     <- as.numeric(meta$group == "bv")                           # group in {0,1}
reads_vec <- as.numeric(scale(log10(rowSums(counts)), center = TRUE, scale = FALSE))  # centered

stan_data_sup <- list(
  N = nrow(counts),
  K = ncol(counts),
  y = y_bin_KN,
  x = x_vec,
  reads = reads_vec
)

sm_sup  <- rstan::stan_model(file = "stan_blogrrall_120625.stan")
fit_sup <- rstan::sampling(sm_sup, data = stan_data_sup,
                           chains = 4, iter = 1000, warmup = 500,
                           seed = 1)

post_beta <- rstan::extract(fit_sup, pars = "beta")$beta  # iterations x K
taxa_sup  <- colnames(counts)

res_sup <- tibble::tibble(
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
      p_gt0 >= pp_thresh_sup ~ "bv",
      p_lt0 >= pp_thresh_sup ~ "healthy",
      TRUE                   ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  dplyr::mutate(
    correct   = (call == ground_truth) & (ground_truth %in% c("bv","healthy")),
    incorrect = (call != "ns") & (ground_truth %in% c("bv","healthy")) & (call != ground_truth)
  )

eval_supervisor <- res_sup_calls |>
  dplyr::filter(ground_truth %in% c("bv","healthy")) |>
  dplyr::summarise(
    TP = sum(correct, na.rm = TRUE),
    FP = sum(incorrect, na.rm = TRUE),
    FN = sum(call == "ns", na.rm = TRUE)
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "supervisor_120625"
  )

# ============================================================
# Side-by-side summary with ORM
# ============================================================
eval_orm <- res2 |>
  dplyr::filter(ground_truth %in% c("bv","healthy")) |>
  dplyr::summarise(
    TP = sum(correct),
    FP = sum(incorrect),
    FN = sum(res == "ns")
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "orm"
  )

dplyr::bind_rows(eval_orm, eval_bayes210825, eval_supervisor) |> print()

# Per-taxon comparison table
compare_calls <- res2 |>
  dplyr::filter(ground_truth %in% c("bv","healthy")) |>
  dplyr::select(taxon, ground_truth, orm_call = res, orm_q = q, orm_est = estimate) |>
  dplyr::left_join(res_bayes_calls |> dplyr::select(taxon, bayes_call = call, bayes_pp = p_gt0, bayes_est = est),
                   by = "taxon") |>
  dplyr::left_join(res_sup_calls   |> dplyr::select(taxon, sup_call = call, sup_pp = p_gt0, sup_est = est),
                   by = "taxon")
print(compare_calls)

# The results seems good with both the ordinal and BayDIPAn model, compared to ORM. But in this case, cant see the main advantage of ordinal model...
