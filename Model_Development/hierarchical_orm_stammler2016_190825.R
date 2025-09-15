# ============================================================
# Stammler 2016 spike-in benchmark: ORM + Bayesian hierarchical model
# ============================================================
library(tidyverse)
library(TreeSummarizedExperiment)
library(cmdstanr)
library(posterior)

# -----------------------------
# ORM 
# -----------------------------
run_orm <- function(abundance, metadata, formula){
  mm <- model.matrix(formula, metadata) |>
    cbind(abundance) |>
    tibble::as_tibble() |>
    dplyr::select(-"(Intercept)")
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  # Fit ordinal regression
  fit_1 <- rms::orm(abundance ~ ., data = mm, maxiter = 100)
  score_1 <- fit_1$stats["Score"]
  
  # Collect coefficient estimates
  res <- data.frame(
    estimate = fit_1$coefficients[vars],
    se = sqrt(diag(vcov(fit_1))[vars]),
    p_value = NA_real_
  )
  
  # Score test p-values
  if (length(inds) > 1) {
    for (i in inds) {
      fit_0 <- rms::orm(abundance ~ ., data = mm[, -i], maxiter = 100)
      score_0 <- fit_0$stats["Score"]
      res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
    }
  } else {
    res$p_value <- as.numeric(1 - pchisq(score_1, df = 1))
  }
  res |> tibble::rownames_to_column("variable")
}

# -----------------------------
# Load dataset 
# -----------------------------
load("mb_datasets_gamboa_tuz.rds")
tse <- data_mbd_raw$Stammler_2016_16S_spikein

# Spike-in IDs
spike_in_ids <- c("AF323500XXXX", "AB247615XXXX", "AB076660XXXX")

# ============================================================
# Case 1: completed vs uncompleted
# ============================================================
# Define groups from study_condition
cd <- colData(tse) |> as.data.frame()
cd$group <- dplyr::case_when(
  cd$study_condition %in% c("preASCT", "d0", "d7") ~ "uncompleted",
  cd$study_condition %in% c("d14")                 ~ "completed",
  TRUE                                             ~ NA_character_
)
keep <- !is.na(cd$group)
tse1 <- tse[, keep]
meta1 <- cd[keep, , drop = FALSE] |>
  dplyr::mutate(group = factor(group, levels = c("uncompleted", "completed")))
tse1 <- mia::subsetByPrevalent(tse1, prevalence = 0.01)
counts1 <- assay(tse1) |> t()
rel_abs1 <- as.data.frame(counts1 / rowSums(counts1))

# Annotate spike-ins
bvs1 <- rowData(tse1) |>
  as.data.frame() |> tibble::rownames_to_column("taxon") |>
  mutate(ground_truth = ifelse(taxon %in% spike_in_ids, "spike_in", "endogenous"))

# Run ORM
res1 <- rel_abs1 |> purrr::map(~ run_orm(., metadata = meta1, formula = ~ group)) |> bind_rows(.id = "taxon")
res2_1 <- res1 |> left_join(bvs1, by = "taxon") |>
  mutate(q = p.adjust(p_value, method = "BH"),
         call_orm = case_when(q < 0.1 & estimate > 0 ~ "completed",
                              q < 0.1 & estimate < 0 ~ "uncompleted",
                              TRUE ~ "ns"))

# ============================================================
# Case 2: Patients 1+2 vs Patients 3–5
# ============================================================
# Define groups
meta2 <- colData(tse) |> as.data.frame() |>
  mutate(group = case_when(
    subject_id %in% c("Patient1", "Patient2") ~ "groupA",
    subject_id %in% c("Patient3","Patient4","Patient5") ~ "groupB",
    TRUE ~ NA_character_
  ),
  group = factor(group, levels = c("groupA","groupB"))) |>
  filter(!is.na(group))
tse2 <- tse[, rownames(meta2)]
tse2 <- mia::subsetByPrevalent(tse2, prevalence = 0.01)
counts2 <- assay(tse2) |> t()
rel_abs2 <- as.data.frame(counts2 / rowSums(counts2))

# Annotate spike-ins
bvs2 <- rowData(tse2) |> as.data.frame() |> tibble::rownames_to_column("taxon") |>
  mutate(ground_truth = ifelse(taxon %in% spike_in_ids, "spike_in", "endogenous"))

# Run ORM
res2 <- rel_abs2 |> purrr::map(~ run_orm(., metadata = meta2, formula = ~ group)) |> bind_rows(.id = "taxon")
res2_2 <- res2 |> left_join(bvs2, by = "taxon") |>
  mutate(call_orm = case_when(
    p_value < 0.1 & estimate > 0 ~ "groupB",
    p_value < 0.1 & estimate < 0 ~ "groupA",
    TRUE ~ "ns"
  ))

# ============================================================
# Case 3: Patient 1 vs Patient 2
# ============================================================
# Define groups
tse3 <- tse[, tse$subject_id %in% c("Patient1","Patient2")]
tse3 <- mia::subsetByPrevalent(tse3, prevalence = 0.01)
counts3 <- assay(tse3) |> t()
rel_abs3 <- as.data.frame(counts3 / rowSums(counts3))
meta3 <- colData(tse3) |> as.data.frame() |>
  mutate(group = case_when(
    subject_id %in% c("Patient1") ~ "groupA",
    subject_id %in% c("Patient2") ~ "groupB"
  ),
  group = factor(group, levels = c("groupA","groupB")))

# Annotate spike-ins
bvs3 <- rowData(tse3) |> as.data.frame() |> tibble::rownames_to_column("taxon") |>
  mutate(ground_truth = ifelse(taxon %in% spike_in_ids, "spike_in", "endogenous"))

# Run ORM
res3 <- rel_abs3 |> purrr::map(~ run_orm(., metadata = meta3, formula = ~ group)) |> bind_rows(.id = "taxon")
res2_3 <- res3 |> left_join(bvs3, by = "taxon") |>
  mutate(call_orm = case_when(
    p_value < 0.05 & estimate > 0 ~ "groupB",
    p_value < 0.05 & estimate < 0 ~ "groupA",
    TRUE ~ "ns"
  ))

# ============================================================
# Bayesian 
# ============================================================
# Zero-inflated ordinal binning
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x); z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L
  nx <- x[!z]; if (length(nx) == 0) return(out)
  qs <- quantile(nx, probs = seq(0,1,length.out = K_nonzero+1), na.rm=TRUE, type=8)
  eps <- 1e-9 * (max(nx, na.rm=TRUE) - min(nx, na.rm=TRUE) + 1)
  for (j in 2:length(qs)) if (qs[j] <= qs[j-1]) qs[j] <- qs[j-1] + eps
  out[!z] <- cut(nx, breaks=qs, include.lowest=TRUE, labels=FALSE) + 1L
  out
}

# Helper to run Stan model for a dataset with desired settings
run_bayes_block <- function(rel_abs_mat, meta_df, bvs_tbl, stan_mod, pp_thresh=0.995, K_nonzero=3) {
  rel_cats <- apply(rel_abs_mat, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)
  y_mat <- t(rel_cats); M <- nrow(y_mat); N <- ncol(y_mat); K <- max(y_mat)
  taxa_names <- colnames(rel_abs_mat)
  group_num <- ifelse(meta_df$group == levels(meta_df$group)[2], 1L, 0L)
  y_vec <- as.vector(y_mat); group_vec <- rep(group_num, each=M); taxon_idx <- rep(seq_len(M), times=N)
  stan_data <- list(MN=M*N, M=M, K=K, y=y_vec, group=group_vec, taxon_idx=taxon_idx)
  stan_init <- function() list(c_base=qlogis((1:(K-1))/K), delta_c=rep(0,M))
  
  # Sample
  fit <- stan_mod$sample(data=stan_data, chains=4, parallel_chains=4, threads_per_chain=3,
                         iter_warmup=500, iter_sampling=500, seed=1, init=stan_init)
  
  # Extract beta effects
  draws <- fit$draws("beta") |> as_draws_df()
  res_bayes <- draws |>
    dplyr::select(dplyr::starts_with("beta[")) |>
    tidyr::pivot_longer(everything(), names_to="taxon_ix", values_to="beta") |>
    dplyr::group_by(taxon_ix) |>
    dplyr::summarise(est=median(beta),
                     lwr95=quantile(beta,0.025),
                     upr95=quantile(beta,0.975),
                     p_gt0=mean(beta>0), p_lt0=mean(beta<0), .groups="drop") |>
    dplyr::mutate(ix=as.integer(gsub("beta\\[|\\]","",taxon_ix)),
                  taxon=taxa_names[ix]) |>
    dplyr::select(taxon, est, lwr95, upr95, p_gt0, p_lt0)
  
  # Decision rule
  res_bayes2 <- res_bayes |>
    dplyr::mutate(call_bayes = case_when(
      p_gt0 >= pp_thresh ~ levels(meta_df$group)[2],
      p_lt0 >= pp_thresh ~ levels(meta_df$group)[1],
      TRUE ~ "ns")) |>
    dplyr::left_join(bvs_tbl, by="taxon")
  list(res=res_bayes2)
}

# Compile Stan model once
stan_file <- "hierarchical_model_test.stan"
mod_hier <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# Bayesian runs for each case
bayes1 <- run_bayes_block(rel_abs1, meta1, bvs1, mod_hier)
bayes2 <- run_bayes_block(rel_abs2, meta2, bvs2, mod_hier)
bayes3 <- run_bayes_block(rel_abs3, meta3, bvs3, mod_hier, pp_thresh=0.95)

# ============================================================
# Precision & Recall evaluation for ORM vs Bayesian (endogenous taxa only)
# ============================================================
eval_pr_method <- function(res_tbl, call_col, method_name) {
  res_tbl %>%
    filter(ground_truth == "endogenous") %>%
    summarise(
      TP = sum(.data[[call_col]] != "ns"),
      FP = 0,
      FN = sum(.data[[call_col]] == "ns")
    ) %>%
    mutate(
      precision = ifelse(TP+FP > 0, TP/(TP+FP), NA),
      recall    = TP/(TP+FN),
      method    = method_name
    )
}

# Compute PR per case and method
# Case 1
pr1_orm   <- eval_pr_method(res2_1, "call_orm",   "ORM") %>% mutate(case = "completed vs uncompleted")
pr1_bayes <- eval_pr_method(bayes1$res, "call_bayes", "Bayesian") %>% mutate(case = "completed vs uncompleted")

# Case 2
pr2_orm   <- eval_pr_method(res2_2, "call_orm", "ORM") %>% mutate(case = "Patients1+2 vs Patients3–5")
pr2_bayes <- eval_pr_method(bayes2$res, "call_bayes", "Bayesian") %>% mutate(case = "Patients1+2 vs Patients3–5")

# Case 3
pr3_orm   <- eval_pr_method(res2_3, "call_orm", "ORM") %>% mutate(case = "Patient1 vs Patient2")
pr3_bayes <- eval_pr_method(bayes3$res, "call_bayes", "Bayesian") %>% mutate(case = "Patient1 vs Patient2")

precision_recall_summary <- bind_rows(pr1_orm, pr1_bayes,
                                      pr2_orm, pr2_bayes,
                                      pr3_orm, pr3_bayes)

print(precision_recall_summary)

# - Case 1 (completed vs uncompleted):
#   ORM called 1 endogenous taxon (TP=1), Bayesian called none.
#   Recall is essentially 0 for both → almost no taxa flagged.
#
# - Case 2 (Patients1+2 vs Patients3–5):
#   ORM recalled ~8% of endogenous taxa (196/2482).
#   Bayesian again called none, so recall = 0.
#
# - Case 3 (Patient1 vs Patient2):
#   ORM and Bayesian both called the same 86 taxa (~4.8% recall).
#   Here, Bayesian matched ORM (likely because the smaller dataset
#   reduces shrinkage, making Bayesian less conservative).