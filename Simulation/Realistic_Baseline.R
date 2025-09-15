library(tidyverse)
library(cmdstanr)
library(posterior)

# Vandeputte dataset
load('mean_vars_vandeputte_250625.rds')
load('cors_vandeputte_250625.rds')

n <- 100
ALPHA <- 0.10
Z90   <- qnorm(1 - ALPHA/2)   # ~1.64485
n_taxa <- nrow(mean_vars)
groups <- rep(0:1, each = n / 2)

# Compile Stan models ONCE 
mod_bziolr <- cmdstan_model("BZIOLR.stan",  cpp_options = list(stan_threads = TRUE))
mod_balor <- cmdstan_model("BALOR.stan",  cpp_options = list(stan_threads = TRUE))

# ---------- Loop over seeds ----------
for (SEED_MAIN in 25:34) {
  set.seed(SEED_MAIN)
  message(">>> Running REALISTIC replicate SEED_MAIN = ", SEED_MAIN)
  
  # ---------------- Data generation ----------------
  # Latent MVN with Vandeputte correlations
  mvn_latent <- MASS::mvrnorm(n, mu = rep(0, n_taxa), Sigma = cors)
  
  # Realistic effects: NO exact zeros, centered negative (sign-only FDR in spirit)
  effects <- rnorm(n_taxa, mean = -0.9, sd = 0.6)
  names(effects) <- paste0("taxon_", 1:n_taxa)
  
  # Apply effects to case group on latent scale
  effects_matrix <- (cbind(replicate(n, effects)) |> t()) * groups
  mvn_latent_plus_effect <- mvn_latent + effects_matrix
  
  # Gamma marginals (Vandeputte params)
  abd_means <- mean_vars$mean
  abd_vars  <- mean_vars$var
  
  gamma_shapes <- abd_means ^ 2 / abd_vars
  m_gamma_shapes <- matrix(rep(gamma_shapes, each = n), nrow = n, byrow = FALSE)
  gamma_scales <- abd_vars / abd_means
  m_gamma_scales <- matrix(rep(gamma_scales, each = n), nrow = n, byrow = FALSE)
  mvu_latent <- pnorm(mvn_latent_plus_effect)
  true_abundances <- qgamma(mvu_latent,
                            shape = m_gamma_shapes,
                            scale = m_gamma_scales)
  
  # Taxon-wise detection biases → closure → Poisson counts
  taxonwise_biases  <- exp(rnorm(n_taxa, mean = 0, sd = 1))
  biased_abundances <- t(taxonwise_biases * t(true_abundances))
  relative_biased_abundances <- biased_abundances / rowSums(biased_abundances)
  library_sizes <- 10 ^ rnorm(n, mean = 4.0, sd = 0.20) |> round()
  lambdas <- library_sizes * relative_biased_abundances
  counts <- rpois(n * n_taxa, lambda = lambdas) |> matrix(nrow = n, ncol = n_taxa)
  colnames(counts) <- names(effects)
  
  # ---------------- DAA baselines (consistent with other scenarios) ----------------
  rel_abs <- counts / rowSums(counts)
  rel_abundances <- rel_abs[, colSums(rel_abs > 0) >= 5, drop = FALSE]
  counts2        <- counts[, colSums(counts > 0) >= 5, drop = FALSE]
  
  run_orm <- function(abundance, metadata, formula) {
    mm <- model.matrix(formula, metadata) |>
      cbind(abundance) |>
      as_tibble() |>
      select(-"(Intercept)")
    inds <- 1:(ncol(mm) - 1)
    vars <- colnames(mm)[inds]
    fit_1 <- rms::orm(abundance ~ ., data = mm)
    score_1 <- fit_1$stats["Score"]
    res_freq <- data.frame(
      estimate = fit_1$coefficients[vars],
      se = sqrt(diag(vcov(fit_1))[vars]),
      p_value = NA_real_
    )
    if (length(inds) > 1) {
      for (i in inds) {
        fit_0 <- rms::orm(abundance ~ ., data = mm[, -i])
        score_0 <- fit_0$stats["Score"]
        res_freq$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
      }
    } else {
      res_freq$p_value <- as.numeric(1 - pchisq(score_1, df = 1))
    }
    res_freq |> rownames_to_column("variable")
  }
  
  meta <- data.frame(group = factor(groups, labels = c("control", "case")))
  
  res_orm <- rel_abundances |>
    as.data.frame() |>
    purrr::map(~ run_orm(., metadata = meta, formula = ~ group)) |>
    bind_rows(.id = 'taxon')
  
  make_physeq <- function(counts, meta){
    otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
    meta_data <- phyloseq::sample_data(meta)
    phyloseq::phyloseq(otu, meta_data)
  }
  
  obj_ancombc <- ANCOMBC::ancombc(
    data = make_physeq(counts2, meta),
    formula = "group",
    p_adj_method = "fdr",
    prv_cut = 0
  )
  
  obj_linda <- LinDA::linda(
    otu.tab = t(counts2),
    meta = meta,
    formula = "~ group"
  )
  
  # Truth table
  true_effects <- tibble(taxon = names(effects), true_effect = effects)
  
  # ---- ANCOM-BC  ----
  res_ancombc <- tibble(
    taxon    = obj_ancombc$res$lfc$taxon,
    estimate = obj_ancombc$res$lfc$groupcase,
    se       = obj_ancombc$res$se$groupcase,
    p        = obj_ancombc$res$p_val$groupcase,
    q        = obj_ancombc$res$q_val$groupcase
  ) |>
    left_join(true_effects, by = "taxon") |>
    mutate(significant = q < ALPHA,
           lwr = estimate - Z90 * se,
           upr = estimate + Z90 * se,
           true_positive  = (true_effect != 0) & significant & (sign(estimate) == sign(true_effect)),
           true_negative  = (true_effect == 0) & !significant,
           false_positive = (true_effect == 0) & significant,                      # (none here if all nonzero)
           false_negative = (true_effect != 0) & !significant,
           method = 'ancombc')
  
  prev_df <- tibble(taxon = colnames(rel_abundances),
                    prevalence = colMeans(rel_abundances > 0))
  
  # ---- ORM ----
  res_orm2 <- res_orm |>
    left_join(true_effects, by = "taxon") |>
    mutate(q = p.adjust(p_value, "fdr"),
           lwr = estimate - Z90 * se,
           upr = estimate + Z90 * se,
           significant = q < ALPHA,
           true_positive  = (true_effect != 0) & significant & (sign(estimate) == sign(true_effect)),
           true_negative  = (true_effect == 0) & !significant,
           false_positive = (true_effect == 0) & significant,
           false_negative = (true_effect != 0) & !significant,
           method = "orm") |>
    left_join(prev_df, by = "taxon")
  
  # ---- LinDA  ----
  res_linda <- obj_linda$output$groupcase |>
    rownames_to_column("taxon") |>
    select(taxon, estimate = log2FoldChange, se = lfcSE, p = pvalue, q = padj) |>
    left_join(true_effects, by="taxon") |>
    mutate(lwr = estimate - Z90 * se,
           upr = estimate + Z90 * se,
           significant = q < ALPHA,
           true_positive  = (true_effect != 0) & significant & (sign(estimate) == sign(true_effect)),
           true_negative  = (true_effect == 0) & !significant,
           false_positive = (true_effect == 0) & significant,
           false_negative = (true_effect != 0) & !significant,
           method = "linda") |>
    left_join(prev_df, by = "taxon")
  
  res_all <- bind_rows(res_orm2, res_ancombc, res_linda) |>
    mutate(taxon = fct_reorder(taxon, -estimate, .na_rm = TRUE))
  
  # ---------------- Ordinalization for Stan (safe cuts; K=4 like other scenarios) ----------------
  K_ord <- 4
  Xcounts <- counts2
  n_use <- nrow(Xcounts); M_use <- ncol(Xcounts); taxa_keep <- colnames(Xcounts)
  
  y_mat <- matrix(1L, nrow = n_use, ncol = M_use)
  pos_idx <- Xcounts > 0
  z_pos <- log1p(Xcounts[pos_idx])
  if (K_ord > 2) {
    qcuts_raw <- stats::quantile(z_pos, probs = c(1/3, 2/3), na.rm = TRUE, type = 8)
    qcuts <- unique(qcuts_raw)
    if (length(qcuts) < 2) qcuts <- sort(qcuts_raw + c(-1e-6, 1e-6))
  } else {
    qcuts <- numeric(0)
  }
  y_mat[pos_idx] <- 1L + as.integer(cut(log1p(Xcounts[pos_idx]),
                                        breaks = c(-Inf, qcuts, Inf),
                                        labels = FALSE, include.lowest = TRUE))
  stopifnot(min(y_mat) >= 1, max(y_mat) <= K_ord)
  
  MN <- n_use * M_use
  stan_y     <- as.integer(as.vector(y_mat))
  stan_taxon <- as.integer(rep(seq_len(M_use), each = n_use))
  stan_group <- as.integer(rep(groups, times = M_use))
  
  stan_data <- list(
    MN = MN, M = M_use, K = K_ord,
    y = stan_y, group = stan_group, taxon_idx = stan_taxon
  )
  
  # ---------------- BZIOLR  ----------------
  fit_bziolr <- mod_bziolr$sample(
    data = stan_data,
    seed = SEED_MAIN, chains = 4, parallel_chains = 4,
    threads_per_chain = 3,
    iter_warmup = 1000, iter_sampling = 1000
  )
  beta_mat_zi <- posterior::as_draws_matrix(fit_bziolr$draws("beta"))
  qs_zi <- apply(beta_mat_zi, 2, quantile, probs = c(ALPHA/2, 1 - ALPHA/2))
  res_bziolr <- tibble(
    taxon = taxa_keep,
    estimate = colMeans(beta_mat_zi),
    se = apply(beta_mat_zi, 2, sd),
    lwr = qs_zi[1,], upr = qs_zi[2,]
  ) |>
    left_join(true_effects, by = "taxon") |>
    mutate(
      significant = (lwr > 0) | (upr < 0),
      true_positive  = (true_effect != 0) & significant & (sign(estimate) == sign(true_effect)),
      true_negative  = (true_effect == 0) & !significant,
      false_positive = (true_effect == 0) & significant,
      false_negative = (true_effect != 0) & !significant,
      method = "bziolr"
    )
  
  # ---------------- BALOR  ----------------
  fit_balor <- mod_balor$sample(
    data = stan_data,
    seed = SEED_MAIN, chains = 4, parallel_chains = 4,
    threads_per_chain = 3,
    iter_warmup = 1000, iter_sampling = 1000
  )
  beta_mat_ba <- posterior::as_draws_matrix(fit_balor$draws("beta"))
  qs_ba <- apply(beta_mat_ba, 2, quantile, probs = c(ALPHA/2, 1 - ALPHA/2))
  res_balor <- tibble(
    taxon = taxa_keep,
    estimate = colMeans(beta_mat_ba),
    se = apply(beta_mat_ba, 2, sd),
    lwr = qs_ba[1,], upr = qs_ba[2,]
  ) |>
    left_join(true_effects, by = "taxon") |>
    mutate(
      significant = (lwr > 0) | (upr < 0),
      true_positive  = (true_effect != 0) & significant & (sign(estimate) == sign(true_effect)),
      true_negative  = (true_effect == 0) & !significant,
      false_positive = (true_effect == 0) & significant,
      false_negative = (true_effect != 0) & !significant,
      method = "balor"
    )
  
  # ---------------- Combine & save ----------------
  res_all2 <- bind_rows(res_all, res_bziolr, res_balor) |>
    mutate(taxon = fct_reorder(taxon, -estimate, .na_rm = TRUE))
  
  RUN_TAG <- sprintf("realistic_N%d_M%d_seed%03d", n, ncol(counts2), SEED_MAIN)
  save(
    SEED_MAIN, n, groups, K_ord,
    counts, counts2, rel_abundances, true_effects, stan_data,
    res_orm2, res_ancombc, res_linda, res_bziolr, res_balor, res_all2,
    beta_mat_zi, beta_mat_ba,
    file = paste0(RUN_TAG, ".RData")
  )
}
