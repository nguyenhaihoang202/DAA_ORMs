# ============================================================
# add_methods_absLOADshift.R
# Extras (DESeq2 / MaAsLin2 / corncob / LDM) for:
# Case B: Absolute-load shift (dense positive effects)
# ============================================================

suppressPackageStartupMessages({
  library(MASS)
  library(tidyverse)
  library(DESeq2)
  library(Maaslin2)
  library(corncob)
  library(LDM)
})

# ---------------- Settings ----------------
ALPHA <- 0.10
Z90   <- qnorm(1 - ALPHA/2)        # 90% CI multiplier
SEEDS <- 25:34
n     <- 100                       # 50 control, 50 case

# ---------------- load Vandeputte data ----------------
load("mean_vars_vandeputte_250625.rds")   # columns: mean, var
load("cors_vandeputte_250625.rds")        # correlation matrix 'cors'
stopifnot(exists("mean_vars"), exists("cors"))

# ---------------- helpers ----------------
sanitize_to_integer_counts <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  rn <- rownames(X); cn <- colnames(X)
  df <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (!is.numeric(v)) v <- suppressWarnings(as.numeric(v))
    v[is.na(v)] <- 0
    v
  })
  Xn <- as.matrix(df)
  Xn[Xn < 0] <- 0
  storage.mode(Xn) <- "integer"
  if (is.null(rn)) rownames(Xn) <- paste0("s", seq_len(nrow(Xn))) else rownames(Xn) <- rn
  if (is.null(cn)) colnames(Xn) <- paste0("taxon_", seq_len(ncol(Xn))) else colnames(Xn) <- cn
  Xn
}

# ---------------- DESeq2 (Wald; Default poscounts) ----------------
run_deseq_wald_default <- function(counts, meta, fm = ~ group) {
  m <- t(counts)                 # features x samples
  mode(m) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = fm)
  sf <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
  DESeq2::sizeFactors(dds) <- sf
  dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("group", "case", "control"),
                         independentFiltering = FALSE, alpha = ALPHA)
  as.data.frame(res) %>%
    tibble::rownames_to_column("taxon") %>%
    dplyr::transmute(
      taxon,
      est = log2FoldChange,
      se  = lfcSE,
      p   = pvalue,
      q   = padj,
      lwr = ifelse(!is.na(est) & !is.na(se), est - Z90 * se, NA_real_),
      upr = ifelse(!is.na(est) & !is.na(se), est + Z90 * se, NA_real_),
      method = "DESeq2 (Wald, Default)"
    )
}

# ---------------- MaAsLin2 (LOG + TSS; LM) ----------------
run_maaslin2_log_tss <- function(counts, meta, fm = ~ group, tr = "LOG", norm = "TSS") {
  fit <- Maaslin2::Maaslin2(
    input_data      = t(counts),
    input_metadata  = meta,
    output          = file.path(tempdir(), "maaslin2_out"),
    fixed_effects   = all.vars(fm),
    normalization   = norm,
    transform       = tr,
    standardize     = FALSE,
    analysis_method = "LM",
    correction      = "BH",
    min_prevalence  = 0,
    plot_heatmap    = FALSE,
    plot_scatter    = FALSE
  )
  res <- fit$results
  need <- c("feature","metadata","value","coef","pval","qval")
  stopifnot(all(need %in% names(res)))
  
  res_f  <- dplyr::filter(res, metadata == "group", value == "case")
  se_vec <- if ("stderr" %in% names(res_f)) res_f$stderr else rep(NA_real_, nrow(res_f))
  
  dplyr::transmute(
    res_f,
    taxon = feature,
    est   = coef,
    se    = se_vec,
    p     = pval,
    q     = qval,
    lwr   = ifelse(!is.na(se_vec), est - Z90 * se_vec, NA_real_),
    upr   = ifelse(!is.na(se_vec), est + Z90 * se_vec, NA_real_),
    method = paste0("MaAsLin2 (", tr, ", ", norm, ")")
  )
}

# ---------------- corncob (single test per taxon) ----------------
run_corncob <- function(counts, meta, fm = ~ group, ev = TRUE,
                        stat = c("LRT", "Wald")) {
  stat <- match.arg(stat)
  
  nvar  <- length(all.vars(fm))
  f_phi <- ifelse(ev, '~ 1', '~ group')
  fm0   <- if (nvar > 1) update.formula(fm, ~ . - group) else ~ 1
  
  obj <- corncob::differentialTest(
    data              = t(counts),      # features x samples
    sample_data       = meta,
    formula           = fm,
    formula_null      = fm0,
    phi.formula       = formula(f_phi),
    phi.formula_null  = formula(f_phi),
    test              = "LRT"
  )
  
  res_list <- vector("list", length = nrow(t(counts)))
  for (i in seq_len(nrow(t(counts)))) {
    mdl <- obj$all_models[[i]]
    if (!is.logical(mdl)) {
      res_list[[i]] <- mdl$coefficients[2, ]
    } else {
      res_list[[i]] <- c('Estimate' = NA_real_, 'Std. Error' = NA_real_,
                         't value'  = NA_real_, 'Pr(>|t|)'  = NA_real_)
    }
  }
  
  coef_tab <- dplyr::bind_rows(res_list) %>%
    dplyr::rename(est = Estimate, se = `Std. Error`) %>%
    dplyr::mutate(taxon = rownames(t(counts)),
                  df    = Inf)
  
  if (stat == "Wald") {
    out <- coef_tab %>%
      dplyr::mutate(p = 2 * (1 - pnorm(abs(est / se))),
                    q = p.adjust(p, "BH"),
                    lwr = NA_real_, upr = NA_real_,
                    method = paste0("corncob (Wald, ev = ", ev, ")"))
  } else {
    out <- coef_tab %>%
      dplyr::mutate(p = as.numeric(obj$p),
                    q = p.adjust(p, "BH"),
                    lwr = NA_real_, upr = NA_real_,
                    method = paste0("corncob (LRT, ev = ", ev, ")"))
  }
  
  dplyr::mutate(out, data_id = meta$data_id[1]) %>%
    dplyr::select(taxon, est, se, df, p, q, lwr, upr, method, data_id)
}

# ---------------- LDM ----------------
run_ldm <- function(counts, meta, fm = ~ group, clr = FALSE){
  meta$otu_table <- counts
  rhs_vars <- all.vars(fm)
  rhs_txt  <- if (length(rhs_vars)) paste(rhs_vars, collapse = " + ") else "1"
  fmla     <- as.formula(paste("otu_table ~", rhs_txt))
  
  obj <- LDM::ldm(
    formula    = fmla,
    data       = meta,
    comp.anal  = clr,
    verbose    = FALSE,
    n.cores    = 1
  )
  
  stopifnot(rownames(obj$beta) %in% c('groupcase','groupcontrol'))
  
  pvals <- as.numeric(obj$p.otu.omni)
  qvals <- as.numeric(obj$q.otu.omni)
  betas <- as.numeric(obj$beta[1, ])
  dirs  <- ifelse(rownames(obj$beta)[1] == 'groupcontrol', -1, 1)
  ests  <- dirs * betas
  
  tibble::tibble(
    taxon   = if (!is.null(names(qvals))) names(qvals) else colnames(obj$beta),
    est     = ests,
    se      = NA_real_,
    df      = Inf,
    p       = pvals,
    q       = qvals,
    lwr     = NA_real_,
    upr     = NA_real_,
    method  = paste0('LDM (CLR = ', clr, ')'),
    data_id = meta$data_id[1]
  )
}

# ---------------- Absolute-load shift simulator ----------------
simulate_counts_absload <- function(seed, n, mean_vars, cors) {
  set.seed(seed)
  n_taxa <- nrow(mean_vars)
  groups <- rep(0:1, each = n / 2)
  
  # 1) Latent MVN without group effect
  mvn_latent <- MASS::mvrnorm(n, mu = rep(0, n_taxa), Sigma = cors)
  
  # 2) Dense directional effects (~60% nonzero, all positive)
  prop_nonzero <- 0.60
  nnz   <- max(1L, round(prop_nonzero * n_taxa))
  nz_idx <- sample.int(n_taxa, nnz, replace = FALSE)
  effects <- numeric(n_taxa)
  magnitudes <- abs(rnorm(nnz, mean = 0.35, sd = 0.15))
  effects[nz_idx] <- magnitudes
  names(effects) <- paste0("taxon_", seq_len(n_taxa))
  
  # Broadcast to cases
  effects_matrix <- (cbind(replicate(n, effects)) |> t()) * groups
  mvn_latent_plus_effect <- mvn_latent + effects_matrix
  
  # 3) Gamma marginals (Vandeputte), per-taxon parameters
  abd_means    <- mean_vars$mean
  abd_vars     <- mean_vars$var
  gamma_shapes <- (abd_means ^ 2) / abd_vars
  gamma_scales <- abd_vars / abd_means
  m_gamma_shapes <- matrix(rep(gamma_shapes, each = n), nrow = n, byrow = FALSE)
  m_gamma_scales <- matrix(rep(gamma_scales, each = n), nrow = n, byrow = FALSE)
  
  uniforms        <- pnorm(mvn_latent_plus_effect)
  true_abundances <- qgamma(uniforms, shape = m_gamma_shapes, scale = m_gamma_scales)
  
  # 4) Biases, closure, library sizes, Poisson sampling
  taxonwise_biases <- exp(rnorm(n_taxa, mean = 0, sd = 1))
  biased_abundances <- t(taxonwise_biases * t(true_abundances))
  rel_biased <- biased_abundances / rowSums(biased_abundances)
  
  library_sizes <- as.integer(round(10 ^ rnorm(n, mean = 4.0, sd = 0.20)))
  lambdas <- library_sizes * rel_biased
  counts <- rpois(n * n_taxa, lambda = as.vector(lambdas)) |>
    matrix(nrow = n, ncol = n_taxa)
  colnames(counts) <- names(effects)
  
  # Keep taxa with prevalence >= 5 samples
  rel <- counts / pmax(1L, rowSums(counts))
  keep <- colSums(rel > 0) >= 5
  counts2 <- counts[, keep, drop = FALSE]
  rel_abundances <- rel[, keep, drop = FALSE]
  
  # QC: absolute-load median ratio (case/control)
  tot_case <- rowSums(true_abundances[groups == 1, , drop = FALSE])
  tot_ctrl <- rowSums(true_abundances[groups == 0, , drop = FALSE])
  load_ratio_med <- median(tot_case / tot_ctrl)
  
  list(
    counts         = sanitize_to_integer_counts(counts),
    counts2        = sanitize_to_integer_counts(counts2),
    rel_abund      = rel_abundances,
    groups         = groups,
    true_effects   = tibble::tibble(taxon = names(effects), true_effect = effects),
    taxa_keep      = colnames(counts2),
    # QC
    effects        = effects,
    nz_idx         = nz_idx,
    load_ratio_med = load_ratio_med
  )
}

# ---------------- main loop ----------------
for (SEED_MAIN in SEEDS) {
  message(">>> Absolute-load shift EXTRAS for SEED = ", SEED_MAIN)
  sim <- simulate_counts_absload(SEED_MAIN, n, mean_vars, cors)
  
  meta <- data.frame(
    group   = factor(sim$groups, labels = c("control", "case")),
    data_id = sprintf("absLOADshift_SUP_seed%03d", SEED_MAIN),
    row.names = paste0("s", seq_len(nrow(sim$counts2)))
  )
  
  # Run methods
  res_deseq2  <- run_deseq_wald_default(sim$counts2, meta)
  res_maaslin <- run_maaslin2_log_tss(sim$counts2, meta)
  res_corncob <- run_corncob(sim$counts2, meta, ev = TRUE, stat = "LRT")
  res_ldm     <- run_ldm(sim$counts2, meta, clr = FALSE)
  
  # Align truth to kept taxa + flags
  truth_keep <- sim$true_effects %>% dplyr::filter(taxon %in% sim$taxa_keep)
  add_perf <- function(df) df %>%
    dplyr::left_join(truth_keep, by = "taxon") %>%
    dplyr::mutate(
      significant    = q < ALPHA,
      tp_dir         = (true_effect != 0) & significant & (sign(est) == sign(true_effect)),
      tn             = (true_effect == 0) & !significant,
      fp_dir         = ((true_effect == 0) & significant) |
        ((true_effect != 0) & significant & (sign(est) != sign(true_effect))),
      fn             = (true_effect != 0) & !significant
    )
  
  res_deseq2  <- add_perf(res_deseq2)
  res_maaslin <- add_perf(res_maaslin)
  res_corncob <- add_perf(res_corncob)
  res_ldm     <- add_perf(res_ldm)
  
  res_extras_only <- dplyr::bind_rows(
    res_deseq2, res_maaslin, res_corncob, res_ldm
  ) %>%
    dplyr::mutate(taxon = forcats::fct_reorder(taxon, -est, .na_rm = TRUE))
  
  # Export sim components for strict merge later
  counts         <- sim$counts
  counts2        <- sim$counts2
  rel_abundances <- sim$rel_abund
  groups         <- sim$groups
  true_effects   <- sim$true_effects
  taxa_keep      <- sim$taxa_keep
  effects        <- sim$effects
  nz_idx         <- sim$nz_idx
  load_ratio_med <- sim$load_ratio_med
  
  RUN_TAG <- sprintf("absLOADshift_EXTRAONLY_SUP_N%d_M%d_seed%03d", n, ncol(sim$counts2), SEED_MAIN)
  save(
    SEED_MAIN, n,
    counts, counts2, rel_abundances, groups, true_effects, taxa_keep,
    res_deseq2, res_maaslin, res_corncob, res_ldm,
    res_extras_only,
    # QC
    effects, nz_idx, load_ratio_med,
    file = paste0(RUN_TAG, ".RData")
  )
}
