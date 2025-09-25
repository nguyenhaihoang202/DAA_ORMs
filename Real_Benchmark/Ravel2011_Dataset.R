# ============================================================
# Ravel 2011 BV — ORM, ANCOM-BC, LinDA, DESeq2, MaAsLin2, corncob (fixed),
# LDM, and Bayesian ordinal models (BALOR, BZIOLR) with consistent calls
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(TreeSummarizedExperiment)
  library(phyloseq)
  library(ANCOMBC)
  library(LinDA)
  library(DESeq2)
  library(Maaslin2)
  library(corncob)
  library(LDM)
  library(cmdstanr)
  library(posterior)
})

# -----------------------------
# Load & prepare Ravel 2011 BV
# -----------------------------
load('mb_datasets_gamboa_tuz.rds')

tse <- data_mbd_raw$Ravel_2011_16S_BV
tse <- tse[, tse$study_condition %in% c('healthy','bacterial_vaginosis')]
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)

# counts: samples x taxa
counts  <- assay(tse) |> t()
rel_abs <- as.data.frame(counts / rowSums(counts))

# GROUP factor (reference = healthy)
meta <- tse |>
  colData() |> as.data.frame() |>
  mutate(group = factor(study_condition,
                        levels = c('healthy','bacterial_vaginosis'),
                        labels = c('healthy','bv')))

# Ground truth (bv-associated → "bv", hv-associated → "healthy")
bvs <- rowData(tse) |> as.data.frame() |>
  rownames_to_column('taxon') |>
  mutate(ground_truth = case_when(
    taxon_annotation == 'bv-associated' ~ 'bv',
    taxon_annotation == 'hv-associated' ~ 'healthy',
    TRUE ~ 'none'
  )) |>
  select(taxon, ground_truth)

# ============================================================
# ORM helper (per-taxon)
# ============================================================
run_orm <- function(abundance, metadata, formula){
  mm <- model.matrix(formula, metadata) |>
    cbind(abundance) |> tibble::as_tibble() |> dplyr::select(-"(Intercept)")
  inds <- 1:(ncol(mm) - 1); vars <- colnames(mm)[inds]
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  score_1 <- fit_1$stats["Score"]
  res <- data.frame(
    estimate = fit_1$coefficients[vars],
    se = sqrt(diag(vcov(fit_1))[vars]),
    p_value = NA_real_
  )
  if(length(inds) > 1){
    for(i in inds){
      fit_0 <- rms::orm(abundance ~ ., data = mm[, -i])
      score_0 <- fit_0$stats["Score"]
      res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
    }
  }else{
    res$p_value <- as.numeric(1 - pchisq(score_1, df = 1))
  }
  res |> tibble::rownames_to_column("variable")
}

# ORM across taxa and evaluation
res <- rel_abs |>
  purrr::map(~ run_orm(., metadata = meta, formula = ~ group)) |>
  bind_rows(.id = 'taxon')

sl <- 0.10  # main alpha
res2 <- res |>
  left_join(bvs, by = 'taxon') |>
  mutate(
    q = p.adjust(p_value, method = 'BH'),
    res = case_when(
      !is.na(q) & q < sl & estimate > 0 ~ 'bv',
      !is.na(q) & q < sl & estimate < 0 ~ 'healthy',
      TRUE ~ 'ns'
    ),
    correct   = res == ground_truth & ground_truth %in% c("bv","healthy"),
    incorrect = res != "ns" & ground_truth %in% c("bv","healthy") & res != ground_truth
  )

eval_orm <- res2 |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(
    TP = sum(correct),
    FP = sum(incorrect),
    FN = sum(res == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "orm")

# ============================================================
# ANCOM-BC
# ============================================================
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)  # taxa x samples
  smp <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, smp)
}
ps <- make_physeq(counts, meta)

obj_ancombc <- ANCOMBC::ancombc(
  data = ps, formula = "group", p_adj_method = "BH", prv_cut = 0
)

pull_group_col <- function(df){
  cand <- intersect(colnames(df), c("groupbv","group_bv"))
  if (length(cand) == 0) cand <- colnames(df)[grep("^group", colnames(df))][1]
  cand
}
col_lfc <- pull_group_col(obj_ancombc$res$lfc)

res_ancombc <- tibble(
  taxon    = obj_ancombc$res$lfc$taxon,
  estimate = obj_ancombc$res$lfc[[col_lfc]],
  se       = obj_ancombc$res$se[[col_lfc]],
  p        = obj_ancombc$res$p_val[[col_lfc]],
  q        = obj_ancombc$res$q_val[[col_lfc]]
) |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      q < sl & estimate > 0 ~ "bv",
      q < sl & estimate < 0 ~ "healthy",
      TRUE ~ "ns"
    ),
    correct   = call == ground_truth & ground_truth %in% c("bv","healthy"),
    incorrect = call != "ns" & ground_truth %in% c("bv","healthy") & call != ground_truth
  )

eval_ancombc <- res_ancombc |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "ancombc")

# ============================================================
# LinDA
# ============================================================
obj_linda <- LinDA::linda(otu.tab = t(counts), meta = meta, formula = "~ group")
coef_name <- if ("groupbv" %in% names(obj_linda$output)) "groupbv" else {
  names(obj_linda$output)[grep("^group", names(obj_linda$output))][1]
}

res_linda <- obj_linda$output[[coef_name]] |>
  rownames_to_column("taxon") |>
  transmute(taxon, estimate = log2FoldChange, se = lfcSE, p = pvalue, q = padj) |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      q < sl & estimate > 0 ~ "bv",
      q < sl & estimate < 0 ~ "healthy",
      TRUE ~ "ns"
    ),
    correct   = call == ground_truth & ground_truth %in% c("bv","healthy"),
    incorrect = call != "ns" & ground_truth %in% c("bv","healthy") & call != ground_truth
  )

eval_linda <- res_linda |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "linda")

# ============================================================
# BALOR (Bayesian ordinal; K=4 on relative abundance with zero bin)
# ============================================================
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x)
  z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L
  nx <- x[!z]; if (length(nx) == 0) return(out)
  qs <- quantile(nx, probs = seq(0, 1, length.out = K_nonzero + 1), na.rm = TRUE, type = 8)
  eps <- 1e-9 * (max(nx, na.rm=TRUE) - min(nx, na.rm=TRUE) + 1)
  for (j in 2:length(qs)) if (qs[j] <= qs[j-1]) qs[j] <- qs[j-1] + eps
  out[!z] <- cut(nx, breaks = qs, include.lowest = TRUE, labels = FALSE) + 1L
  out
}

K_nonzero <- 3
rel_cats  <- apply(rel_abs, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)
y_mat     <- t(rel_cats)                                      # taxa x samples
M <- nrow(y_mat); N <- ncol(y_mat); K_ord <- max(y_mat)
taxa_names <- colnames(rel_abs)

# Direction: +beta means higher in BV
group_num  <- ifelse(meta$group == "bv", 1L, 0L)

stan_data_balor <- list(
  MN = M * N, M = M, K = K_ord,
  y = as.vector(y_mat),
  group = rep(group_num, each = M),
  taxon_idx = rep(seq_len(M), times = N)
)

mod_balor <- cmdstan_model("BALOR.stan", cpp_options = list(stan_threads = TRUE))
fit_balor <- mod_balor$sample(
  data = stan_data_balor,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1
)

res_balor <- fit_balor$draws("beta") |> as_draws_df() |>
  dplyr::select(dplyr::starts_with("beta[")) |>
  tidyr::pivot_longer(everything(), names_to = "taxon_ix", values_to = "beta") |>
  group_by(taxon_ix) |>
  summarise(
    est = median(beta),
    lwr95 = quantile(beta, 0.025),
    upr95 = quantile(beta, 0.975),
    p_gt0 = mean(beta > 0),
    p_lt0 = mean(beta < 0),
    .groups = "drop"
  ) |>
  mutate(ix = as.integer(gsub("beta\\[|\\]", "", taxon_ix)),
         taxon = taxa_names[ix]) |>
  select(taxon, est, lwr95, upr95, p_gt0, p_lt0)

PP_THRESH <- 1 - sl/2
res_balor_calls <- res_balor |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      p_gt0 >= PP_THRESH ~ "bv",
      p_lt0 >= PP_THRESH ~ "healthy",
      TRUE ~ "ns"
    )
  )

eval_balor <- res_balor_calls |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(TP = sum(call == ground_truth & call != "ns"),
            FP = sum(call != ground_truth & call != "ns"),
            FN = sum(call == "ns")) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "balor")

# ============================================================
# BZIOLR (zero-inflated ordinal; same coding)
# ============================================================
stan_data_bziolr <- list(
  MN = M * N, M = M, K = K_ord,
  y = as.vector(y_mat),
  group = rep(group_num, each = M),
  taxon_idx = rep(seq_len(M), times = N)
)

mod_bziolr <- cmdstan_model("BZIOLR.stan", cpp_options = list(stan_threads = TRUE))
fit_bziolr <- mod_bziolr$sample(
  data = stan_data_bziolr,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1
)

res_bziolr <- fit_bziolr$draws("beta") |> as_draws_df() |>
  dplyr::select(dplyr::starts_with("beta[")) |>
  tidyr::pivot_longer(everything(), names_to = "taxon_ix", values_to = "beta") |>
  group_by(taxon_ix) |>
  summarise(
    est = median(beta),
    lwr95 = quantile(beta, 0.025),
    upr95 = quantile(beta, 0.975),
    p_gt0 = mean(beta > 0),
    p_lt0 = mean(beta < 0),
    .groups = "drop"
  ) |>
  mutate(ix = as.integer(gsub("beta\\[|\\]", "", taxon_ix)),
         taxon = taxa_names[ix]) |>
  select(taxon, est, lwr95, upr95, p_gt0, p_lt0)

res_bziolr_calls <- res_bziolr |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      p_gt0 >= PP_THRESH ~ "bv",
      p_lt0 >= PP_THRESH ~ "healthy",
      TRUE ~ "ns"
    )
  )

eval_bziolr <- res_bziolr_calls |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(TP = sum(call == ground_truth & call != "ns"),
            FP = sum(call != ground_truth & call != "ns"),
            FN = sum(call == "ns")) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "bziolr")

# ============================================================
# Frequentist extras: DESeq2, MaAsLin2, corncob (FIXED), LDM
# ============================================================
ALPHA  <- sl
ZALPHA <- qnorm(1 - ALPHA/2)

sanitize_to_integer_counts <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  mode(X) <- "numeric"
  X[is.na(X)] <- 0
  X[X < 0] <- 0
  storage.mode(X) <- "integer"
  X
}

# DESeq2
run_deseq_wald_default <- function(counts, meta, fm = ~ group,
                                   contrast = c("group", "bv", "healthy")) {
  m <- t(counts) |> sanitize_to_integer_counts()
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = fm)
  sf  <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
  DESeq2::sizeFactors(dds) <- sf
  dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
  res <- DESeq2::results(dds, contrast = contrast, independentFiltering = FALSE, alpha = ALPHA)
  as.data.frame(res) |> rownames_to_column("taxon") |>
    transmute(
      taxon, estimate = log2FoldChange, se = lfcSE, p = pvalue, q = padj,
      lwr = ifelse(!is.na(estimate) & !is.na(se), estimate - ZALPHA * se, NA_real_),
      upr = ifelse(!is.na(estimate) & !is.na(se), estimate + ZALPHA * se, NA_real_),
      method = "deseq2"
    )
}

# MaAsLin2
run_maaslin2_log_tss <- function(counts, meta, fm = ~ group) {
  fit <- Maaslin2::Maaslin2(
    input_data      = t(counts),
    input_metadata  = meta,
    output          = file.path(tempdir(), "maaslin2_out"),
    fixed_effects   = all.vars(fm),
    normalization   = "TSS",
    transform       = "LOG",
    standardize     = FALSE,
    analysis_method = "LM",
    correction      = "BH",
    min_prevalence  = 0,
    plot_heatmap    = FALSE,
    plot_scatter    = FALSE
  )
  res <- fit$results
  stopifnot(all(c("feature","metadata","value","coef","pval","qval") %in% names(res)))
  res_f  <- dplyr::filter(res, metadata == "group", value == "bv")
  se_vec <- if ("stderr" %in% names(res_f)) res_f$stderr else NA_real_
  transmute(res_f,
            taxon = feature, estimate = coef, se = se_vec,
            p = pval, q = qval,
            lwr = ifelse(!is.na(se), estimate - ZALPHA * se, NA_real_),
            upr = ifelse(!is.na(se), estimate + ZALPHA * se, NA_real_),
            method = "maaslin2")
}

# corncob (LRT) — robust coef for direction + fallback sign
run_corncob <- function(counts, meta, fm = ~ group, ev = TRUE) {
  stopifnot(nrow(counts) == nrow(meta))
  if (!identical(rownames(counts), rownames(meta))) {
    common <- intersect(rownames(counts), rownames(meta))
    counts <- counts[common, , drop = FALSE]
    meta   <- meta[common, , drop = FALSE]
  }
  nvar  <- length(all.vars(fm))
  fm0   <- if (nvar > 1) update.formula(fm, ~ . - group) else ~ 1
  phi_f <- if (isTRUE(ev)) ~ 1 else ~ group
  
  obj <- corncob::differentialTest(
    formula          = fm,
    formula_null     = fm0,
    phi.formula      = phi_f,
    phi.formula_null = phi_f,
    data             = t(counts),   # features x samples
    sample_data      = meta,
    test             = "LRT"
  )
  
  taxa <- colnames(counts)
  ests <- numeric(length(taxa)); ses <- rep(NA_real_, length(taxa))
  relA <- counts / rowSums(counts)
  g_bv <- meta$group == "bv"; g_h  <- meta$group == "healthy"
  
  for (i in seq_along(taxa)) {
    mdl <- obj$all_models[[i]]
    took <- FALSE
    if (inherits(mdl, "bbdml")) {
      mu_coef <- try(corncob::coef(mdl, "mu"), silent = TRUE)
      if (!inherits(mu_coef, "try-error")) {
        idx <- which(grepl("^group", names(mu_coef)))
        if (length(idx) > 0) {
          ests[i] <- unname(mu_coef[idx[1]])  # + means higher in BV
          took <- TRUE
        }
      }
      if (took) {
        sm <- try(summary(mdl), silent = TRUE)
        if (!inherits(sm, "try-error") && "coefficients" %in% names(sm)) {
          cf <- as.data.frame(sm$coefficients)
          rn <- rownames(cf)
          r  <- which(grepl("^mu\\.", rn) & grepl("^group", sub("^mu\\.", "", rn)))
          if (length(r) > 0 && "Std. Error" %in% colnames(cf)) ses[i] <- cf$`Std. Error`[r[1]]
        }
      }
    }
    if (!took) {
      ests[i] <- mean(relA[g_bv, i], na.rm = TRUE) - mean(relA[g_h, i], na.rm = TRUE)
    }
  }
  
  tibble(
    taxon = taxa, estimate = ests, se = ses,
    p = as.numeric(obj$p), q = p.adjust(as.numeric(obj$p), "BH"),
    lwr = NA_real_, upr = NA_real_, method = "corncob_lrt"
  )
}

# LDM
run_ldm <- function(counts, meta, fm = ~ group, clr = FALSE){
  meta2 <- meta; meta2$otu_table <- counts
  rhs_vars <- all.vars(fm); rhs_txt <- if (length(rhs_vars)) paste(rhs_vars, collapse = " + ") else "1"
  fmla <- as.formula(paste("otu_table ~", rhs_txt))
  obj <- LDM::ldm(formula = fmla, data = meta2, comp.anal = clr, verbose = FALSE, n.cores = 1)
  stopifnot(rownames(obj$beta)[1] %in% c('groupbv','grouphealthy'))
  pvals <- as.numeric(obj$p.otu.omni)
  qvals <- as.numeric(obj$q.otu.omni)
  betas <- as.numeric(obj$beta[1, ])
  dirs  <- ifelse(rownames(obj$beta)[1] == 'grouphealthy', -1, 1) # + means BV
  ests  <- dirs * betas
  tibble(
    taxon = if (!is.null(names(qvals))) names(qvals) else colnames(obj$beta),
    estimate = ests, se = NA_real_, p = pvals, q = qvals,
    lwr = NA_real_, upr = NA_real_, method = "ldm"
  )
}

# Run extras
counts_int <- sanitize_to_integer_counts(counts)
res_deseq2  <- run_deseq_wald_default(counts_int, meta)
res_maaslin <- run_maaslin2_log_tss(counts_int, meta)
res_corncob <- run_corncob(counts_int, meta, ev = TRUE)
res_ldm     <- run_ldm(counts_int, meta, clr = FALSE)

# Attach truth + calls (same rule as others)
attach_truth_and_call <- function(df, truth_tbl, alpha = ALPHA, pos = "bv", neg = "healthy") {
  df |>
    left_join(truth_tbl, by = "taxon") |>
    mutate(
      q = ifelse(is.na(q) & !is.na(p), p.adjust(p, "BH"), q),
      call = case_when(
        q < alpha & estimate > 0 ~ pos,
        q < alpha & estimate < 0 ~ neg,
        TRUE ~ "ns"
      )
    )
}

res_deseq2  <- attach_truth_and_call(res_deseq2,  bvs)
res_maaslin <- attach_truth_and_call(res_maaslin, bvs)
res_corncob <- attach_truth_and_call(res_corncob, bvs)
res_ldm     <- attach_truth_and_call(res_ldm,     bvs)

eval_from_calls <- function(df, gt_levels = c("bv","healthy"), model_name = "") {
  df2 <- dplyr::filter(df, ground_truth %in% gt_levels)
  tibble(
    TP = sum(df2$call == df2$ground_truth & df2$call != "ns"),
    FP = sum(df2$call != df2$ground_truth & df2$call != "ns"),
    FN = sum(df2$call == "ns")
  ) |>
    mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = model_name)
}

eval_deseq2  <- eval_from_calls(res_deseq2,  model_name = "deseq2")
eval_maaslin <- eval_from_calls(res_maaslin, model_name = "maaslin2")
eval_corncob <- eval_from_calls(res_corncob, model_name = "corncob_lrt")
eval_ldm     <- eval_from_calls(res_ldm,     model_name = "ldm")

# ============================================================
# Side-by-side summaries
# ============================================================
side_by_side <- dplyr::bind_rows(
  eval_orm, eval_ancombc, eval_linda,
  eval_deseq2, eval_maaslin, eval_corncob, eval_ldm,
  eval_balor, eval_bziolr
) %>% dplyr::select(model, TP, FP, FN, precision, recall)

print(side_by_side)

# Guarded renames for any prebuilt 'calls_all'
if (exists("calls_all")) {
  nms <- names(calls_all)
  if ("bayes4_call" %in% nms) names(calls_all)[match("bayes4_call", nms)] <- "balor_call"
  if ("zi_call"    %in% nms) names(calls_all)[match("zi_call",    nms)] <- "bziolr_call"
}

# ============================================================
# Alpha-grid comparisons (0.05, 0.10, 0.15, 0.20)
# ============================================================
eval_freq_at_alpha <- function(df, alpha, positive = "bv", negative = "healthy") {
  stopifnot(all(c("taxon","ground_truth","estimate") %in% names(df)))
  q_vals <- if ("p" %in% names(df)) {
    p.adjust(df$p, method = "BH")
  } else if ("p_value" %in% names(df)) {
    p.adjust(df$p_value, method = "BH")
  } else if ("q" %in% names(df)) {
    df$q
  } else stop("Need p, p_value, or q in df")
  calls <- df %>%
    mutate(q = q_vals,
           call = case_when(
             q < alpha & estimate > 0 ~ positive,
             q < alpha & estimate < 0 ~ negative,
             TRUE ~ "ns"
           )) %>%
    filter(ground_truth %in% c("bv","healthy"))
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

eval_bayes_at_alpha <- function(df, alpha, positive = "bv", negative = "healthy") {
  stopifnot(all(c("taxon","ground_truth","p_gt0","p_lt0") %in% names(df)))
  PP_THRESH_A <- 1 - alpha/2
  calls <- df %>%
    mutate(
      call = case_when(
        p_gt0 >= PP_THRESH_A ~ positive,
        p_lt0 >= PP_THRESH_A ~ negative,
        TRUE ~ "ns"
      )
    ) %>%
    filter(ground_truth %in% c("bv","healthy"))
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

df_orm     <- res %>% left_join(bvs, by = "taxon")
df_ancombc <- res_ancombc
df_linda   <- res_linda
df_deseq2  <- res_deseq2
df_maaslin <- res_maaslin
df_corncob <- res_corncob
df_ldm     <- res_ldm
df_balor   <- res_balor %>% left_join(bvs, by = "taxon")
df_bziolr  <- res_bziolr %>% left_join(bvs, by = "taxon")

alphas <- c(0.05, 0.10, 0.15, 0.20)
accumulate_method <- function(method, eval_fun, df_arg) {
  purrr::map_dfr(alphas, ~ eval_fun(df_arg, .x)) %>% mutate(model = method, .before = 1)
}

grid_orm     <- accumulate_method("orm",         function(df,a) eval_freq_at_alpha(df_orm, a), df_orm)
grid_ancombc <- accumulate_method("ancombc",     function(df,a) eval_freq_at_alpha(df_ancombc, a), df_ancombc)
grid_linda   <- accumulate_method("linda",       function(df,a) eval_freq_at_alpha(df_linda, a), df_linda)
grid_deseq2  <- accumulate_method("deseq2",      function(df,a) eval_freq_at_alpha(df_deseq2, a), df_deseq2)
grid_maaslin <- accumulate_method("maaslin2",    function(df,a) eval_freq_at_alpha(df_maaslin, a), df_maaslin)
grid_corncob <- accumulate_method("corncob_lrt", function(df,a) eval_freq_at_alpha(df_corncob, a), df_corncob)
grid_ldm     <- accumulate_method("ldm",         function(df,a) eval_freq_at_alpha(df_ldm, a), df_ldm)
grid_balor   <- accumulate_method("balor",       function(df,a) eval_bayes_at_alpha(df_balor, a), df_balor)
grid_bziolr  <- accumulate_method("bziolr",      function(df,a) eval_bayes_at_alpha(df_bziolr, a), df_bziolr)

side_by_side_grid <- bind_rows(
  grid_orm, grid_ancombc, grid_linda,
  grid_deseq2, grid_maaslin, grid_corncob, grid_ldm,
  grid_balor, grid_bziolr
) %>% select(model, alpha, TP, FP, FN, precision, recall)

print(side_by_side_grid)
