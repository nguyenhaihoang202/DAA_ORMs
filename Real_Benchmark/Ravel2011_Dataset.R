# ============================================================
# Ravel 2011 BV benchmark: ORM, ANCOM-BC, LinDA, Laplace-Ordinal, ZI-Ordinal
# ============================================================
library(tidyverse)
library(TreeSummarizedExperiment)

# -----------------------------
# ORM helper
# -----------------------------
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

# -----------------------------
# Load & prepare Ravel 2011 BV
# -----------------------------
load('mb_datasets_gamboa_tuz.rds')

tse <- data_mbd_raw$Ravel_2011_16S_BV
tse <- tse[, tse$study_condition %in% c('healthy','bacterial_vaginosis')]
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)

counts  <- assay(tse) |> t()                         # samples x taxa (integers)
rel_abs <- as.data.frame(counts / rowSums(counts))   # samples x taxa

meta <- tse |> colData() |> as.data.frame() |>
  mutate(group = factor(study_condition,
                        levels = c('healthy','bacterial_vaginosis'),
                        labels = c('healthy','bv')))

bvs <- rowData(tse) |> as.data.frame() |>
  rownames_to_column('taxon') |>
  mutate(ground_truth = case_when(
    taxon_annotation == 'bv-associated' ~ 'bv',
    taxon_annotation == 'hv-associated' ~ 'healthy',
    TRUE ~ 'none'
  )) |>
  select(taxon, ground_truth)

# -----------------------------
# Frequentist ORM 
# -----------------------------
res <- rel_abs |>
  purrr::map(~ run_orm(., metadata = meta, formula = ~ group)) |>
  bind_rows(.id = 'taxon')

sl <- 0.10
res2 <- res |>
  left_join(bvs, by = 'taxon') |>
  mutate(
    q = p.adjust(p_value, method = 'BH'),
    res = case_when(
      q < sl & estimate > 0 ~ 'bv',
      q < sl & estimate < 0 ~ 'healthy',
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
library(phyloseq)
library(ANCOMBC)

make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)  # taxa x samples
  smp <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, smp)
}

ps <- make_physeq(counts, meta)
obj_ancombc <- ANCOMBC::ancombc(
  data = ps, formula = "group",
  p_adj_method = "BH", prv_cut = 0
)

# Pull of group contrast column (e.g., "groupbv")
pull_group_col <- function(df){
  cand <- intersect(colnames(df), c("groupbv","group_bv"))
  if (length(cand) == 0) cand <- colnames(df)[grep("^group", colnames(df))][1]
  cand
}

col_lfc <- pull_group_col(obj_ancombc$res$lfc)
res_ancombc <- tibble(
  taxon = obj_ancombc$res$lfc$taxon,
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
  summarise(
    TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "ancombc")

# ============================================================
# LinDA 
# ============================================================
library(LinDA)

obj_linda <- LinDA::linda(
  otu.tab = t(counts),  # taxa x samples
  meta    = meta,
  formula = "~ group"
)

# Coefficient table for 'groupbv'
coef_name <- if ("groupbv" %in% names(obj_linda$output)) "groupbv" else {
  names(obj_linda$output)[grep("^group", names(obj_linda$output))][1]
}

res_linda <- obj_linda$output[[coef_name]] |>
  rownames_to_column("taxon") |>
  transmute(
    taxon,
    estimate = log2FoldChange,
    se = lfcSE,
    p  = pvalue,
    q  = padj
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

eval_linda <- res_linda |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(
    TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "linda")

# ============================================================
# Laplace-ordinal
# ============================================================
library(cmdstanr)
library(posterior)

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

# Ordinalize to K = 4 categories
K_nonzero <- 3
rel_cats  <- apply(rel_abs, 2, safe_quartile_bins_zero, K_nonzero = K_nonzero)
y_mat     <- t(rel_cats)                                      # taxa x samples
M <- nrow(y_mat); N <- ncol(y_mat); K_ord <- max(y_mat)
taxa_names <- colnames(rel_abs)
group_num  <- ifelse(meta$group == "bv", 1L, 0L)

stan_data_bayes4 <- list(
  MN = M * N, M = M, K = K_ord,
  y = as.vector(y_mat),
  group = rep(group_num, each = M),
  taxon_idx = rep(seq_len(M), times = N)
)

mod_bayes4 <- cmdstan_model("bayes_model_210825.stan", cpp_options = list(stan_threads = TRUE))
fit_bayes4 <- mod_bayes4$sample(
  data = stan_data_bayes4,
  chains = 4, parallel_chains = 4, threads_per_chain = 4,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1
)

res_bayes4 <- fit_bayes4$draws("beta") |> as_draws_df() |>
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

pp_thresh_bayes4 <- 0.975
res_bayes4_calls <- res_bayes4 |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      p_gt0 >= pp_thresh_bayes4 ~ "bv",
      p_lt0 >= pp_thresh_bayes4 ~ "healthy",
      TRUE ~ "ns"
    ),
    correct   = call == ground_truth & ground_truth %in% c("bv","healthy"),
    incorrect = call != "ns" & ground_truth %in% c("bv","healthy") & call != ground_truth
  )

eval_bayes4 <- res_bayes4_calls |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(
    TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "Laplace-ordinal")

# ============================================================
# ZI-Ordinal (revised) 
# ============================================================
stan_data_zi <- list(
  MN = M * N, M = M, K = K_ord,
  y = as.vector(y_mat),
  group = rep(group_num, each = M),
  taxon_idx = rep(seq_len(M), times = N)
)

mod_zi <- cmdstan_model("ZI_model_revised.stan", cpp_options = list(stan_threads = TRUE))
fit_zi <- mod_zi$sample(
  data = stan_data_zi,
  chains = 4, parallel_chains = 4, threads_per_chain = 4,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1
)

res_zi <- fit_zi$draws("beta") |> as_draws_df() |>
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

pp_thresh_zi <- 0.975
res_zi_calls <- res_zi |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      p_gt0 >= pp_thresh_zi ~ "bv",
      p_lt0 >= pp_thresh_zi ~ "healthy",
      TRUE ~ "ns"
    ),
    correct   = call == ground_truth & ground_truth %in% c("bv","healthy"),
    incorrect = call != "ns" & ground_truth %in% c("bv","healthy") & call != ground_truth
  )

eval_zi <- res_zi_calls |>
  filter(ground_truth %in% c("bv","healthy")) |>
  summarise(
    TP = sum(correct), FP = sum(incorrect), FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "zi-ordinal")

# ============================================================
# Final side-by-side comparison table (ordered: orm, ancombc, linda, laplace-ordinal, zi-ordinal)
# ============================================================
# Fix model names in the summary table
eval_bayes4 <- eval_bayes4 %>% dplyr::mutate(model = "balor")
eval_zi     <- eval_zi     %>% dplyr::mutate(model = "ziolr")

# Rebuild the final table in the desired order
side_by_side <- dplyr::bind_rows(
  eval_orm,
  eval_ancombc,
  eval_linda,
  eval_bayes4,
  eval_zi
) %>% dplyr::select(model, TP, FP, FN, precision, recall)

print(side_by_side)

# Fix per-taxon column names in the calls table
calls_all <- calls_all %>%
  dplyr::rename(
    balor_call = bayes4_call,
    ziolr_call = zi_call
  )

# --- helper: rename an object if it exists ---
rename_obj <- function(old, new, env = .GlobalEnv) {
  if (exists(old, envir = env, inherits = FALSE)) {
    assign(new, get(old, envir = env), envir = env)
    rm(list = old, envir = env)
    message(sprintf("Renamed %s -> %s", old, new))
  } else {
    message(sprintf("Object %s not found. Skipped.", old))
  }
}

# --- results and eval objects ---
rename_obj("res_bayes4",       "res_balor")
rename_obj("res_bayes4_calls", "res_balor_calls")
rename_obj("eval_bayes4",      "eval_balor")

rename_obj("res_zi",           "res_ziolr")
rename_obj("res_zi_calls",     "res_ziolr_calls")
rename_obj("eval_zi",          "eval_ziolr")

# --- optional: fits and data lists, if you want consistent names too ---
rename_obj("fit_bayes4",       "fit_balor")
rename_obj("stan_data_bayes4", "stan_data_balor")
rename_obj("fit_zi",           "fit_ziolr")
rename_obj("stan_data_zi",     "stan_data_ziolr")

# --- fix side_by_side model labels if it already exists ---
if (exists("side_by_side")) {
  if ("model" %in% names(side_by_side)) {
    side_by_side <- side_by_side %>%
      dplyr::mutate(model = dplyr::recode(model,
                                          "Laplace-ordinal" = "balor",
                                          "bayes4"          = "balor",
                                          "zi-ordinal"      = "ziolr",
                                          "zi"              = "ziolr"
      ))
  } else {
    # or rebuild using the renamed eval_* objects if available
    needed <- c("eval_orm","eval_ancombc","eval_linda","eval_balor","eval_ziolr")
    if (all(vapply(needed, exists, logical(1)))) {
      side_by_side <- dplyr::bind_rows(
        eval_orm, eval_ancombc, eval_linda, eval_balor, eval_ziolr
      ) %>% dplyr::select(model, TP, FP, FN, precision, recall)
    }
  }
}

# --- fix per taxon calls table if already built ---
if (exists("calls_all")) {
  nms <- names(calls_all)
  if ("bayes4_call" %in% nms) names(calls_all)[match("bayes4_call", nms)] <- "balor_call"
  if ("zi_call" %in% nms)     names(calls_all)[match("zi_call", nms)]     <- "ziolr_call"
}

# ======================================================================
# Add-on: Comparable evaluation grid for Ravel BV (α = 0.05, 0.10, 0.15, 0.20)
# No refits. Reuses res, res_ancombc, res_linda, res_bayes4/res_zi (→ balor/ziolr)
# ======================================================================

suppressPackageStartupMessages({ library(dplyr); library(purrr); library(tidyr) })

# --- helper (defined if missing) ---
if (!exists("rename_obj")) {
  rename_obj <- function(old, new, env = .GlobalEnv) {
    if (exists(old, envir = env, inherits = FALSE)) {
      assign(new, get(old, envir = env), envir = env)
      rm(list = old, envir = env)
      message(sprintf("Renamed %s -> %s", old, new))
    }
  }
}

# --- ensure BALOR / ZIOLR names (no refit) ---
rename_obj("res_bayes4", "res_balor")
rename_obj("res_zi",     "res_ziolr")

# --- safety & constants ---
stopifnot(all(c("bvs","meta","res") %in% ls()))
GT <- c("bv","healthy")            # ground-truth labels in this dataset
alphas <- c(0.05, 0.10, 0.15, 0.20)

# --- frequentist evaluator (ORM, ANCOM-BC, LinDA) ---
eval_freq_at_alpha <- function(df, alpha, positive = "bv", negative = "healthy") {
  stopifnot(all(c("taxon","ground_truth","estimate") %in% names(df)))
  # Recompute BH q-values from p if available; else use existing q
  q_vals <- if ("p" %in% names(df)) {
    p.adjust(df$p, method = "BH")
  } else if ("p_value" %in% names(df)) {
    p.adjust(df$p_value, method = "BH")
  } else if ("q" %in% names(df)) {
    df$q
  } else stop("Need p, p_value, or q in the input data frame.")
  calls <- df %>%
    mutate(q = q_vals,
           call = case_when(
             q < alpha & estimate > 0 ~ positive,
             q < alpha & estimate < 0 ~ negative,
             TRUE ~ "ns"
           )) %>%
    filter(ground_truth %in% GT)
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# --- Bayesian evaluator (BALOR, ZIOLR): PP_THRESH = 1 - alpha/2 ---
# Requires p_gt0 and p_lt0 columns (present in your res_bayes4/res_zi code)
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
    filter(ground_truth %in% GT)
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# --- assemble method-specific frames the evaluators expect ---
# ORM: rejoin to compute q per alpha
df_orm     <- res %>% dplyr::left_join(bvs, by = "taxon")
# ANCOM-BC & LinDA were already joined to bvs in your code
df_ancombc <- if (exists("res_ancombc")) res_ancombc else NULL
df_linda   <- if (exists("res_linda"))   res_linda   else NULL
# BALOR/ZIOLR: ensure present and join to bvs
if (!exists("res_balor")) stop("res_balor not found. Use the rename above or load it.")
if (!exists("res_ziolr")) stop("res_ziolr not found. Use the rename above or load it.")
df_balor   <- res_balor %>% dplyr::left_join(bvs, by = "taxon")
df_ziolr   <- res_ziolr %>% dplyr::left_join(bvs, by = "taxon")

# --- evaluate at the alpha grid (Ravel method order) ---
accumulate_method <- function(method, eval_fun, df_arg) {
  purrr::map_dfr(alphas, ~ eval_fun(df_arg, .x)) %>%
    mutate(model = method, .before = 1)
}

grid_orm     <- accumulate_method("orm",     function(df,a) eval_freq_at_alpha(df_orm, a), df_orm)
grid_ancombc <- if (!is.null(df_ancombc)) accumulate_method("ancombc", function(df,a) eval_freq_at_alpha(df_ancombc, a), df_ancombc) else tibble()
grid_linda   <- if (!is.null(df_linda))   accumulate_method("linda",   function(df,a) eval_freq_at_alpha(df_linda, a), df_linda)   else tibble()
grid_balor   <- accumulate_method("balor",  function(df,a) eval_bayes_at_alpha(df_balor, a), df_balor)
grid_ziolr   <- accumulate_method("ziolr",  function(df,a) eval_bayes_at_alpha(df_ziolr, a), df_ziolr)

side_by_side_grid <- dplyr::bind_rows(
  grid_orm, grid_ancombc, grid_linda, grid_balor, grid_ziolr
) %>% dplyr::select(model, alpha, TP, FP, FN, precision, recall)

print(side_by_side_grid)

# ======================================================================
# Add-on: Ravel 2011 — add DESeq2, MaAsLin2, corncob, LDM (no refits)
#          + comparable alpha grid with all methods
# ======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(Maaslin2)
  library(corncob)
  library(LDM)
})

# Safety: need these from the upper script
stopifnot(all(c("counts","meta","bvs","res") %in% ls()))
# Use your frequentist alpha if already set; default to 0.10
ALPHA  <- get0("sl", ifnotfound = 0.10)
ZALPHA <- qnorm(1 - ALPHA/2)

# ---------------- helpers (reused style) ----------------
sanitize_to_integer_counts <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  mode(X) <- "numeric"
  X[is.na(X)] <- 0
  X[X < 0] <- 0
  storage.mode(X) <- "integer"
  X
}

# DESeq2 (Wald; poscounts)
run_deseq_wald_default <- function(counts, meta, fm = ~ group,
                                   contrast = c("group", "bv", "healthy")) {
  m <- t(counts) # features x samples
  m <- sanitize_to_integer_counts(m)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = fm)
  sf <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
  DESeq2::sizeFactors(dds) <- sf
  dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
  res <- DESeq2::results(dds, contrast = contrast,
                         independentFiltering = FALSE, alpha = ALPHA)
  as.data.frame(res) |>
    tibble::rownames_to_column("taxon") |>
    transmute(
      taxon,
      estimate = log2FoldChange,
      se       = lfcSE,
      p        = pvalue,
      q        = padj,
      lwr      = ifelse(!is.na(estimate) & !is.na(se), estimate - ZALPHA * se, NA_real_),
      upr      = ifelse(!is.na(estimate) & !is.na(se), estimate + ZALPHA * se, NA_real_),
      method   = "deseq2"
    )
}

# MaAsLin2 (LOG + TSS; LM)
run_maaslin2_log_tss <- function(counts, meta, fm = ~ group) {
  fit <- Maaslin2::Maaslin2(
    input_data      = t(counts),  # features x samples
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
  res_f <- dplyr::filter(res, metadata == "group", value == "bv")
  se_vec <- if ("stderr" %in% names(res_f)) res_f$stderr else NA_real_
  transmute(res_f,
            taxon = feature, estimate = coef, se = se_vec,
            p = pval, q = qval,
            lwr = ifelse(!is.na(se), estimate - ZALPHA * se, NA_real_),
            upr = ifelse(!is.na(se), estimate + ZALPHA * se, NA_real_),
            method = "maaslin2")
}

# corncob (single test per taxon; LRT)
run_corncob <- function(counts, meta, fm = ~ group, ev = TRUE, stat = c("LRT","Wald")) {
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
  # pull coef/SE for group effect (row 2) when available
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
  coef_tab <- bind_rows(res_list) |>
    rename(estimate = Estimate, se = `Std. Error`) |>
    mutate(taxon = rownames(t(counts)))
  out <- coef_tab |>
    mutate(p = as.numeric(obj$p),
           q = p.adjust(p, "BH"),
           lwr = NA_real_, upr = NA_real_,
           method = "corncob_lrt")
  select(out, taxon, estimate, se, p, q, lwr, upr, method)
}

# LDM
run_ldm <- function(counts, meta, fm = ~ group, clr = FALSE){
  meta2 <- meta
  meta2$otu_table <- counts
  rhs_vars <- all.vars(fm)
  rhs_txt  <- if (length(rhs_vars)) paste(rhs_vars, collapse = " + ") else "1"
  fmla     <- as.formula(paste("otu_table ~", rhs_txt))
  obj <- LDM::ldm(formula = fmla, data = meta2, comp.anal = clr, verbose = FALSE, n.cores = 1)
  stopifnot(rownames(obj$beta)[1] %in% c('groupbv','grouphealthy'))
  pvals <- as.numeric(obj$p.otu.omni)
  qvals <- as.numeric(obj$q.otu.omni)
  betas <- as.numeric(obj$beta[1, ])
  dirs  <- ifelse(rownames(obj$beta)[1] == 'grouphealthy', -1, 1)
  ests  <- dirs * betas
  tibble(
    taxon = if (!is.null(names(qvals))) names(qvals) else colnames(obj$beta),
    estimate = ests, se = NA_real_, p = pvals, q = qvals,
    lwr = NA_real_, upr = NA_real_,
    method = "ldm"
  )
}

# ---------------- run the extras on Ravel ----------------
counts_int <- sanitize_to_integer_counts(counts)

res_deseq2  <- run_deseq_wald_default(counts_int, meta)
res_maaslin <- run_maaslin2_log_tss(counts_int, meta)
res_corncob <- run_corncob(counts_int, meta, ev = TRUE, stat = "LRT")
res_ldm     <- run_ldm(counts_int, meta, clr = FALSE)

# Attach ground truth + calls (same decision rule as your other methods)
attach_truth_and_call <- function(df, truth_tbl, alpha = ALPHA,
                                  pos = "bv", neg = "healthy") {
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

# Extend your side-by-side table (keeps your original order first)
side_by_side <- dplyr::bind_rows(
  get0("side_by_side", ifnotfound = tibble()),
  eval_deseq2, eval_maaslin, eval_corncob, eval_ldm
) %>% dplyr::select(model, TP, FP, FN, precision, recall)

print(side_by_side)

# ---------------- comparable α-grid with ALL methods ----------------
suppressPackageStartupMessages({ library(purrr) })

# (re)name Bayesian summaries if you still have the old names
rename_obj <- function(old, new, env = .GlobalEnv) {
  if (exists(old, envir = env, inherits = FALSE)) {
    assign(new, get(old, envir = env), envir = env)
    rm(list = old, envir = env)
    message(sprintf("Renamed %s -> %s", old, new))
  }
}
rename_obj("res_bayes4", "res_balor")
rename_obj("res_zi",     "res_ziolr")

# Safety & constants
GT     <- c("bv","healthy")
alphas <- c(0.05, 0.10, 0.15, 0.20)

# Frequentist evaluator (ORM, ANCOM-BC, LinDA, + new methods)
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
    filter(ground_truth %in% GT)
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# Bayesian evaluator (BALOR, ZIOLR): PP_THRESH = 1 - alpha/2
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
    filter(ground_truth %in% GT)
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# Assemble method-specific frames the evaluators expect
df_orm     <- res %>% left_join(bvs, by = "taxon")
df_ancombc <- if (exists("res_ancombc")) res_ancombc else NULL
df_linda   <- if (exists("res_linda"))   res_linda   else NULL
df_deseq2  <- res_deseq2
df_maaslin <- res_maaslin
df_corncob <- res_corncob
df_ldm     <- res_ldm

if (!exists("res_balor")) stop("res_balor not found. Use the rename above or load it.")
if (!exists("res_ziolr")) stop("res_ziolr not found. Use the rename above or load it.")
df_balor   <- res_balor %>% left_join(bvs, by = "taxon")
df_ziolr   <- res_ziolr %>% left_join(bvs, by = "taxon")

accumulate_method <- function(method, eval_fun, df_arg) {
  purrr::map_dfr(alphas, ~ eval_fun(df_arg, .x)) %>% mutate(model = method, .before = 1)
}

grid_orm     <- accumulate_method("orm",      function(df,a) eval_freq_at_alpha(df_orm, a), df_orm)
grid_ancombc <- if (!is.null(df_ancombc)) accumulate_method("ancombc", function(df,a) eval_freq_at_alpha(df_ancombc, a), df_ancombc) else tibble()
grid_linda   <- if (!is.null(df_linda))   accumulate_method("linda",   function(df,a) eval_freq_at_alpha(df_linda, a), df_linda)   else tibble()
grid_deseq2  <- accumulate_method("deseq2",   function(df,a) eval_freq_at_alpha(df_deseq2, a), df_deseq2)
grid_maaslin <- accumulate_method("maaslin2", function(df,a) eval_freq_at_alpha(df_maaslin, a), df_maaslin)
grid_corncob <- accumulate_method("corncob_lrt", function(df,a) eval_freq_at_alpha(df_corncob, a), df_corncob)
grid_ldm     <- accumulate_method("ldm",      function(df,a) eval_freq_at_alpha(df_ldm, a), df_ldm)
grid_balor   <- accumulate_method("balor",    function(df,a) eval_bayes_at_alpha(df_balor, a), df_balor)
grid_ziolr   <- accumulate_method("ziolr",    function(df,a) eval_bayes_at_alpha(df_ziolr, a), df_ziolr)

side_by_side_grid <- bind_rows(
  grid_orm, grid_ancombc, grid_linda,
  grid_deseq2, grid_maaslin, grid_corncob, grid_ldm,
  grid_balor, grid_ziolr
) %>% select(model, alpha, TP, FP, FN, precision, recall)

print(side_by_side_grid)

