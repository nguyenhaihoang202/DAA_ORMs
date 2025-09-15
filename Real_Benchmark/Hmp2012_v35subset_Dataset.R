# ============================================================
# Frequentist ORM: HMP_2012 Gingival V35 SUBSET (Gamboa-Tuz). Here I run 1000 iter warmup and 1000 iter sampling; take 4 hours for the ordinal model.
# Should run 500 + 500 for testing
# ============================================================
library(tidyverse)
library(TreeSummarizedExperiment)
install.packages('microbiome')
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
# ZI-Ordinal (revised) via cmdstanr  
# Uses the same ordinal y_mat / group / taxon_idx created above
# ============================================================

# Rebuild (or reuse) the data list for the ZI model
stan_data_zi <- list(
  MN = M * N,
  M  = M,
  K  = K_ord,
  y  = as.vector(y_mat),
  group = rep(group_num, each = M),
  taxon_idx = rep(seq_len(M), times = N)
)

# Compile + sample
mod_zi <- cmdstan_model("ZI_model_revised.stan",
                        cpp_options = list(stan_threads = TRUE))

fit_zi <- mod_zi$sample(
  data = stan_data_zi,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 1000, iter_sampling = 1000,   
  seed = 1
)

# Summarise per-taxon beta and call
res_zi <- fit_zi$draws("beta") |> as_draws_df() |>
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

pp_thresh_zi <- 0.975
res_zi_calls <- res_zi |>
  dplyr::mutate(
    call = case_when(
      p_gt0 >= pp_thresh_zi ~ "supragingival",
      p_lt0 >= pp_thresh_zi ~ "subgingival",
      TRUE                  ~ "ns"
    )
  ) |>
  dplyr::left_join(bvs, by = "taxon") |>
  tidyr::drop_na()

eval_zi <- res_zi_calls |>
  dplyr::summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns" & ground_truth %in% c("subgingival","supragingival"))
  ) |>
  dplyr::mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    model     = "zi_ordinal_revised"
  )

# ============================================================
# UPDATED side-by-side summary (ORM first, then Bayes)
# ============================================================
side_by_side <- dplyr::bind_rows(
  eval_orm,
  eval_bayes210825,    # Bayes4 (cmdstanr)
  eval_zi              # ZI-ordinal (cmdstanr)
)
print(side_by_side)

# ======================================================================
# Add-on block: ANCOM-BC, LinDA, rename helper, Bayes call repair, and comparison
# Paste below your existing code and run from here to the end
# ======================================================================

suppressPackageStartupMessages({
  library(phyloseq); library(ANCOMBC); library(LinDA); library(dplyr); library(tidyr); library(purrr)
})

# Safety: required objects from the upper script
stopifnot(all(c("counts","meta","bvs","res2","eval_orm") %in% ls()))

# Use your frequentist alpha if already set; default to 0.10
sl <- get0("sl", ifnotfound = 0.10)
# For Bayes, align nominal level by default
PP_THRESH <- get0("PP_THRESH", ifnotfound = 1 - sl/2)

# -----------------------------
# ANCOM-BC
# -----------------------------
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

pull_group_col <- function(df){
  cand <- intersect(colnames(df),
                    c("groupsupragingival","group_supragingival","group.supragingival","groupcase","groupbv"))
  if (length(cand) == 0) cand <- colnames(df)[grep("^group", colnames(df))][1]
  cand
}
col_lfc <- pull_group_col(obj_ancombc$res$lfc)

res_ancombc <- tibble::tibble(
  taxon    = obj_ancombc$res$lfc$taxon,
  estimate = obj_ancombc$res$lfc[[col_lfc]],
  se       = obj_ancombc$res$se[[col_lfc]],
  p        = obj_ancombc$res$p_val[[col_lfc]],
  q        = obj_ancombc$res$q_val[[col_lfc]]
) |>
  left_join(bvs, by = "taxon") |>
  mutate(
    call = case_when(
      q < sl & estimate > 0 ~ "supragingival",
      q < sl & estimate < 0 ~ "subgingival",
      TRUE ~ "ns"
    )
  )

eval_ancombc <- res_ancombc |>
  filter(ground_truth %in% c("subgingival","supragingival")) |>
  summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "ancombc")

# -----------------------------
# LinDA
# -----------------------------
obj_linda <- LinDA::linda(
  otu.tab = t(counts),  # taxa x samples
  meta    = meta,
  formula = "~ group"
)

coef_name <- if ("groupsupragingival" %in% names(obj_linda$output)) "groupsupragingival" else {
  hits <- names(obj_linda$output)[grep("^group", names(obj_linda$output))]
  hits[1]
}

res_linda <- obj_linda$output[[coef_name]] |>
  tibble::rownames_to_column("taxon") |>
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
      q < sl & estimate > 0 ~ "supragingival",
      q < sl & estimate < 0 ~ "subgingival",
      TRUE ~ "ns"
    )
  )

eval_linda <- res_linda |>
  filter(ground_truth %in% c("subgingival","supragingival")) |>
  summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "linda")

# -----------------------------
# Rename helper for Bayesian objects (no refitting)
# -----------------------------
rename_obj <- function(old, new, env = .GlobalEnv) {
  if (exists(old, envir = env, inherits = FALSE)) {
    assign(new, get(old, envir = env), envir = env)
    rm(list = old, envir = env)
    message(sprintf("Renamed %s -> %s", old, new))
  }
}

# Map any old names to the new ones
rename_obj("res_bayes4",       "res_balor")
rename_obj("res_bayes4_calls", "res_balor_calls")
rename_obj("eval_bayes4",      "eval_balor")
rename_obj("res_bayes",        "res_balor")
rename_obj("res_bayes_calls",  "res_balor_calls")
rename_obj("eval_bayes210825", "eval_balor")

rename_obj("res_zi",           "res_ziolr")
rename_obj("res_zi_calls",     "res_ziolr_calls")
rename_obj("eval_zi",          "eval_ziolr")

# -----------------------------
# Build/repair BALOR and ZIOLR calls (no refits)
# Works with either p_gt0/p_lt0 or CI columns (lwr/upr or lwr95/upr95)
# -----------------------------
make_calls_generic <- function(res_df, bvs, pos = "supragingival", neg = "subgingival",
                               PP_THRESH = NULL) {
  x <- left_join(res_df, bvs, by = "taxon")
  stopifnot("ground_truth" %in% names(x))
  if (all(c("p_gt0","p_lt0") %in% names(x))) {
    if (is.null(PP_THRESH)) stop("PP_THRESH must be provided when using p_gt0/p_lt0.")
    x <- x %>%
      mutate(
        call = case_when(
          p_gt0 >= PP_THRESH ~ pos,
          p_lt0 >= PP_THRESH ~ neg,
          TRUE ~ "ns"
        )
      )
  } else if (all(c("lwr","upr") %in% names(x))) {
    x <- x %>% mutate(call = case_when(lwr > 0 ~ pos, upr < 0 ~ neg, TRUE ~ "ns"))
  } else if (all(c("lwr95","upr95") %in% names(x))) {
    x <- x %>% mutate(call = case_when(lwr95 > 0 ~ pos, upr95 < 0 ~ neg, TRUE ~ "ns"))
  } else {
    stop("Need p_gt0/p_lt0 or CI columns (lwr/upr or lwr95/upr95) to call BALOR/ZIOLR.")
  }
  x
}

if (!exists("res_balor")) stop("res_balor not found. Load your saved BALOR results.")
if (!exists("res_ziolr")) stop("res_ziolr not found. Load your saved ZIOLR results.")

# Only rebuild 'call' if missing; otherwise keep your existing calls
if (!exists("res_balor_calls") || !"call" %in% names(res_balor_calls)) {
  res_balor_calls <- make_calls_generic(res_balor, bvs, "supragingival", "subgingival", PP_THRESH)
}
if (!exists("res_ziolr_calls") || !"call" %in% names(res_ziolr_calls)) {
  res_ziolr_calls <- make_calls_generic(res_ziolr, bvs, "supragingival", "subgingival", PP_THRESH)
}

# Evals computed directly from 'call' (no need for 'correct/incorrect' columns)
eval_balor <- res_balor_calls |>
  filter(ground_truth %in% c("subgingival","supragingival")) |>
  summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "balor")

eval_ziolr <- res_ziolr_calls |>
  filter(ground_truth %in% c("subgingival","supragingival")) |>
  summarise(
    TP = sum(call == ground_truth & call != "ns"),
    FP = sum(call != ground_truth & call != "ns"),
    FN = sum(call == "ns")
  ) |>
  mutate(precision = TP/(TP+FP), recall = TP/(TP+FN), model = "ziolr")

# ============================================================
# Comparable significance levels for HMP (α = 0.05, 0.10, 0.15, 0.20)
# No refits. Reuses res, res_ancombc, res_linda, res_balor, res_ziolr.
# ============================================================

# Frequentist evaluator (ORM, ANCOM-BC, LinDA)
eval_freq_at_alpha <- function(df, alpha, positive = "supragingival", negative = "subgingival") {
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
    filter(ground_truth %in% c("subgingival","supragingival"))
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# Bayesian evaluator (requires p_gt0/p_lt0)
eval_bayes_at_alpha <- function(df, alpha, positive = "supragingival", negative = "subgingival") {
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
    filter(ground_truth %in% c("subgingival","supragingival"))
  tibble(
    alpha = alpha,
    TP = sum(calls$call == calls$ground_truth & calls$call != "ns"),
    FP = sum(calls$call != calls$ground_truth & calls$call != "ns"),
    FN = sum(calls$call == "ns"),
    precision = TP/(TP+FP),
    recall    = TP/(TP+FN)
  )
}

# Build method-wise data frames the evaluators expect
stopifnot(exists("res"))
df_orm <- res %>% left_join(bvs, by = "taxon")
df_ancombc <- if (exists("res_ancombc")) res_ancombc else NULL
df_linda   <- if (exists("res_linda"))   res_linda   else NULL
df_balor <- res_balor %>% left_join(bvs, by = "taxon")
df_ziolr <- res_ziolr %>% left_join(bvs, by = "taxon")

alphas <- c(0.05, 0.10, 0.15, 0.20)

accumulate_method <- function(method, eval_fun, df_list) {
  purrr::map_dfr(alphas, ~ eval_fun(df_list, .x)) %>% mutate(model = method, .before = 1)
}

grid_orm      <- accumulate_method("orm",      function(df,a) eval_freq_at_alpha(df_orm, a), df_orm)
grid_ancombc  <- if (!is.null(df_ancombc)) accumulate_method("ancombc", function(df,a) eval_freq_at_alpha(df_ancombc, a), df_ancombc) else tibble()
grid_linda    <- if (!is.null(df_linda))   accumulate_method("linda",   function(df,a) eval_freq_at_alpha(df_linda, a), df_linda)   else tibble()
grid_balor    <- accumulate_method("balor",   function(df,a) eval_bayes_at_alpha(df_balor, a), df_balor)
grid_ziolr    <- accumulate_method("ziolr",   function(df,a) eval_bayes_at_alpha(df_ziolr, a), df_ziolr)

# Combined grid in Ravel-style method order
side_by_side_grid <- bind_rows(
  grid_orm, grid_ancombc, grid_linda, grid_balor, grid_ziolr
) %>% select(model, alpha, TP, FP, FN, precision, recall)

print(side_by_side_grid)

# ======================================================================
# Add-on: HMP Gingival — add DESeq2, MaAsLin2, corncob, LDM (no refits)
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
stopifnot(all(c("counts","meta","bvs","res","eval_orm") %in% ls()))
# Use your frequentist alpha if already set; default to 0.10
ALPHA  <- get0("sl", ifnotfound = 0.10)
ZALPHA <- qnorm(1 - ALPHA/2)

# ---------------- helpers ----------------
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
                                   contrast = c("group", "supragingival", "subgingival")) {
  m <- t(counts) # features x samples
  m <- sanitize_to_integer_counts(m)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = fm)
  sf  <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
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
  res_f  <- dplyr::filter(res, metadata == "group", value == "supragingival")
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
    test              = stat
  )
  
  res_list <- vector("list", length = nrow(t(counts)))
  for (i in seq_len(nrow(t(counts)))) {
    mdl <- obj$all_models[[i]]
    if (!is.logical(mdl)) {
      coefs <- as.data.frame(mdl$coefficients)
      if (nrow(coefs) >= 2) {
        res_list[[i]] <- tibble(
          estimate = coefs[2, "Estimate", drop=TRUE],
          se       = coefs[2, "Std. Error", drop=TRUE]
        )
      } else {
        res_list[[i]] <- tibble(estimate = NA_real_, se = NA_real_)
      }
    } else {
      res_list[[i]] <- tibble(estimate = NA_real_, se = NA_real_)
    }
  }
  
  coef_tab <- bind_rows(res_list) %>%
    mutate(taxon = rownames(t(counts)))
  
  out <- coef_tab %>%
    mutate(
      p = as.numeric(obj$p),
      q = p.adjust(p, "BH"),
      lwr = NA_real_, upr = NA_real_,
      method = paste0("corncob_", tolower(stat))
    )
  
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
  stopifnot(rownames(obj$beta)[1] %in% c('groupsupragingival','groupsubgingival'))
  pvals <- as.numeric(obj$p.otu.omni)
  qvals <- as.numeric(obj$q.otu.omni)
  betas <- as.numeric(obj$beta[1, ])
  dirs  <- ifelse(rownames(obj$beta)[1] == 'groupsubgingival', -1, 1) # ensure + is supragingival
  ests  <- dirs * betas
  tibble(
    taxon = if (!is.null(names(qvals))) names(qvals) else colnames(obj$beta),
    estimate = ests, se = NA_real_, p = pvals, q = qvals,
    lwr = NA_real_, upr = NA_real_,
    method = "ldm"
  )
}

# ---------------- run the extras on HMP gingival ----------------
counts_int  <- sanitize_to_integer_counts(counts)

res_deseq2  <- run_deseq_wald_default(counts_int, meta)
res_maaslin <- run_maaslin2_log_tss(counts_int, meta)
res_corncob <- run_corncob(counts_int, meta, ev = TRUE, stat = "LRT")
res_ldm     <- run_ldm(counts_int, meta, clr = FALSE)

# Attach ground truth + calls (same decision rule as your other methods)
attach_truth_and_call <- function(df, truth_tbl, alpha = ALPHA,
                                  pos = "supragingival", neg = "subgingival") {
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

eval_from_calls <- function(df, gt_levels = c("subgingival","supragingival"), model_name = "") {
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

# Extend your side-by-side table (keeps the original order first)
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
rename_obj("res_bayes",  "res_balor")
rename_obj("res_bayes4", "res_balor")
rename_obj("res_zi",     "res_ziolr")

# Safety & constants
GT     <- c("subgingival","supragingival")
alphas <- c(0.05, 0.10, 0.15, 0.20)

# Frequentist evaluator (ORM, ANCOM-BC, LinDA, + new methods)
eval_freq_at_alpha <- function(df, alpha, positive = "supragingival", negative = "subgingival") {
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
eval_bayes_at_alpha <- function(df, alpha, positive = "supragingival", negative = "subgingival") {
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

grid_orm     <- accumulate_method("orm",         function(df,a) eval_freq_at_alpha(df_orm, a), df_orm)
grid_ancombc <- if (!is.null(df_ancombc)) accumulate_method("ancombc", function(df,a) eval_freq_at_alpha(df_ancombc, a), df_ancombc) else tibble()
grid_linda   <- if (!is.null(df_linda))   accumulate_method("linda",    function(df,a) eval_freq_at_alpha(df_linda, a), df_linda)   else tibble()
grid_deseq2  <- accumulate_method("deseq2",      function(df,a) eval_freq_at_alpha(df_deseq2, a), df_deseq2)
grid_maaslin <- accumulate_method("maaslin2",    function(df,a) eval_freq_at_alpha(df_maaslin, a), df_maaslin)
grid_corncob <- accumulate_method("corncob_lrt", function(df,a) eval_freq_at_alpha(df_corncob, a), df_corncob)
grid_ldm     <- accumulate_method("ldm",         function(df,a) eval_freq_at_alpha(df_ldm, a), df_ldm)
grid_balor   <- accumulate_method("balor",       function(df,a) eval_bayes_at_alpha(df_balor, a), df_balor)
grid_ziolr   <- accumulate_method("ziolr",       function(df,a) eval_bayes_at_alpha(df_ziolr, a), df_ziolr)

side_by_side_grid <- bind_rows(
  grid_orm, grid_ancombc, grid_linda,
  grid_deseq2, grid_maaslin, grid_corncob, grid_ldm,
  grid_balor, grid_ziolr
) %>% select(model, alpha, TP, FP, FN, precision, recall)

print(side_by_side_grid)
