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
  library(rms)
})
# ============================================================
# Runtime logger for frequentist baselines (HMP Gingival V35 subset)
# Methods: ORM, ANCOM-BC2, LinDA, DESeq2, MaAsLin2, corncob (LRT), LDM
# Output: runtime_log_frequentist.csv (one row per method run)
# ============================================================
# -----------------------------
# Robust ORM runner 
# -----------------------------
run_orm <- function(abundance, metadata, formula){
  # Build design without intercept and with explicit outcome column
  X <- model.matrix(formula, metadata) |> as.data.frame()
  if ("(Intercept)" %in% names(X)) X <- X[, setdiff(names(X), "(Intercept)"), drop = FALSE]
  df <- cbind(X, abundance = as.numeric(abundance))
  vars <- setdiff(names(df), "abundance")
  # Empty/degenerate design => return NAs
  if (length(vars) == 0) {
    return(tibble(variable = character(), estimate = numeric(), se = numeric(), p_value = numeric()))
  }
  # Full model
  fit1 <- try(rms::orm(abundance ~ ., data = df, maxiter = 100), silent = TRUE)
  out  <- tibble(variable = vars, estimate = NA_real_, se = NA_real_, p_value = NA_real_)
  if (inherits(fit1, "try-error")) return(out)
  
  # Extract slopes that actually exist in the fit (exclude thresholds)
  vc   <- try(vcov(fit1), silent = TRUE)
  if (inherits(vc, "try-error")) return(out)
  all_coef <- coef(fit1)
  slope_nms <- intersect(names(all_coef), vars)
  if (length(slope_nms)) {
    out$estimate[match(slope_nms, out$variable)] <- all_coef[slope_nms]
    out$se[match(slope_nms, out$variable)]       <- sqrt(diag(as.matrix(vc))[slope_nms])
  }
  
  # LR p-values by dropping each var
  score1 <- fit1$stats["Score"]
  for (i in seq_along(vars)) {
    var_i <- vars[i]
    df0   <- df[, setdiff(names(df), var_i), drop = FALSE]
    fit0  <- try(rms::orm(abundance ~ ., data = df0, maxiter = 100), silent = TRUE)
    if (!inherits(fit0, "try-error")) {
      score0 <- fit0$stats["Score"]
      out$p_value[i] <- as.numeric(1 - pchisq(score1 - score0, df = 1))
    }
  }
  out
}

# -----------------------------
# Load gingival SUBSET dataset
# -----------------------------
load("mb_datasets_gamboa_tuz.rds")
tse <- data_mbd_raw$HMP_2012_16S_gingival_V35_subset
tse <- tse[, tse$body_subsite %in% c("subgingival_plaque", "supragingival_plaque")]

# Extract counts & metadata
counts  <- assay(tse) |> t()
rel_abs <- as.data.frame(counts / rowSums(counts))

meta <- tse |>
  colData() |>
  as.data.frame() |>
  mutate(group = factor(body_subsite,
                        levels = c("subgingival_plaque","supragingival_plaque"),
                        labels = c("subgingival","supragingival")))

# Ground truth table (for downstream evaluation; not used for timing)
taxonomy_df <- rowData(tse) |> as.data.frame() |> tibble::rownames_to_column("taxon")
bvs <- taxonomy_df %>%
  filter(taxon_annotation %in% c("aerobic","anaerobic")) %>%
  transmute(taxon,
            ground_truth = if_else(taxon_annotation == "aerobic","supragingival","subgingival"))

# -----------------------------
# Timing helpers & metadata
# -----------------------------
rt_log_path <- "runtime_log_frequentist_hmp.csv"
dataset_id  <- "HMP_2012_gingival_V35_subset"

append_csv <- function(df, path){
  write.table(df, path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path))
}
time_it <- function(expr){
  gc(); t0 <- Sys.time(); val <- eval.parent(substitute(expr)); t1 <- Sys.time()
  list(val = val, sec = as.numeric(difftime(t1, t0, units = "secs")))
}
sanitize_to_integer_counts <- function(X){
  if (!is.matrix(X)) X <- as.matrix(X)
  mode(X) <- "numeric"; X[is.na(X)] <- 0; X[X < 0] <- 0
  storage.mode(X) <- "integer"; X
}
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  smp <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, smp)
}
log_row <- function(method, seconds, N, M, notes = NA_character_){
  tibble(dataset = dataset_id, method = method, total_sec = seconds,
         N = N, M = M, timestamp = as.character(Sys.time()), notes = notes)
}

N <- nrow(counts); M <- ncol(counts)

# -----------------------------
# TIMING RUNS (each wrapped in try to keep going on errors)
# -----------------------------

# ORM (per-taxon fits + LR p-values)
try({
  tim <- time_it({
    rel_abs %>%
      purrr::map(~ run_orm(., metadata = meta, formula = ~ group)) %>%
      bind_rows(.id = "taxon")
  })
  append_csv(log_row("ORM", tim$sec, N, M, "rms::orm per taxon; LR p via drop-1"), rt_log_path)
}, silent = TRUE)

# ANCOM-BC2
try({
  ps  <- make_physeq(counts, meta)
  tim <- time_it(ANCOMBC::ancombc(data = ps, formula = "group",
                                  p_adj_method = "BH", prv_cut = 0))
  append_csv(log_row("ANCOM-BC2", tim$sec, N, M, "phyloseq input; defaults"), rt_log_path)
}, silent = TRUE)

# LinDA
try({
  tim <- time_it(LinDA::linda(otu.tab = t(counts), meta = meta, formula = "~ group"))
  append_csv(log_row("LinDA", tim$sec, N, M, "TSS+LM (package defaults)"), rt_log_path)
}, silent = TRUE)

# DESeq2 (Wald; poscounts)
try({
  tim <- time_it({
    m   <- sanitize_to_integer_counts(t(counts))              # features x samples
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = ~ group)
    sf  <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
    DESeq2::sizeFactors(dds) <- sf
    dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
    DESeq2::results(dds, contrast = c("group","supragingival","subgingival"),
                    independentFiltering = FALSE, alpha = 0.10)
  })
  append_csv(log_row("DESeq2", tim$sec, N, M, "poscounts; Wald"), rt_log_path)
}, silent = TRUE)

# MaAsLin2 (LOG + TSS; LM)
try({
  tim <- time_it({
    Maaslin2::Maaslin2(
      input_data      = t(counts),   # features x samples
      input_metadata  = meta,
      output          = file.path(tempdir(), "maaslin2_out"),
      fixed_effects   = "group",
      normalization   = "TSS",
      transform       = "LOG",
      standardize     = FALSE,
      analysis_method = "LM",
      correction      = "BH",
      min_prevalence  = 0,
      plot_heatmap    = FALSE,
      plot_scatter    = FALSE
    )
  })
  append_csv(log_row("MaAsLin2", tim$sec, N, M, "TSS+LOG; LM; BH"), rt_log_path)
}, silent = TRUE)

# corncob (LRT; mean ~ group, phi ~ 1)
try({
  tim <- time_it({
    corncob::differentialTest(
      data              = t(counts),   # features x samples
      sample_data       = meta,
      formula           = ~ group,
      formula_null      = ~ 1,
      phi.formula       = ~ 1,
      phi.formula_null  = ~ 1,
      test              = "LRT"
    )
  })
  append_csv(log_row("corncob_LRT", tim$sec, N, M, "phi ~ 1; LRT"), rt_log_path)
}, silent = TRUE)

# LDM (compositional analysis OFF to mirror earlier call)
try({
  tim <- time_it({
    meta2 <- meta
    meta2$otu_table <- counts
    LDM::ldm(formula = otu_table ~ group, data = meta2,
             comp.anal = FALSE, verbose = FALSE, n.cores = 1)
  })
  append_csv(log_row("LDM", tim$sec, N, M, "comp.anal = FALSE; n.cores = 1"), rt_log_path)
}, silent = TRUE)

message("Runtime log written to: ", normalizePath(rt_log_path))

# ============================================================
# Runtime logger — Ravel 2011 (BV) dataset
# Methods: ORM, ANCOM-BC2, LinDA, DESeq2, MaAsLin2, corncob (LRT), LDM
# Appends to: runtime_log_frequentist.csv
# ============================================================

# ---------- helpers ----------
rt_log_path <- "runtime_log_frequentist.csv"   # <- change if expect a different file name
append_csv <- function(df, path){
  write.table(df, path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path))
}
time_it <- function(expr){
  gc(); t0 <- Sys.time(); val <- eval.parent(substitute(expr)); t1 <- Sys.time()
  list(val = val, sec = as.numeric(difftime(t1, t0, units = "secs")))
}
sanitize_to_integer_counts <- function(X){
  if (!is.matrix(X)) X <- as.matrix(X)
  mode(X) <- "numeric"; X[is.na(X)] <- 0; X[X < 0] <- 0
  storage.mode(X) <- "integer"; X
}
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)  # taxa x samples
  smp <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, smp)
}
log_row <- function(dataset_id, method, seconds, N, M, notes = NA_character_){
  tibble(dataset = dataset_id, method = method, total_sec = seconds,
         N = N, M = M, timestamp = as.character(Sys.time()), notes = notes)
}

# Robust ORM (fixed)
run_orm <- function(abundance, metadata, formula){
  X <- model.matrix(formula, metadata) |> as.data.frame()
  if ("(Intercept)" %in% names(X)) X <- X[, setdiff(names(X), "(Intercept)"), drop = FALSE]
  df <- cbind(X, abundance = as.numeric(abundance))
  vars <- setdiff(names(df), "abundance")
  if (length(vars) == 0) {
    return(tibble(variable = character(), estimate = numeric(), se = numeric(), p_value = numeric()))
  }
  fit1 <- try(rms::orm(abundance ~ ., data = df, maxiter = 100), silent = TRUE)
  out  <- tibble(variable = vars, estimate = NA_real_, se = NA_real_, p_value = NA_real_)
  if (inherits(fit1, "try-error")) return(out)
  vc <- try(vcov(fit1), silent = TRUE); if (inherits(vc, "try-error")) return(out)
  all_coef <- coef(fit1)
  slope_nms <- intersect(names(all_coef), vars)
  if (length(slope_nms)) {
    out$estimate[match(slope_nms, out$variable)] <- all_coef[slope_nms]
    out$se[match(slope_nms, out$variable)]       <- sqrt(diag(as.matrix(vc))[slope_nms])
  }
  score1 <- fit1$stats["Score"]
  for (i in seq_along(vars)) {
    var_i <- vars[i]
    df0   <- df[, setdiff(names(df), var_i), drop = FALSE]
    fit0  <- try(rms::orm(abundance ~ ., data = df0, maxiter = 100), silent = TRUE)
    if (!inherits(fit0, "try-error")) {
      score0 <- fit0$stats["Score"]
      out$p_value[i] <- as.numeric(1 - pchisq(score1 - score0, df = 1))
    }
  }
  out
}

# ---------- Load & prepare Ravel 2011 ----------
load("mb_datasets_gamboa_tuz.rds")

tse <- data_mbd_raw$Ravel_2011_16S_BV
tse <- tse[, tse$study_condition %in% c("healthy","bacterial_vaginosis")]
# prevalence filter (uses mia without attaching)
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)

counts <- assay(tse) |> t()                         # samples x taxa
meta   <- colData(tse) |> as.data.frame() |>
  mutate(group = factor(study_condition,
                        levels = c("healthy","bacterial_vaginosis"),
                        labels = c("healthy","bv")))

# align samples between counts and meta
common <- intersect(rownames(counts), rownames(meta))
meta   <- meta[common, , drop = FALSE]
counts <- counts[common, , drop = FALSE]

# define identifiers used in logging (these were missing before)
dataset_id <- "Ravel_2011_16S_BV"
N <- nrow(counts); M <- ncol(counts)

message(sprintf("Dataset: %s | N=%d, M=%d | groups: %s", dataset_id, N, M, paste(levels(meta$group), collapse=" vs ")))

# ---------- TIMING RUNS ----------
# ORM
try({
  rel_abs <- as.data.frame(counts / rowSums(counts))
  tim <- time_it({
    rel_abs %>% purrr::map(~ run_orm(., metadata = meta, formula = ~ group)) %>%
      bind_rows(.id = "taxon")
  })
  append_csv(log_row(dataset_id, "ORM", tim$sec, N, M, "rms::orm per taxon; drop-1 LR p"), rt_log_path)
  message("Logged ORM runtime")
}, silent = TRUE)

# ANCOM-BC2
try({
  ps  <- make_physeq(counts, meta)
  tim <- time_it(ANCOMBC::ancombc(data = ps, formula = "group",
                                  p_adj_method = "BH", prv_cut = 0))
  append_csv(log_row(dataset_id, "ANCOM-BC2", tim$sec, N, M, "phyloseq input; defaults"), rt_log_path)
  message("Logged ANCOM-BC2 runtime")
}, silent = TRUE)

# LinDA
try({
  tim <- time_it(LinDA::linda(otu.tab = t(counts), meta = meta, formula = "~ group"))
  append_csv(log_row(dataset_id, "LinDA", tim$sec, N, M, "TSS+LM (defaults)"), rt_log_path)
  message("Logged LinDA runtime")
}, silent = TRUE)

# DESeq2 (Wald; poscounts)  — FIXED contrast to match levels: "bv" vs "healthy"
try({
  tim <- time_it({
    m   <- sanitize_to_integer_counts(t(counts))              # features x samples
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = ~ group)
    sf  <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
    DESeq2::sizeFactors(dds) <- sf
    dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
    DESeq2::results(dds, contrast = c("group","bv","healthy"),
                    independentFiltering = FALSE, alpha = 0.10)
  })
  append_csv(log_row(dataset_id, "DESeq2", tim$sec, N, M, "poscounts; Wald"), rt_log_path)
  message("Logged DESeq2 runtime")
}, silent = TRUE)

# MaAsLin2 (LOG + TSS; LM)
try({
  tim <- time_it({
    Maaslin2::Maaslin2(
      input_data      = t(counts),   # features x samples
      input_metadata  = meta,
      output          = file.path(tempdir(), "maaslin2_ravel_out"),
      fixed_effects   = "group",
      normalization   = "TSS",
      transform       = "LOG",
      standardize     = FALSE,
      analysis_method = "LM",
      correction      = "BH",
      min_prevalence  = 0,
      plot_heatmap    = FALSE,
      plot_scatter    = FALSE
    )
  })
  append_csv(log_row(dataset_id, "MaAsLin2", tim$sec, N, M, "TSS+LOG; LM; BH"), rt_log_path)
  message("Logged MaAsLin2 runtime")
}, silent = TRUE)

# corncob (LRT; mean ~ group, phi ~ 1)
try({
  tim <- time_it({
    corncob::differentialTest(
      data              = t(counts),   # features x samples
      sample_data       = meta,
      formula           = ~ group,
      formula_null      = ~ 1,
      phi.formula       = ~ 1,
      phi.formula_null  = ~ 1,
      test              = "LRT"
    )
  })
  append_csv(log_row(dataset_id, "corncob_LRT", tim$sec, N, M, "phi ~ 1; LRT"), rt_log_path)
  message("Logged corncob runtime")
}, silent = TRUE)

# LDM (comp.anal = FALSE)
try({
  tim <- time_it({
    meta2 <- meta
    meta2$otu_table <- counts
    LDM::ldm(formula = otu_table ~ group, data = meta2,
             comp.anal = FALSE, verbose = FALSE, n.cores = 1)
  })
  append_csv(log_row(dataset_id, "LDM", tim$sec, N, M, "comp.anal = FALSE; n.cores = 1"), rt_log_path)
  message("Logged LDM runtime")
}, silent = TRUE)

message("Runtime log updated at: ", normalizePath(rt_log_path, mustWork = FALSE))

# ============================================================
# Add dataset name to legend
# ============================================================
library(tidyverse)
library(ggsci)

rt <- read_csv("runtime_log_frequentist.csv")

# Build dataset label including name, N and M
rt <- rt %>%
  mutate(method = recode(method,
                         "ORM"      = "ORM",
                         "ANCOM-BC2"= "ANCOM-BC",
                         "LinDA"    = "LinDA",
                         "DESeq2"   = "DESeq2",
                         "MaAsLin2" = "MaAsLin2",
                         "corncob_LRT"  = "corncob",
                         "LDM"      = "LDM"),
         dataset = paste0(dataset, "\n(N=", N, "; N_taxa=", M, ")"),
         dataset = factor(dataset),
         dataset = fct_reorder(dataset, N * M))

# Compute geometric means per method
rt_gm <- rt %>%
  group_by(method) %>%
  dplyr::summarise(m = exp(mean(log(total_sec))), .groups="drop")

# Reorder methods
method_levels <- c("ORM","ANCOM-BC","LinDA","DESeq2","MaAsLin2","corncob","LDM")
rt <- rt %>% mutate(method = factor(method, levels = method_levels))

rt_gm <- rt_gm %>% mutate(method = factor(method, levels = method_levels))

# Plot
p1 <- ggplot(rt, aes(method, total_sec, color = dataset)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_point(data = rt_gm, aes(y = m), color = "black", size = 7, shape = 124) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.01", "0.1", "1", "10", "100", "1000"),
                name = "Run-time (s)") +
  scale_x_discrete(limits = rev) +
  ggsci::scale_color_jco(name = "Dataset") +
  coord_flip() +
  theme_light() +
  theme(axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.text  = element_text(size = 8),
        legend.position = "bottom")

# Save
date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0("fig_4.4_runtimeA", ".png"),
       plot = p1, width = 170, height = 150, dpi = 300, unit = "mm", bg = "white")
