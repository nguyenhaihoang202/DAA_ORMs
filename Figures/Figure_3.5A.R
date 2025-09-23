# ============================================================
# Runtime bench: Frequentist baselines, Realistic scenario (10 runs)
# Methods: ORM, ANCOM-BC, LinDA, DESeq2, MaAsLin2(LOG+TSS),
#          corncob (LRT, ev=TRUE), LDM (CLR=FALSE)
# Output: tidy CSV/RDS with per-method runtime (seconds) per replicate
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(phyloseq)
  library(ANCOMBC)
  library(LinDA)
  library(DESeq2)
  library(Maaslin2)
  library(corncob)
  library(LDM)
})

# ---------------- global settings ----------------
ALPHA <- 0.10
SEEDS <- 25:34             # 10 runs
n     <- 100               # 50 control, 50 case

# ---------------- load Vandeputte anchors ----------------
load("mean_vars_vandeputte_250625.rds")  # 'mean_vars' (mean, var)
load("cors_vandeputte_250625.rds")       # 'cors'
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

simulate_counts_once_realistic <- function(seed, n, mean_vars, cors) {
  set.seed(seed)
  n_taxa <- nrow(mean_vars)
  groups <- rep(0:1, each = n / 2)
  
  # Latent MVN + realistic effects (no exact zeros; centered negative)
  mvn_latent <- MASS::mvrnorm(n, mu = rep(0, n_taxa), Sigma = cors)
  effects <- rnorm(n_taxa, mean = -0.9, sd = 0.6)
  names(effects) <- paste0("taxon_", seq_len(n_taxa))
  
  effects_matrix <- (cbind(replicate(n, effects)) |> t()) * groups
  mvn_latent_plus_effect <- mvn_latent + effects_matrix
  
  abd_means <- mean_vars$mean
  abd_vars  <- mean_vars$var
  
  gamma_shapes <- abd_means ^ 2 / abd_vars
  m_gamma_shapes <- matrix(rep(gamma_shapes, each = n), nrow = n, byrow = FALSE)
  gamma_scales <- abd_vars / abd_means
  m_gamma_scales <- matrix(rep(gamma_scales, each = n), nrow = n, byrow = FALSE)
  
  uniforms <- pnorm(mvn_latent_plus_effect)
  true_abundances <- qgamma(uniforms, shape = m_gamma_shapes, scale = m_gamma_scales)
  
  # Detection bias → closure → Poisson counts
  n_taxa <- ncol(true_abundances)
  taxonwise_biases <- exp(rnorm(n_taxa, mean = 0, sd = 1))
  biased_abundances <- t(taxonwise_biases * t(true_abundances))
  rel_biased <- biased_abundances / rowSums(biased_abundances)
  
  library_sizes <- as.integer(round(10 ^ rnorm(n, mean = 4.0, sd = 0.20)))
  lambdas <- library_sizes * rel_biased
  counts <- rpois(n * n_taxa, lambda = lambdas) |>
    matrix(nrow = n, ncol = n_taxa)
  colnames(counts) <- names(effects)
  
  # Filters (relative>0; counts>0; prevalence >= 5)
  rel_abs <- counts / pmax(1L, rowSums(counts))
  rel_abundances <- rel_abs[, colSums(rel_abs > 0) >= 5, drop = FALSE]
  counts2 <- counts[, colSums(counts > 0) >= 5, drop = FALSE]
  
  # Ensure integer & sample rownames present
  counts  <- sanitize_to_integer_counts(counts)
  counts2 <- sanitize_to_integer_counts(counts2)
  
  list(
    counts         = counts,
    counts2        = counts2,
    rel_abund      = rel_abundances,
    groups         = groups
  )
}

# ---------------- method wrappers (timing only) ----------------

# ORM on relative abundances
run_orm <- function(rel_abund, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")))
  invisible(
    as.data.frame(rel_abund) |>
      purrr::walk(function(abundance) {
        mm <- model.matrix(~ group, meta) |>
          cbind(abundance) |>
          tibble::as_tibble() |>
          dplyr::select(-"(Intercept)")
        fit_1 <- rms::orm(abundance ~ ., data = mm)
        if (ncol(mm) > 1) {
          fit_0 <- rms::orm(abundance ~ ., data = mm[, -1, drop = FALSE])
          invisible(fit_1$stats["Score"] - fit_0$stats["Score"])
        }
      })
  )
}

# ANCOM-BC 
run_ancombc <- function(counts2, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = rownames(counts2))
  phy <- phyloseq(
    otu_table(t(counts2), taxa_are_rows = TRUE),
    sample_data(meta)
  )
  invisible(
    ANCOMBC::ancombc2(
      data = phy,
      fix_formula   = "group",
      p_adj_method  = "BH",
      prv_cut       = 0
    )
  )
}

# LinDA on counts2
run_linda <- function(counts2, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = rownames(counts2))
  invisible(
    LinDA::linda(
      otu.tab = t(counts2),
      meta    = meta,
      formula = "~ group"
    )
  )
}

# DESeq2 (Wald, poscounts) on counts2
run_deseq2 <- function(counts2, groups) {
  m <- t(counts2); mode(m) <- "integer"
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = colnames(m))   # must match sample IDs
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = meta, design = ~ group)
  sf <- DESeq2::estimateSizeFactorsForMatrix(m, type = "poscounts")
  DESeq2::sizeFactors(dds) <- sf
  invisible(DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE))
}

# MaAsLin2 (LOG + TSS) on counts2
run_maaslin2 <- function(counts2, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = rownames(counts2))  # must match columns of t(counts2)
  tmpdir <- file.path(tempdir(), paste0("maaslin2_", as.integer(runif(1,1e6,9e6))))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  invisible(
    Maaslin2::Maaslin2(
      input_data      = t(counts2),
      input_metadata  = meta,
      output          = tmpdir,
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
  )
}

# corncob (LRT, ev=TRUE) on counts2
run_corncob <- function(counts2, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = rownames(counts2))
  invisible(
    corncob::differentialTest(
      data              = t(counts2),
      sample_data       = meta,
      formula           = ~ group,
      formula_null      = ~ 1,
      phi.formula       = ~ 1,
      phi.formula_null  = ~ 1,
      test              = "LRT"
    )
  )
}

# LDM (CLR = FALSE) on counts2
run_ldm <- function(counts2, groups) {
  meta <- data.frame(group = factor(groups, labels = c("control","case")),
                     row.names = rownames(counts2))
  # LDM formula interface expects the OTU table as a matrix column in 'data'
  meta$otu_table <- counts2   # samples x taxa, rownames must be sample IDs
  invisible(
    LDM::ldm(
      formula   = otu_table ~ group,
      data      = meta,
      comp.anal = FALSE,
      verbose   = FALSE,
      n.cores   = 1
    )
  )
}

# ---------------- timing harness ----------------
time_one <- function(expr) {
  t <- system.time(force(expr))["elapsed"]
  as.numeric(t)
}

results <- list()

for (SEED_MAIN in SEEDS) {
  message(">>> Timing replicate SEED = ", SEED_MAIN)
  sim <- simulate_counts_once_realistic(SEED_MAIN, n, mean_vars, cors)
  
  # Keep sizes to annotate later
  n_samples <- nrow(sim$counts2)
  n_taxa    <- ncol(sim$counts2)
  
  # run + time with error safety
  run_safe <- function(fun, label) {
    t_sec <- NA_real_; err <- NA_character_
    gc()
    t_sec <- tryCatch(time_one(fun()),
                      error = function(e) { err <<- conditionMessage(e); NA_real_ })
    tibble(method = label, elapsed_sec = t_sec, error = err)
  }
  
  # Build a list of closures so each gets the sim data
  jobs <- list(
    run_safe(function() run_orm(sim$rel_abund, sim$groups),    "ORM"),
    run_safe(function() run_ancombc(sim$counts2, sim$groups),  "ANCOM-BC"),
    run_safe(function() run_linda(sim$counts2, sim$groups),    "LinDA"),
    run_safe(function() run_deseq2(sim$counts2, sim$groups),   "DESeq2"),
    run_safe(function() run_maaslin2(sim$counts2, sim$groups), "MaAsLin2 (LOG+TSS)"),
    run_safe(function() run_corncob(sim$counts2, sim$groups),  "corncob (LRT, ev=TRUE)"),
    run_safe(function() run_ldm(sim$counts2, sim$groups),      "LDM (CLR=FALSE)")
  )
  
  df <- bind_rows(jobs) |>
    mutate(seed = SEED_MAIN,
           n_samples = n_samples,
           n_taxa = n_taxa,
           scenario = "Realistic")
  
  results[[as.character(SEED_MAIN)]] <- df
}

rt <- bind_rows(results) |>
  relocate(scenario, seed, method, elapsed_sec, error, n_samples, n_taxa)

# ---------------- save ----------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile_csv <- paste0("runtime_frequentist_realistic_", stamp, ".csv")
outfile_rds <- paste0("runtime_frequentist_realistic_", stamp, ".rds")

readr::write_csv(rt, outfile_csv)
saveRDS(rt, outfile_rds)

message("Saved: ", outfile_csv, " and ", outfile_rds)

# ============================================================
# Figure 3.5A — Runtime (Frequentist baselines, Realistic scenario)
# Style: points = runs; black tick = geometric-mean sec per method
# Legend labels include N samples and N_taxa
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
})

# ---- load the latest timing RDS from your harness ----
rds_file <- sort(Sys.glob("runtime_frequentist_realistic_*.rds"), decreasing = TRUE)[1]
stopifnot(length(rds_file) == 1)
rt_raw <- readRDS(rds_file)

# ---- tidy & label (match earlier style) ----
rt <- rt_raw %>%
  transmute(
    method,
    total_sec = as.numeric(elapsed_sec),
    seed,
    N = as.integer(n_samples),
    M = as.integer(n_taxa),
    scenario
  ) %>%
  # recode method display names to your canon
  mutate(method = recode(method,
                         "ORM"                      = "ORM",
                         "ANCOM-BC"                = "ANCOM-BC",
                         "LinDA"                   = "LinDA",
                         "DESeq2"                  = "DESeq2",
                         "MaAsLin2 (LOG+TSS)"      = "MaAsLin2",
                         "corncob (LRT, ev=TRUE)"  = "corncob",
                         "LDM (CLR=FALSE)"         = "LDM",
                         .default = method))

# compose a single dataset label (median M across seeds avoids clutter)
N_lab <- if (length(unique(rt$N)) == 1) unique(rt$N) else round(median(rt$N, na.rm = TRUE))
M_lab <- round(median(rt$M, na.rm = TRUE))
dataset_label <- paste0("N=", N_lab)

rt <- rt %>%
  mutate(
    dataset = factor(dataset_label),
    # keep only positive times for log/geo-mean; drop NAs
    total_sec = ifelse(is.finite(total_sec) & total_sec > 0, total_sec, NA_real_)
  )

# ---- geometric mean per method (black tick) ----
rt_gm <- rt %>%
  filter(!is.na(total_sec)) %>%
  group_by(method) %>%
  summarise(m = exp(mean(log(total_sec))), .groups = "drop")

# ---- fixed method order ----
method_levels <- c("ORM","ANCOM-BC","LinDA","DESeq2","MaAsLin2","corncob","LDM")
rt    <- rt    %>% mutate(method = factor(method, levels = method_levels))
rt_gm <- rt_gm %>% mutate(method = factor(method, levels = method_levels))

# ---- plot ----
p1 <- ggplot(rt, aes(method, total_sec, color = dataset)) +
  geom_point(size = 1, alpha = 0.6, na.rm = TRUE) +
  geom_point(data = rt_gm, aes(y = m), color = "black", size = 7, shape = 124) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.01", "0.1", "1", "10", "100", "1000"),
                name = "Run-time (s)") +
  # mirror your previous flip: reverse x discrete then coord_flip
  scale_x_discrete(limits = rev) +
  ggsci::scale_color_jco(name = "Dataset") +
  coord_flip() +
  theme_light() +
  theme(axis.title.y   = element_blank(),
        axis.text.x    = element_text(size = 8),
        axis.title.x   = element_text(size = 10),
        legend.text    = element_text(size = 8),
        legend.position= "bottom",
        plot.caption   = element_text(size = 8))

# ---- save ----
ggsave(filename = "fig_3.5_runtimeA.png",
       plot = p1, width = 170, height = 150, dpi = 300, units = "mm", bg = "white")
message("Saved: fig_3.5_runtimeA.png")
