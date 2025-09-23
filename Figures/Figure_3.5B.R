# ============================================================
# Realistic scenario — timing grid (BALOR & BZIOLR only)
# N ∈ {50, 100, 150}, M fixed at 191, 10 seeds each (25:34)
# Saves: runtime_log_realistic_scaling.csv (+ keeps Stan CSVs)
# ============================================================

library(cmdstanr)
library(tidyverse)
library(MASS)

# --------- Config (keep these fixed for fairness) ----------
Ns                  <- c(50, 100, 150)
M_fixed             <- 191
seeds               <- 25:34                          # 10 seeds
chains              <- 4
iter_warmup         <- 1000
iter_sampling       <- 1000
threads_per_chain   <- 3
K_ord               <- 4                               # ordinal bins
log_path            <- "runtime_log_realistic_scaling.csv"
out_dir             <- "cmdstan_out/realistic_scaling"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------- Load Vandeputte parameters (means/vars, correlations) ----------
load("mean_vars_vandeputte_250625.rds")  
load("cors_vandeputte_250625.rds")       

stopifnot(nrow(mean_vars) >= M_fixed, ncol(cors) >= M_fixed, nrow(cors) >= M_fixed)

# Subset to exactly M_fixed taxa to keep M constant across N
mean_vars <- mean_vars[seq_len(M_fixed), , drop = FALSE]
cors      <- cors[seq_len(M_fixed), seq_len(M_fixed), drop = FALSE]

# --------- Compile models once (exclude compile time from logs) ----------
mod_balor <- cmdstan_model("BALOR.stan",  cpp_options = list(stan_threads = TRUE))
mod_bziolr <- cmdstan_model("BZIOLR.stan",  cpp_options = list(stan_threads = TRUE))

# --------- Helpers ----------
append_csv <- function(df, path) {
  write.table(df, path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path))
}

ordinalize_counts <- function(counts, K, groups) {
  # No prevalence filtering — we keep M fixed
  X <- counts
  n_use <- nrow(X); M_use <- ncol(X)
  
  # Start with everything in bin 1 (for zeros)
  y_mat <- matrix(1L, nrow = n_use, ncol = M_use)
  pos_idx <- X > 0
  z_pos <- log1p(X[pos_idx])
  
  if (K > 2 && length(z_pos) >= K - 1) {
    # K bins total = 1 for zeros + (K-1) quantile bins for positives
    probs <- seq(1/(K-1), 1 - 1/(K-1), length.out = K - 2)
    qcuts_raw <- stats::quantile(z_pos, probs = probs, na.rm = TRUE, type = 8)
    qcuts <- unique(qcuts_raw)
    # Ensure strictly increasing cutpoints
    if (length(qcuts) < (K - 2)) {
      eps <- seq_len(K - 2) * 1e-6
      qcuts <- sort(qcuts_raw + eps)
    }
  } else if (K > 2) {
    stop("Not enough positive values to form ", K - 1, " bins.")
  } else {
    qcuts <- numeric(0)
  }
  
  if (length(z_pos) > 0) {
    y_mat[pos_idx] <- 1L + as.integer(cut(
      log1p(X[pos_idx]),
      breaks = c(-Inf, qcuts, Inf),
      labels = FALSE,
      include.lowest = TRUE
    ))
  }
  
  stopifnot(min(y_mat) >= 1, max(y_mat) <= K)
  
  MN <- n_use * M_use
  stan_y     <- as.integer(as.vector(y_mat))
  stan_taxon <- as.integer(rep(seq_len(M_use), each = n_use))
  stan_group <- as.integer(rep(groups, times = M_use))
  
  list(
    MN = MN, M = M_use, K = K,
    y = stan_y, group = stan_group, taxon_idx = stan_taxon
  )
}


simulate_realistic <- function(N, mean_vars, cors, seed) {
  set.seed(seed)
  M  <- nrow(mean_vars)
  g  <- rep(0:1, length.out = N)   # balanced groups
  # latent MVN with correlations
  mvn_latent <- MASS::mvrnorm(N, mu = rep(0, M), Sigma = cors)
  
  # centered negative effects (as in current realistic code)
  effects <- rnorm(M, mean = -0.9, sd = 0.6)
  mvn_latent_plus <- mvn_latent + tcrossprod(g, effects)
  
  # Gamma marginals (taxon-wise), replicated across samples
  abd_means <- mean_vars$mean
  abd_vars  <- mean_vars$var
  gamma_shapes <- abd_means^2 / abd_vars
  gamma_scales <- abd_vars / abd_means
  m_gamma_shapes <- matrix(rep(gamma_shapes, each = N), nrow = N, byrow = FALSE)
  m_gamma_scales <- matrix(rep(gamma_scales, each = N), nrow = N, byrow = FALSE)
  
  # Uniformization + qgamma
  mvu_latent <- pnorm(mvn_latent_plus)
  true_abund <- qgamma(mvu_latent, shape = m_gamma_shapes, scale = m_gamma_scales)
  
  # Detection biases → closure → library sizes → Poisson counts
  taxonwise_biases  <- exp(rnorm(M, mean = 0, sd = 1))
  biased_abund      <- t(taxonwise_biases * t(true_abund))
  rel_abund         <- biased_abund / rowSums(biased_abund)
  library_sizes     <- round(10 ^ rnorm(N, mean = 4.0, sd = 0.20))
  lambdas           <- library_sizes * rel_abund
  counts            <- matrix(rpois(N * M, lambda = as.vector(lambdas)), nrow = N, ncol = M)
  colnames(counts)  <- paste0("taxon_", seq_len(M))
  
  list(counts = counts, groups = g)
}

time_run <- function(fit, meta) {
  tm <- fit$time()  # $warmup, $sample, $total (seconds)
  c(warmup_sec = unname(tm$warmup), sample_sec = unname(tm$sample), total_sec = unname(tm$total), meta)
}

# --------- Main loop (N grid × seeds) ----------
for (N in Ns) {
  for (seed in seeds) {
    message(sprintf(">>> Timing realistic scenario: N=%d, M=%d, seed=%d", N, M_fixed, seed))
    
    sim   <- simulate_realistic(N, mean_vars, cors, seed)
    counts <- sim$counts
    groups <- sim$groups
    
    stan_data <- ordinalize_counts(counts, K_ord, groups)
    zero_rate <- mean(counts == 0)
    
    # Run tag + output dir for this cell
    run_tag  <- sprintf("realistic_N%d_M%d_seed%03d", N, M_fixed, seed)
    out_cell <- file.path(out_dir, run_tag)
    dir.create(out_cell, recursive = TRUE, showWarnings = FALSE)
    
    #---- bziolr ----
    fit_zi <- mod_bziolr$sample(
      data = stan_data,
      seed = seed,
      chains = chains, parallel_chains = chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      threads_per_chain = threads_per_chain,
      output_dir = out_cell
    )
    row_zi <- time_run(fit_zi, meta = c(
      model = "bziolr", N = N, M = M_fixed, seed = seed,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      threads = threads_per_chain, K = K_ord, zero_rate = zero_rate,
      tag = run_tag, timestamp = as.character(Sys.time())
    ))
    
    # ---- BALOR ----
    fit_ba <- mod_balor$sample(
      data = stan_data,
      seed = seed,
      chains = chains, parallel_chains = chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      threads_per_chain = threads_per_chain,
      output_dir = out_cell
    )
    row_ba <- time_run(fit_ba, meta = c(
      model = "BALOR", N = N, M = M_fixed, seed = seed,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      threads = threads_per_chain, K = K_ord, zero_rate = zero_rate,
      tag = run_tag, timestamp = as.character(Sys.time())
    ))
    
    # Append both rows to the CSV log (create if absent)
    df <- as.data.frame(rbind(row_zi, row_ba), stringsAsFactors = FALSE)
    append_csv(df, log_path)
  }
}

message("Done. Runtime log saved to: ", normalizePath(log_path))
message("Stan outputs saved under: ", normalizePath(out_dir))

# ============================================================
# Figure 3.5B — Runtime (Bayesian models, Realistic scenario)
# Style: points = runs (colored by N); colored tick = geo-mean HOURS per method × N
# NOTE: y-axis is linear hours (not log seconds)
# Input preference: runtime_log_realistic_scaling*.xlsx (fallback *.csv)
# Methods shown: BALOR, BZIOLR
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(readxl)
})

# ---- load the latest Bayesian timing file (xlsx preferred; csv fallback) ----
xlsx_file <- sort(c(Sys.glob("runtime_log_realistic_scaling*.xlsx"),
                    Sys.glob("runtime_log_realistic_scaling.xlsx")), decreasing = TRUE)[1]
csv_file  <- sort(c(Sys.glob("runtime_log_realistic_scaling*.csv"),
                    Sys.glob("runtime_log_realistic_scaling.csv")), decreasing = TRUE)[1]

stopifnot(length(xlsx_file) + length(csv_file) >= 1)
if (!is.na(xlsx_file) && nzchar(xlsx_file) && file.exists(xlsx_file)) {
  rt_raw <- readxl::read_excel(xlsx_file)
} else {
  rt_raw <- suppressMessages(readr::read_csv(csv_file, show_col_types = FALSE))
}

# ---- tidy & label  ----
sec_candidates <- c("elapsed_sec","total_sec","secs","seconds","elapsed_seconds","wall_time_sec")
hr_candidates  <- c("elapsed_hr","total_hr","hours","elapsed_hours","total_hours")

col_first <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

sec_col <- col_first(rt_raw, sec_candidates)
hr_col  <- col_first(rt_raw, hr_candidates)

rt0 <- rt_raw %>%
  mutate(
    .method_raw = if ("method" %in% names(.)) method else if ("model" %in% names(.)) .data$model else NA_character_,
    .seed_raw   = if ("seed"   %in% names(.)) seed   else NA,
    .N_raw      = if ("N"      %in% names(.)) N      else if ("n" %in% names(.)) .data$n else NA,
    .M_raw      = if ("M"      %in% names(.)) M      else if ("n_taxa" %in% names(.)) .data$n_taxa else NA,
    .scenario   = if ("scenario" %in% names(.)) scenario else "Realistic"
  )

# seconds -> hours (or use existing hours)
if (!is.null(hr_col)) {
  rt0 <- rt0 %>% mutate(total_hours = as.numeric(.data[[hr_col]]))
} else if (!is.null(sec_col)) {
  rt0 <- rt0 %>% mutate(total_hours = as.numeric(.data[[sec_col]]) / 3600)
} else {
  stop("No runtime column found: expected seconds or hours.")
}

# Canonical method names (BALOR / BZIOLR)
rt <- rt0 %>%
  transmute(
    method = case_when(
      str_detect(toupper(.method_raw %||% ""), "BZIOLR|ZIOLR|BZ?IOLR") ~ "BZIOLR",
      str_detect(toupper(.method_raw %||% ""), "BALOR|BAYES.*ORD|ALOR") ~ "BALOR",
      TRUE ~ as.character(.method_raw)
    ),
    total_hours = ifelse(is.finite(total_hours) & total_hours > 0, total_hours, NA_real_),
    seed = suppressWarnings(as.integer(.seed_raw)),
    N    = suppressWarnings(as.integer(.N_raw)),
    M    = suppressWarnings(as.integer(.M_raw)),
    scenario = .scenario
  ) %>%
  filter(method %in% c("BALOR","BZIOLR")) %>%
  filter(!is.na(N))

# ---- dataset label per N (include N_taxa median per N) ----
M_by_N <- rt %>%
  group_by(N) %>%
  summarise(M_lab = round(median(M, na.rm = TRUE)), .groups = "drop")

rt <- rt %>%
  mutate(dataset = paste0("N = ", N),
         dataset = factor(dataset, levels = paste0("N = ", c(50,100,150))))

# ---- geometric mean per method × N (in HOURS) ----
rt_gm <- rt %>%
  filter(!is.na(total_hours) & total_hours > 0) %>%
  group_by(method, dataset) %>%
  summarise(m_hours = exp(mean(log(total_hours))), .groups = "drop")

# ---- fixed method order ----
method_levels <- c("BALOR","BZIOLR")
rt    <- rt    %>% mutate(method = factor(method, levels = method_levels))
rt_gm <- rt_gm %>% mutate(method = factor(method, levels = method_levels))

# ---- plot (linear HOURS; color by N/dataset) ----
p1 <- ggplot(rt, aes(method, total_hours, color = dataset)) +
  geom_point(size = 1, alpha = 0.6, na.rm = TRUE) +
  # colored tick = geo-mean per method × N
  geom_point(data = rt_gm, aes(x = method, y = m_hours, color = dataset),
             inherit.aes = FALSE, size = 7, shape = 124) +
  scale_y_continuous(name = "Run-time (h)", breaks = seq(0,5,0.5)) +
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
ggsave(filename = "fig_3.5B_runtime_bayesian.png",
       plot = p1, width = 170, height = 150, dpi = 300, units = "mm", bg = "white")

message("Saved: fig_3.5B_runtime_bayesian.png")

