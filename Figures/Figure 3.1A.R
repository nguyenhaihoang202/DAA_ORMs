# ================================================================
# Symmetric Scenario
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(rlang)
  library(patchwork)   # for layout and outer border
})

# --------------------------- Config ------------------------------
ALPHA_MAIN <- 0.1  

# Methods
method_levels <- rev(c("BALOR","BZIOLR","ORM","ANCOM-BC","LinDA",
                   "DESeq2","MaAsLin2","corncob","LDM"))

# Color palette 
method_colors <- c(
  BALOR      = "#045275",
  BZIOLR     = "#089099",
  ORM        = "#7CCBA2",
  "ANCOM-BC" = "#FCDE9C",
  LinDA      = "#F0746E",
  DESeq2     = "#7C1D6F",
  MaAsLin2   = "#ED85B0",
  corncob    = "#A3A500FF",
  LDM        = "#EEB849"
)

# ------------------------- Helpers -------------------------------
standardize_method <- function(x) {
  case_when(
    str_detect(x, regex("^balor$",   TRUE)) ~ "BALOR",
    str_detect(x, regex("^ziolr$",   TRUE)) ~ "BZIOLR",
    str_detect(x, regex("^orm$",     TRUE)) ~ "ORM",
    str_detect(x, regex("^ancombc",  TRUE)) ~ "ANCOM-BC",
    str_detect(x, regex("^linda",    TRUE)) ~ "LinDA",
    str_detect(x, regex("^deseq2",   TRUE)) ~ "DESeq2",
    str_detect(x, regex("^maaslin2", TRUE)) ~ "MaAsLin2",
    str_detect(x, regex("^corncob",  TRUE)) ~ "corncob",
    str_detect(x, regex("^ldm",      TRUE)) ~ "LDM",
    TRUE ~ x
  )
}

safe_load_env <- function(path) { e <- new.env(parent = emptyenv()); load(path, envir = e); e }

normalize_cols <- function(df, Z = qnorm(1 - ALPHA_MAIN/2)) {
  df <- as_tibble(df)
  if (!"taxon"  %in% names(df)) stop("res_merged missing 'taxon'")
  if (!"method" %in% names(df)) stop("res_merged missing 'method'")
  if (!"seed"   %in% names(df)) stop("res_merged missing 'seed'")
  
  if (!"estimate" %in% names(df)) {
    if ("est" %in% names(df))       df <- df %>% mutate(estimate = .data$est) %>% select(-est)
    else if ("beta" %in% names(df)) df$estimate <- df$beta
    else                            df$estimate <- NA_real_
  }
  if (!"se"  %in% names(df)) df$se  <- NA_real_
  if (!"p"   %in% names(df)) df$p   <- NA_real_
  if (!"q"   %in% names(df)) df$q   <- NA_real_
  if (!"lwr" %in% names(df)) df$lwr <- ifelse(!is.na(df$estimate) & !is.na(df$se),
                                              df$estimate - Z * df$se, NA_real_)
  if (!"upr" %in% names(df)) df$upr <- ifelse(!is.na(df$estimate) & !is.na(df$se),
                                              df$estimate + Z * df$se, NA_real_)
  df %>%
    mutate(
      taxon   = as.character(taxon),
      method  = standardize_method(as.character(method)),
      dir     = sign(estimate),
      seed    = as.integer(seed),
      data_id = sprintf("seed%03d", as.integer(seed))
    ) %>%
    relocate(seed, data_id, taxon, method, estimate, se, lwr, upr, p, q, dir)
}

round_up_to <- function(x, step) step * ceiling(x / step)

# ----------------------- File discovery --------------------------
# Symmetric scenario files
merged_files <- c(
  Sys.glob("sym_FORCED_ADD_ALL_N*_M*_seed*.RData")
) %>% unique() %>% sort()

base_files <- Sys.glob("sym_FORCED_N*_M*_seed*.RData") %>% sort()
stopifnot(length(merged_files) > 0, length(base_files) > 0)

seed_of <- function(x) as.integer(str_match(basename(x), "seed(\\d{3})\\.RData$")[,2])

merged_map <- tibble(file = merged_files, seed = seed_of(merged_files)) %>% filter(!is.na(seed))
base_map   <- tibble(file = base_files,   seed = seed_of(base_files))   %>% filter(!is.na(seed))

seeds <- intersect(merged_map$seed, base_map$seed) %>% sort()
stopifnot(length(seeds) > 0)
message("Seeds detected: ", paste(seeds, collapse = ", "))

# ------------------------ Load & harmonize ------------------------
res_all_list <- list()
truth_list   <- list()
bayes_draws  <- list()

for (sd in seeds) {
  e_m <- safe_load_env(merged_map$file[merged_map$seed == sd])
  e_b <- safe_load_env(base_map$file[base_map$seed   == sd])
  stopifnot(exists("res_merged",   envir = e_m))
  stopifnot(exists("true_effects", envir = e_b))
  
  # harmonized results
  df <- normalize_cols(e_m$res_merged) %>%
    filter(method %in% method_levels)
  res_all_list[[as.character(sd)]] <- df
  
  # ground truth
  te <- e_b$true_effects %>%
    as_tibble() %>%
    transmute(taxon = as.character(taxon),
              true_effect = as.numeric(true_effect)) %>%
    mutate(seed = sd, data_id = sprintf("seed%03d", sd))
  truth_list[[as.character(sd)]] <- te
  
  # posterior draws for BALOR/BZIOLR (for CI-exclusion)
  if (exists("beta_mat_zi", envir = e_b)) {
    taxa_zi <- if (exists("res_ziolr", envir = e_b)) e_b$res_ziolr$taxon else colnames(e_b$counts2)
    stopifnot(ncol(e_b$beta_mat_zi) == length(taxa_zi))
    bayes_draws[[paste0("BZIOLR_", sd)]] <- list(seed = sd, method = "BZIOLR",
                                                 taxa = taxa_zi, draws = e_b$beta_mat_zi)
  }
  if (exists("beta_mat_ba", envir = e_b)) {
    taxa_ba <- if (exists("res_balor", envir = e_b)) e_b$res_balor$taxon else colnames(e_b$counts2)
    stopifnot(ncol(e_b$beta_mat_ba) == length(taxa_ba))
    bayes_draws[[paste0("BALOR_", sd)]] <- list(seed = sd, method = "BALOR",
                                                taxa = taxa_ba, draws = e_b$beta_mat_ba)
  }
}

res_all <- bind_rows(res_all_list) %>%
  mutate(method = factor(method, levels = method_levels),
         taxon  = as.character(taxon))

truth_all <- bind_rows(truth_list) %>%
  mutate(seed = as.integer(seed), taxon = as.character(taxon))

res_all <- res_all %>%
  select(-matches("^true_effect($|\\.)")) %>%
  left_join(truth_all %>% select(seed, taxon, true_effect), by = c("seed","taxon")) %>%
  relocate(true_effect, .after = taxon)

stopifnot("true_effect" %in% names(res_all))

# Tag scenario if missing
if (!"scenario" %in% names(res_all)) {
  res_all <- res_all %>% mutate(scenario = "Symmetric")
}

# -------------------- Bayes significance grid ---------------------
bayes_sig_tbl <- purrr::map_dfr(bayes_draws, function(xx){
  est <- colMeans(xx$draws)
  dir <- sign(est)
  qs  <- apply(xx$draws, 2, quantile, probs = c(ALPHA_MAIN/2, 1 - ALPHA_MAIN/2))
  tibble(
    seed    = xx$seed,
    data_id = sprintf("seed%03d", xx$seed),
    method  = xx$method,
    taxon   = xx$taxa,
    alpha   = ALPHA_MAIN,
    sgn     = (qs[1, ] > 0) | (qs[2, ] < 0),
    dir     = dir
  )
})

# ------------------ Significance helper --------------------------
mark_signif <- function(df, alpha) {
  need <- c("seed","data_id","taxon","method","estimate","q","dir","true_effect","scenario")
  stopifnot(all(need %in% names(df)))
  
  add <- bayes_sig_tbl %>%
    filter(alpha == alpha) %>%
    select(seed, method, taxon, sgn_bayes = sgn, dir_bayes = dir)
  
  df %>%
    select(all_of(need)) %>%
    left_join(add, by = c("seed","method","taxon")) %>%
    mutate(
      sgn     = if_else(method %in% c("BZIOLR","BALOR"),
                        sgn_bayes,
                        if_else(is.na(q), NA, q < alpha, missing = NA)),
      dir_use = if_else(method %in% c("BZIOLR","BALOR"), dir_bayes, dir)
    )
}

# --------------- Per-seed metrics at α = ALPHA_MAIN --------------
SCEN <- "Symmetric"
res_sym <- res_all %>% filter(scenario == SCEN)

perf_sym <- res_sym %>%
  mark_signif(ALPHA_MAIN) %>%
  filter(!is.na(sgn)) %>%
  mutate(
    tp_dir = (!is.na(true_effect)) & (true_effect != 0) & sgn & sign(estimate) == sign(true_effect),
    fp_dir = (!is.na(true_effect)) & (
      ((true_effect == 0) & sgn) |
        ((true_effect != 0) & sgn & sign(estimate) != sign(true_effect))
    )
  ) %>%
  group_by(method, seed) %>%
  summarise(
    power = sum(tp_dir) / pmax(1, sum(true_effect != 0, na.rm = TRUE)),
    fdr   = ifelse(sum(sgn) > 0, sum(fp_dir) / sum(sgn), 0),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = method_levels))

# ---------------- Axis limits ------------------------
# FDR: ticks every 0.1
fdr_upper  <- 1
fdr_breaks <- seq(0, fdr_upper, by = 0.1)

# Power: ticks every 0.2
power_upper  <- 1
power_breaks <- seq(0, power_upper, by = 0.1)

# ---------------- Build the two panels -----
p_fdr <- ggplot(perf_sym, aes(x = fdr, y = method, fill = method, color = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.28, width = 0.6, linewidth = 0.35) +
  geom_jitter(height = 0.16, width = 0, size = 1.1, alpha = 0.85) +
  # ---- target FDR line ----
geom_vline(xintercept = ALPHA_MAIN, linetype = "dashed",
           linewidth = 0.5, color = "grey30") +
  # -------------------------
scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, fdr_upper), breaks = fdr_breaks) +
  labs(x = "FDR", y = NULL) +
  coord_cartesian(clip = "off") +  # so the label above the panel isn’t clipped
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_text(size = 9))


p_pow <- ggplot(perf_sym, aes(x = power, y = method, fill = method, color = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.28, width = 0.6, linewidth = 0.35) +
  geom_jitter(height = 0.16, width = 0, size = 1.1, alpha = 0.85) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, power_upper), breaks = power_breaks) +
  labs(x = "Power", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_text(size = 9))

# -------- Top title cell with its own box --------
p_title <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Symmetric Scenario",
           fontface = "bold", size = 4.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = NA, color = "black", linewidth = 0.6) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
  theme_void()

# ----------------- Arrange  ---------
combined <- p_title / (p_fdr | p_pow)

combined <- combined +
  plot_layout(heights = c(0.18, 1)) +  
  plot_annotation(
    theme = theme(
      plot.margin = margin(6, 6, 6, 6)
    )
  )

# ------------------------------ Save ------------------------------
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig4_symmetric_combined.png", combined, width = 10, height = 5, dpi = 300)
message("Saved: figures/fig4_symmetric_combined.png")
