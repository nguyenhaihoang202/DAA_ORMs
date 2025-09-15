# =============================================================================
# Test script - Bayesian vs frequentist DAA on the example microbiome dataset
# Pipeline:
#   1) Load relative abundances and encode into ordinal categories with explicit zero bin
#   2) Vectorize data for Stan + prepare initial values
#   3) Fit Bayesian ordered-logit model using cmdstanr with parallel computing
#   4) Extract Bayesian taxon-level effects and sanity-check centering
#   5) Fit per-taxon frequentist ordered logit (rms::orm) as a baseline
#   6) Join and visualize effect estimates with 95 percent intervals
# =============================================================================
library(tidyverse)
library(cmdstanr)
library(posterior)
library(rms)

# ----------------------------- 1) Load and prep -------------------------------
load("example_mb_data.rds")

# Subset only taxon columns into a tibble
rel_abundances <- mb_data |> select(contains('taxon_'))

# Map relative abundances to ordinal categories with an explicit zero bin
# Output is integer in {1, ..., K} where:
#   1 = zero bin
#   2..(K_nonzero+1) = quantile bins among nonzero observations
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x)
  z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L # zeros -> category 1
  nx <- x[!z] # nonzero values
  if (length(nx) == 0) return(out)  # all zero
  qs <- quantile(nx, probs = seq(0, 1, length.out = K_nonzero + 1),
                 na.rm = TRUE, type = 8)
  # Enforce strictly increasing breaks to protect cut
  eps <- 1e-9 * (max(nx, na.rm=TRUE) - min(nx, na.rm=TRUE) + 1)
  for (j in 2:length(qs)) if (qs[j] <= qs[j-1]) qs[j] <- qs[j-1] + eps
  # Bin nonzeros into 2..(K_nonzero+1)
  out[!z] <- cut(nx, breaks = qs, include.lowest = TRUE, labels = FALSE) + 1L
  out
}
# Apply binning to each taxon column 
rel_cats <- apply(rel_abundances, 2, safe_quartile_bins_zero)

# Binary group coding 
group <- ifelse(mb_data$group == 'case', 1, 0)

# --------------------------- 2) Shape data for Stan ---------------------------
# Matrix: taxa x samples
y_mat <- t(as.matrix(rel_cats))   
M <- nrow(y_mat)
N <- ncol(y_mat)
K <- max(y_mat)

# Vectorize for Stan (obs: taxa 1..M for sample 1, then taxa 1..M for sample 2, ...)
y_vec     <- as.vector(y_mat)          # column-major: sample 1 (all taxa), sample 2, ...
group_vec <- rep(group, each = M)      # sample's group repeated for its M taxa
taxon_idx <- rep(1:M, times = N)       # 1..M per sample
MN <- M * N

stan_data <- list(
  MN = MN,
  M = M,
  K = K,
  y = y_vec,
  group = group_vec,
  taxon_idx = taxon_idx
)

# ---------------------------- 3) Stan initial values --------------------------
# Minimal init for the new model:
#  - c_taxon: per-taxon ordered cutpoints initialized to equally spaced logits
#  - tau, nu: safe starting values for asymmetric Laplace prior
stan_init <- function() {
  list(
    c_taxon = replicate(M, qlogis((1:(K-1))/K), simplify = FALSE),
    tau = 1,
    nu = 0.5
  )
}

# ------------------------------ 4) Compile and run ----------------------------
# Compile the new Stan model, enable threads for parallel reduce_sum
mod_hier <- cmdstan_model("bayes_model_210825.stan",
                          cpp_options = list(stan_threads = TRUE))

# Run the Bayesian model for this dataset
start_time_hier <- Sys.time()
fit_hier <- mod_hier$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 3,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 1,
  init = stan_init
)
end_time_hier <- Sys.time()
cat("Hierarchical model runtime:", end_time_hier - start_time_hier, "\n")

# ----------------------- 5) Extract QA variables from Stan --------------------
# beta_centered_mean is beta[m] minus the within-draw mean, for comparability across taxa
beta_cmean_mat   <- fit_hier$draws("beta_centered_mean") |> posterior::as_draws_matrix()
beta_center_mean <- fit_hier$draws("beta_center_mean")   |> as.vector()

# QA: the within-draw mean of centered betas should be near zero
summary(rowMeans(beta_cmean_mat))
hist(beta_center_mean, breaks = 40,
     main = "Per-draw mean center used in Stan (QA)",
     xlab = "beta_center_mean")

# --------------------------- 6) Summarize Bayesian betas ----------------------
# Extract posterior draws for beta[m] and compute medians, 95 percent intervals, and tail probs
draws_hier <- fit_hier$draws(variables = c("beta")) |> as_draws_df()

results_hier <- draws_hier %>%
  select(starts_with("beta[")) %>%
  rename_all(~ colnames(rel_abundances)) |> 
  pivot_longer(everything(), names_to = "taxon", values_to = "beta") %>%
  group_by(taxon) %>%
  summarise(
    est   = median(beta),
    lwr95 = quantile(beta, 0.025),
    upr95 = quantile(beta, 0.975),
    p_lt_0 = mean(beta < 0),
    p_gt_0 = mean(beta > 0),
    .groups = "drop"
  ) 
print(results_hier, n = 10)

# ------------------------- 7) Frequentist per-taxon orm -----------------------
run_orm <- function(abundance, metadata, formula) {
  # Build model matrix from metadata based on 'formula', then append response
  mm <- model.matrix(formula, metadata) |> cbind(abundance) |> as_tibble() |> select(-"(Intercept)")
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  score_1 <- fit_1$stats["Score"]
  res_freq <- data.frame(
    estimate = fit_1$coefficients[vars],
    se = sqrt(diag(vcov(fit_1))[vars]),
    p_value = NA_real_
  )
  # Score test by dropping each predictor one at a time
  if (length(inds) > 1) {
    for (i in inds) {
      fit_0 <- rms::orm(abundance ~ ., data = mm[, -i])
      score_0 <- fit_0$stats["Score"]
      res_freq$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
    }
  } else {
    res_freq$p_value <- as.numeric(1 - pchisq(score_1, df = 1))
  }
  return(res_freq |> rownames_to_column("variable"))
}

# Metadata used to build design matrix for frequentist models
meta <- data.frame(group = factor(mb_data$group, levels = c("control", "case")))
res_freq <- rel_abundances |> 
  map(~ run_orm(., metadata = mb_data, formula = ~ group)) |> 
  bind_rows(.id = 'taxon') |> 
  mutate(
    ci_lwr_95 = estimate - qnorm(.975) * se,
    ci_upr_95 = estimate + qnorm(.975) * se
  ) |> 
  as_tibble()

res_freq

# ----------------------------- 8) Join and plot -------------------------------
# Align frequentist and Bayesian summaries by taxon
res_freq_clean <- res_freq |>
  select(taxon, 
         frequentist_estimate = estimate,
         frequentist_lwr = ci_lwr_95,
         frequentist_upr = ci_upr_95)

res_hier_clean <- results_hier |>
  select(taxon,
         hierarchical_estimate = est,
         hierarchical_lwr = lwr95,
         hierarchical_upr = upr95)

comparison_table <- res_freq_clean |>
  left_join(res_hier_clean,    by = "taxon") 

# Long format with columns: model, estimate, lwr, upr
comparison_long <- comparison_table |>
  pivot_longer(
    cols = -taxon,
    names_to   = c("model", ".value"),
    names_pattern = "(frequentist|hierarchical)_(estimate|lwr|upr)"
  )

# Panel plot of estimates with 95 percent intervals for both approaches
ggplot(comparison_long, aes(x = estimate, y = taxon, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ model) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Parameter Value", y = NULL,
       title = "Estimates and 95% intervals (incl. centered Bayesian effects)",
       color = "Model") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Trimming to remove extreme frequentist intervals that dominate the axis
comparison_long_trimmed <- comparison_long |>
  mutate(ci_width = upr - lwr) |>
  # Only trim the frequentist intervals (keep all Bayesian, hierarchical intervals)
  filter(model != "frequentist" | ci_width <= 20)

# Plot with trimmed frequentist outliers for clearer scales
ggplot(comparison_long_trimmed, aes(x = estimate, y = taxon, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ model) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Parameter Value", y = NULL,
       title = "Compare Estimates and 95% Intervals by Model (Trimmed Frequentist Outliers)",
       color = "Model") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


# "Significant" taxa------------------------------------------------------------
res_freq |> 
  mutate(q = p.adjust(p_value, method = "BH")) |>
  filter(q < 0.10) |>
  arrange(p_value)

results_hier |> 
  mutate(q = 2 * pmin(p_lt_0, p_gt_0)) |> 
  filter(q < .10) |>
  arrange(q)


# =============================================================================
# 9) Use REAL counts (no pseudocount). Align taxa with rel_abundances.
# =============================================================================
# If your counts are in a separate file in wide format with taxon_* columns:
# e.g. example_mb_data_counts.rds -> object 'mb_counts' (samples x taxon_* cols)
load("example_mb_data_counts.rds")   # adjust name; object should be a data.frame or matrix
# If the object is named differently, assign it to counts_mat:
counts_mat <- as.matrix(mb_data)   # replace 'mb_counts' with your object name

# Sanity checks and alignment
stopifnot(nrow(counts_mat) == nrow(mb_data))              # same samples
common_taxa <- intersect(colnames(counts_mat), colnames(rel_abundances))
stopifnot(length(common_taxa) > 0)

# Restrict and reorder both to the same taxa
counts_mat        <- counts_mat[, common_taxa, drop = FALSE]
rel_abundances2   <- rel_abundances[, common_taxa, drop = FALSE]
taxa_names        <- common_taxa

# Integer counts required
storage.mode(counts_mat) <- "integer"

# Prevalence filter applied once, then reused downstream
keep_taxa <- colSums(counts_mat > 0) >= 5
counts2   <- counts_mat[, keep_taxa, drop = FALSE]
rel_abundances2 <- rel_abundances2[, keep_taxa, drop = FALSE]
taxa_names <- colnames(rel_abundances2)


# Build phyloseq for ANCOM-BC
make_physeq <- function(counts, meta_df){
  otu  <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  samp <- phyloseq::sample_data(meta_df)
  phyloseq::phyloseq(otu, samp)
}
meta_df <- data.frame(group = factor(mb_data$group, levels = c("control","case")))

obj_ancombc <- ANCOMBC::ancombc(
  data = make_physeq(counts2, meta_df),
  formula = "group",
  p_adj_method = "BH",
  prv_cut = 0
)

res_ancombc <- tibble(
  taxon    = rownames(obj_ancombc$res$lfc),
  estimate = obj_ancombc$res$lfc[,"groupcase"],
  se       = obj_ancombc$res$se[,"groupcase"],
  p_value  = obj_ancombc$res$p_val[,"groupcase"],
  q        = obj_ancombc$res$q_val[,"groupcase"]
) |>
  mutate(
    lwr95 = estimate - 1.96 * se,
    upr95 = estimate + 1.96 * se,
    method = "ancombc"
  )

obj_linda <- LinDA::linda(
  otu.tab = t(counts2),
  meta    = meta_df,
  formula = "~ group"
)
res_linda <- obj_linda$output$groupcase |>
  as.data.frame() |>
  rownames_to_column("taxon") |>
  transmute(
    taxon,
    estimate = log2FoldChange,   # LinDA reports log2 fold change
    se       = lfcSE,
    p_value  = pvalue,
    q        = padj,
    lwr95    = estimate - 1.96 * se,
    upr95    = estimate + 1.96 * se,
    method   = "linda"
  )

# =============================================================================
# 10) Add ZI-ordinal Bayesian model (separate Stan file). Reuse same categories.
# =============================================================================
# Use the same ordinal categories you computed earlier for Bayes4:
# rel_cats is samples x taxa for ALL taxa; subset to kept taxa for consistency
rel_cats2 <- rel_cats[, keep_taxa, drop = FALSE]
y_mat_zi  <- t(as.matrix(rel_cats2))   # taxa x samples
M_zi <- nrow(y_mat_zi); N_zi <- ncol(y_mat_zi); K_zi <- max(y_mat_zi)

# Group coding (0 or 1). If you prefer centered coding, replace with {-0.5, 0.5}
group_01 <- ifelse(mb_data$group == "case", 1L, 0L)

stan_data_zi <- list(
  MN        = M_zi * N_zi,
  M         = M_zi,
  K         = K_zi,
  y         = as.integer(as.vector(y_mat_zi)),   # vec by sample then taxon (matches above)
  taxon_idx = as.integer(rep(1:M_zi, times = N_zi)),
  group     = as.integer(rep(group_01, each = M_zi))
  # add any extra fields required by your zi-ordinal Stan file if needed
)

# Compile and fit ZI-ordinal
# Update the path if your file is named differently
zi_stan_file <- "zi_model_revised.stan"
mod_zi <- cmdstan_model(zi_stan_file, cpp_options = list(stan_threads = TRUE))

start_time_zi <- Sys.time()
fit_zi <- mod_zi$sample(
  data = stan_data_zi,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 1000, iter_sampling = 1000,
  seed = 1
)
end_time_zi <- Sys.time()
cat("ZI-ordinal runtime:", end_time_zi - start_time_zi, "\n")

beta_zi <- fit_zi$draws("beta") |> posterior::as_draws_matrix()
stopifnot(ncol(beta_zi) == M_zi)

res_zi <- tibble(
  taxon    = taxa_names,                         # same kept taxa
  estimate = colMeans(beta_zi),
  lwr95    = apply(beta_zi, 2, quantile, 0.025),
  upr95    = apply(beta_zi, 2, quantile, 0.975),
  p_lt_0   = colMeans(beta_zi < 0),
  p_gt_0   = colMeans(beta_zi > 0),
  method   = "zi_ordinal"
)

# =============================================================================
# 11) Harmonize all results and apply matched significance rules
#      Bayesian: CI excludes 0 at 90 percent for q < 0.10, and 95 percent for q < 0.05
# =============================================================================

# Existing results
res_orm <- res_freq |>
  transmute(
    taxon,
    estimate = estimate, lwr95 = ci_lwr_95, upr95 = ci_upr_95,
    p_value = p_value, q = p.adjust(p_value, method = "BH"),
    method = "orm"
  )

# Bayes4 results already computed as results_hier
res_bayes4 <- results_hier |>
  transmute(
    taxon,
    estimate = est, lwr95 = lwr95, upr95 = upr95,
    p_lt_0 = p_lt_0, p_gt_0 = p_gt_0,
    method = "bayes4"
  )

# Function to add significance flags for both families
add_sig_flags <- function(df){
  if ("q" %in% names(df)) {
    df |>
      mutate(
        sig_q10 = q < 0.10,
        sig_q05 = q < 0.05
      )
  } else {
    # Bayesian: CI exclusion and tail probability q_bayes = 2*min(p_lt_0, p_gt_0)
    stopifnot(all(c("lwr95","upr95") %in% names(df)))
    df |>
      mutate(
        # 95 percent CI already present
        sigCI95 = (lwr95 > 0) | (upr95 < 0),
        # 90 percent CI from posterior draws if available or by shrinking bounds
        # We have draws for Bayes models, so compute 90 percent properly when we build df.
        q_bayes = 2 * pmin(coalesce(p_lt_0, 0), coalesce(p_gt_0, 0)),
        sig90_rule = NA    # will be filled when we attach 90 percent intervals
      )
  }
}

# For Bayes4 and ZI we also need 90 percent CI to match q < 0.10
# Compute directly from draws
make_bayes_ci90 <- function(beta_mat, taxa, method_label){
  tibble(
    taxon = taxa,
    lwr90 = apply(beta_mat, 2, quantile, 0.05),
    upr90 = apply(beta_mat, 2, quantile, 0.95),
    method = method_label
  )
}
ci90_bayes4 <- make_bayes_ci90(
  fit_hier$draws("beta") |> posterior::as_draws_matrix(),
  taxa = colnames(rel_abundances2),
  method_label = "bayes4"
)
ci90_zi <- make_bayes_ci90(beta_zi, taxa = taxa_names, method_label = "zi_ordinal")

res_bayes4 <- res_bayes4 |>
  left_join(ci90_bayes4, by = c("taxon","method")) |>
  mutate(sigCI90 = (lwr90 > 0) | (upr90 < 0))
res_zi <- res_zi |>
  left_join(ci90_zi, by = c("taxon","method")) |>
  mutate(sigCI90 = (lwr90 > 0) | (upr90 < 0))

res_all <- bind_rows(
  res_orm |> add_sig_flags(),
  res_ancombc |> add_sig_flags(),
  res_linda   |> add_sig_flags(),
  res_bayes4,
  res_zi
)

# =============================================================================
# 12) Summary tables that obey supervisor’s matching-rule advice
# =============================================================================
# Counts of discoveries at matched levels
summary_sig <- bind_rows(
  # Frequentists
  res_all |> filter(method %in% c("orm","ancombc","linda")) |>
    summarise(
      method = first(method),
      n_sig_q10 = sum(sig_q10, na.rm = TRUE),
      n_sig_q05 = sum(sig_q05, na.rm = TRUE),
      .by = method
    ),
  # Bayes
  res_all |> filter(method %in% c("bayes4","zi_ordinal")) |>
    summarise(
      method = first(method),
      n_sig_q10 = sum(sigCI90, na.rm = TRUE),  # 90 percent CI excludes 0
      n_sig_q05 = sum((lwr95 > 0) | (upr95 < 0), na.rm = TRUE),  # 95 percent CI excludes 0
      .by = method
    )
) |>
  arrange(method)

print(summary_sig)

# Interval width comparison to note "Bayesian intervals narrower"
ci_widths <- res_all |>
  mutate(ciw95 = upr95 - lwr95) |>
  group_by(method) |>
  summarise(
    median_ciw95 = median(ciw95, na.rm = TRUE),
    p25_ciw95    = quantile(ciw95, 0.25, na.rm = TRUE),
    p75_ciw95    = quantile(ciw95, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
print(ci_widths)

# =============================================================================
# 13) Plots: per-method intervals and bar chart of significant counts
# =============================================================================
# Panel of intervals per method (scales may differ across methods)
plot_df <- res_all |>
  mutate(taxon = factor(taxon, levels = taxa_names)) |>
  select(method, taxon, estimate, lwr95, upr95)

ggplot(plot_df, aes(x = estimate, y = taxon)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_errorbarh(aes(xmin = lwr95, xmax = upr95), height = 0, alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ method, scales = "free_x") +
  labs(x = "Effect scale (method specific)", y = NULL,
       title = "Per-taxon estimates with 95 percent intervals") +
  theme_light() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Bar chart of discoveries at matched thresholds
sig_counts_long <- summary_sig |>
  tidyr::pivot_longer(cols = c(n_sig_q10, n_sig_q05),
                      names_to = "threshold", values_to = "n_sig") |>
  mutate(threshold = recode(threshold,
                            n_sig_q10 = "q < 0.10  vs  CI90",
                            n_sig_q05 = "q < 0.05  vs  CI95"))

ggplot(sig_counts_long, aes(x = method, y = n_sig, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.7)) +
  labs(x = NULL, y = "Number of significant taxa",
       title = "Matched significance: frequentist q vs Bayesian credible interval") +
  theme_minimal()

# Optional: table of Bayesian tail-prob q_bayes for reference
bayes_tail_pvals <- res_all |>
  filter(method %in% c("bayes4","zi_ordinal")) |>
  transmute(method, taxon, q_bayes = 2 * pmin(p_lt_0, p_gt_0))
head(bayes_tail_pvals)
