# =============================================================================
# Test script - Bayesian vs frequentist DAA on the example microbiome dataset
# Pipeline:
#   1) Load relative abundances and encode into zero-inflated ordinal categories
#   2) Vectorize data for Stan + prepare initial values
#   3) Fit hierarchical zero-inflated ordered-logit model using cmdstanr with parallel computing (~20 mins)
#   4) Extract Bayesian taxon-level effects and sanity-check centering
#   5) Fit per-taxon frequentist ordered logit (rms::orm) as a baseline
#   6) Join and visualize effect estimates with 95 percent intervals
# Notes:
#   - Zero category encodes structural zeros, nonzero values are binned by quantiles
#   - Group is coded 0 control, 1 case, consistent with the Stan file
# =============================================================================
library(tidyverse)
library(cmdstanr)
library(posterior)
library(rms)

# Load & prep data
load("example_mb_data.rds")
# Subset only taxon columns into a tibble
rel_abundances <- mb_data |> select(contains('taxon_'))

# Map relative abundances to ordinal categories with an explicit zero bin
# Output is integer in {1, ..., K} where:
#   1 = structural zero
#   2..(K_nonzero+1) = quantile bins among nonzero observations
safe_quartile_bins_zero <- function(x, K_nonzero = 3) {
  x <- as.numeric(x)
  z <- !is.na(x) & x == 0
  out <- integer(length(x)); out[z] <- 1L # structural zeros -> category 1
  nx <- x[!z] # nonzero values
  if (length(nx) == 0) return(out)  # all zero
  # Equi-probability cutpoints among nonzeros
  qs <- quantile(nx, probs = seq(0, 1, length.out = K_nonzero + 1), na.rm = TRUE, type = 8)
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


# Minimal init that helps the sampler by seeding global cutpoints near equal-probability logits
# delta_c starts at zero shift for all taxa
stan_init <- function() {
  list(
    c_base  = qlogis((1:(K-1))/K),  
    delta_c = rep(0, M)             
    # Other parameters use cmdstanr random inits
  )
}

# Compile hierarchical Stan model, enable threads for parallel reduce_sum
mod_hier <- cmdstan_model("hierarchical_model_test.stan", cpp_options = list(stan_threads = TRUE))

# Run the hierarchical model once for this test dataset
start_time_hier <- Sys.time()
fit_hier <- mod_hier$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 3,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1,
  init = stan_init,
  refresh = 0,
  show_messages = FALSE
)
end_time_hier <- Sys.time()
cat("Hierarchical model runtime:", end_time_hier - start_time_hier, "\n")

# beta_centered_mean is beta[m] minus the within-draw mean, for comparability across taxa
beta_cmean_mat   <- fit_hier$draws("beta_centered_mean") |> posterior::as_draws_matrix()
beta_center_mean <- fit_hier$draws("beta_center_mean")   |> as.vector()

# QA: the within-draw mean of centered betas should be near zero
summary(rowMeans(beta_cmean_mat))
hist(beta_center_mean, breaks = 40,
     main = "Per-draw mean center used in Stan (QA)",
     xlab = "beta_center_mean")


# Extract posterior draws for beta[m] and compute medians, 95 percent intervals, and tail probs
draws_hier <- fit_hier$draws(variables = c("beta")) |> as_draws_df()
# Each beta[1], beta[2], ... is for a taxon

results_hier <- draws_hier %>%
  select(starts_with("beta[")) %>%
  pivot_longer(everything(), names_to = "taxon", values_to = "beta") %>%
  group_by(taxon) %>%
  summarise(
    est = median(beta),
    lwr95 = quantile(beta, 0.025),
    upr95 = quantile(beta, 0.975),
    p_lt_0 = mean(beta < 0),
    p_gt_0 = mean(beta > 0)
  ) %>%
  mutate(taxon = paste0("taxon_", row_number()))

print(results_hier)

# Fit rms::orm for a single taxon abundance vs group
# Returns coefficient for the group effect, its SE, and a score-test p-value
run_orm <- function(abundance, metadata, formula) {
  # Build model matrix from metadata based on 'formula', then append response
  mm <- model.matrix(formula, metadata) |> cbind(abundance) |> as_tibble() |> select(-"(Intercept)")
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  score_1 <- fit_1$stats["Score"]
  res_freq <- data.frame(estimate = fit_1$coefficients[vars],
                         se = sqrt(diag(vcov(fit_1))[vars]), p_value = NA)
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
  )

print(res_freq)

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
  facet_wrap(~ model, scales = "free_x") +
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
  facet_wrap(~ model, scales = "free_x") +
  labs(x = "Parameter Value", y = NULL,
       title = "Compare Estimates and 95% Intervals by Model (Trimmed Frequentist Outliers)",
       color = "Model") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
