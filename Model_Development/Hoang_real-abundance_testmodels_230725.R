# Load the libraries
library(tidyverse)
library(cmdstanr)
library(posterior)
library(rms)

# A function to transform real numbers into ordinal categories
create_cats <- function(y) {
  factor(y, levels = sort(unique(y)), labels = 1:length(unique(y))) |>
    as.integer()
}

# Prepare data
load("example_mb_data.rds")

rel_abundances <- mb_data |> select(contains('taxon_'))

rel_cats <- apply(rel_abundances, 2, create_cats)

group <- ifelse(mb_data$group == 'case', 1, 0)

group_0 <- as.numeric(group == 0)
group_1 <- as.numeric(group == 1)

# Convert each taxon column to ordinal categories
rel_cats <- apply(rel_abundances, 2, create_cats)
n_cats <- apply(rel_cats, 2, max)

# Generate shared cutpoints from equal cumulative probabilities
K_max <- max(n_cats)
p <- rep(1/K_max, K_max)
cumprobs <- cumsum(p)[1:(K_max - 1)]
quantiles <- qlogis(cumprobs)

# --- Simplified Bayesian model ---
mod <- cmdstan_model("stan_ordinal_simplified_hoang_230625.stan", cpp_options = list(stan_threads = TRUE))

# Track execution time
start_time <- Sys.time() 

# Run model for each taxon
results <- map2_dfr(
  .x = as.data.frame(rel_cats),
  .y = n_cats,
  .id = "taxon",
  .f = function(y, K) {
    # Stan data per taxon
    d <- list(
      N = length(y),
      K = K,
      y = y,
      group_0 = group_0,
      group_1 = group_1,
      c = quantiles[1:(K - 1)]
    )
    
    # Sample from model
    fit <- mod$sample(
      data = d,
      chains = 4,
      parallel_chains = 4,
      threads_per_chain = 3,
      iter_warmup = 500, # Reduce to 500 for faster
      iter_sampling = 500, # Reduce to 500 for faster,
      refresh = 0,
      show_messages = FALSE,
      seed = 1
    )
    
    end_time <- Sys.time()
    execution_time <- end_time - start_time
    print(paste("Execution time:", execution_time))
    
    # Extract difference beta_1 - beta_0
    draws <- fit$draws(variables = c("beta_0", "beta_1")) |> as_draws_df()
    diff <- draws$beta_1 - draws$beta_0
    
    tibble(
      est = median(diff),
      lwr95 = quantile(diff, 0.025),
      upr95 = quantile(diff, 0.975),
      p_lt_0 = mean(diff < 0),
      p_gt_0 = mean(diff > 0)
    )
  }
)

print(results)

# 15 minutes to run 

# --- Improved Bayesian model ---
mod_imp <- cmdstan_model("stan_ordinal_improved_tests_parameterization.stan", cpp_options = list(stan_threads = TRUE))

# Track execution time
start_time_imp <- Sys.time()

# Run model for each taxon
results_imp <- map2_dfr(
  .x = as.data.frame(rel_cats),
  .y = n_cats,
  .id = "taxon",
  .f = function(y, K) {
    # Stan data per taxon
    d_imp <- list(
      N = length(y),
      K = K,
      y = y,
      group = group,  # pass the binary group vector (0 = control, 1 = case)
      c = quantiles[1:(K - 1)]
    )
    
    # Sample from improved model
    fit_imp <- mod_imp$sample(
      data = d_imp,
      chains = 4,
      parallel_chains = 4,
      threads_per_chain = 3,
      iter_warmup = 500,
      iter_sampling = 500,
      seed = 1,
      refresh = 0,
      show_messages = FALSE
    )
    end_time_imp <- Sys.time()
    print(paste("Improved model runtime:", end_time_imp - start_time_imp))
    
    # Extract and summarize beta (group effect)
    draws_imp <- fit_imp$draws(variables = c("beta")) |> as_draws_df()
    beta_samps <- draws_imp$beta
    
    tibble(
      est = median(beta_samps),
      lwr95 = quantile(beta_samps, 0.025),
      upr95 = quantile(beta_samps, 0.975),
      p_lt_0 = mean(beta_samps < 0),
      p_gt_0 = mean(beta_samps > 0)
    )
  }
)

results_imp
# 

# --- Previously Improved Bayesian model ---
mod_imp2 <- cmdstan_model("stan_ordinal_improved_2nd_hoang_230725.stan", cpp_options = list(stan_threads = TRUE))

# Track execution time
start_time_imp2 <- Sys.time()

# Run model for each taxon
results_imp2 <- map2_dfr(
  .x = as.data.frame(rel_cats),
  .y = n_cats,
  .id = "taxon",
  .f = function(y, K) {
    # Stan data per taxon
    d_imp2 <- list(
      N = length(y),
      K = K,
      y = y,
      group_0 = group_0,
      group_1 = group_1,
      c = quantiles[1:(K - 1)]
    )
    
    # Sample from improved model
    fit_imp2 <- mod_imp2$sample(
      data = d_imp2,
      chains = 4,
      parallel_chains = 4,
      threads_per_chain = 3,
      iter_warmup = 500,
      iter_sampling = 500,
      seed = 1,
      refresh = 0,
      show_messages = FALSE
    )
    end_time_imp2 <- Sys.time()
    print(paste("Improved model runtime:", end_time_imp2 - start_time_imp2))
    # Extract and summarize beta_1 - beta_0
    draws_imp2 <- fit_imp2$draws(variables = c("beta_0", "beta_1")) |> as_draws_df()
    diff_imp2 <- draws_imp2$beta_1 - draws_imp2$beta_0
    
    tibble(
      est = median(diff_imp2),
      lwr95 = quantile(diff_imp2, 0.025),
      upr95 = quantile(diff_imp2, 0.975),
      p_lt_0 = mean(diff_imp2 < 0),
      p_gt_0 = mean(diff_imp2 > 0)
    )
  }
)

results_imp2

# --- Frequentist model using rms::orm ---
run_orm <- function(abundance, metadata, formula) {
  mm <- model.matrix(formula, metadata) |> cbind(abundance) |> as_tibble() |> select(-"(Intercept)")
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  score_1 <- fit_1$stats["Score"]
  res_freq <- data.frame(estimate = fit_1$coefficients[vars],
                         se = sqrt(diag(vcov(fit_1))[vars]), p_value = NA)
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
meta <- data.frame(group = factor(mb_data$group, levels = c("control", "case")))
res_freq <- rel_abundances |> 
  map(~ run_orm(., metadata = mb_data, formula = ~ group)) |> 
  bind_rows(.id = 'taxon') |> 
  mutate(
    ci_lwr_95 = estimate - qnorm(.975) * se,
    ci_upr_95 = estimate + qnorm(.975) * se
  )

print(res_freq)


# --- Clean and rename each result set ---

# Frequentist results
res_freq_clean <- res_freq |>
  select(taxon, 
         frequentist_estimate = estimate,
         frequentist_lwr = ci_lwr_95,
         frequentist_upr = ci_upr_95)

# Simplified Bayesian results
res_simplified_clean <- results |>
  select(taxon,
         simplified_estimate = est,
         simplified_lwr = lwr95,
         simplified_upr = upr95)

# Improved Bayesian results
res_imp_clean <- results_imp |>
  select(taxon,
         improved_estimate = est,
         improved_lwr = lwr95,
         improved_upr = upr95)

# Improved Bayesian results (between)
res_imp2_clean <- results_imp2 |>
  select(taxon,
         improved2_estimate = est,
         improved2_lwr = lwr95,
         improved2_upr = upr95)

# --- Join all three into one comparison table ---
comparison_table <- res_freq_clean |>
  left_join(res_simplified_clean, by = "taxon") |>
  left_join(res_imp_clean, by = "taxon") |>
  left_join(res_imp2_clean, by = "taxon")

comparison_table <- comparison_table |>
  select(
    taxon,
    frequentist_estimate, simplified_estimate, improved_estimate, improved2_estimate,
    frequentist_lwr,      simplified_lwr,      improved_lwr, improved2_lwr,
    frequentist_upr,      simplified_upr,      improved_upr, improved2_upr
  )
comparison_table

# We generate a plot to compare the Bayesian and frequentist versions in term of credible and confidence intervals. 
# The y-axis in each model represent for each taxon. I decided to hide this so the plot will be easier to look at (since we are comparing CI, no need the taxon name).
# Reshape to long format
comparison_long <- comparison_table |>
  pivot_longer(cols = -taxon,
               names_to = c("model", ".value"),
               names_pattern = "(frequentist|simplified|improved|improved2)_(estimate|lwr|upr)")

# Plot
ggplot(comparison_long, aes(x = estimate, y = taxon, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ model, scales = "free_x") +
  labs(x = "Parameter Value", y = NULL,
       title = "Compare Estimates and 95% Intervals by Model",
       color = "Model") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # hides taxon names
    axis.ticks.y = element_blank()  # hides tick marks on y-axis
  )

# Since we observed some outliers from the frequentist model, I decided to trim these outliers, so it will be easier to compare the interval of those models.
# Add CI width and filter out outliers
comparison_long_trimmed <- comparison_long |>
  mutate(ci_width = upr - lwr) |>
  filter(ci_width <= 20)  # adjust as needed

# Plot without outlier CIs
ggplot(comparison_long_trimmed, aes(x = estimate, y = taxon, color = model)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ model, scales = "free_x") +
  labs(x = "Parameter Value", y = NULL,
       title = "Compare Estimates and 95% Intervals by Model (Trimmed)",
       color = "Model") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # hides taxon names
    axis.ticks.y = element_blank()  # hides tick marks on y-axis
  )

# Findings
# First, it takes me 15 minutes in the simplified models (even though this is a very simple model with 2 parameters.
# Although each Stan fit reports only ~1.2s of raw sampler time for the 4 parallel chains, the real-time takes 3-4s.
# This incurs overhead that is not counted in the chain timer. Those milliseconds accumulate when the loop is executed 256 times.
# Secondly, for the improved model (added 2 scalar parameters), it takes 18 minutes to run, which is acceptable.
# When comparing these 3 CI, The frequentist orm fits occasionally produce huge confidence bounds (+/- 400 on log-odd scale).
# This happens whenever a taxon is rare or absent in one group, so we need to trim these extreme cases.
# The Bayesian models avoid this by placing weakly‑informative priors on the coefficients
# They pull the estimates back toward zero and yield finite, more stable intervals.

# When comparing these models:
# Frequentist (red) points cluster near zero, but the error bars are still the widest of the three.
# Improved Bayesian (green) intervals are the narrowest (maybe the best among 3 until now?), thanks to the extra global shift/scale.
# Simplified Bayesian (blue) sits between the two: intervals are regular and almost always overlap zero.
# So, I think the Bayesian approach delivers more stable interferences, and the improved model relatively improve the simplified model.