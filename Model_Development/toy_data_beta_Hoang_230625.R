# Load the libraries
library(tidyverse)
library(cmdstanr)
library(posterior)
library(rms)

# Generate toy data with beta-distributed latent variable and true differences between the group

# Set seed 
set.seed(1)

# Main parameters
N <- 1000 # Number of observations
K <- 10 # Number of classes, pre-defined 

# Group variable and indicator variable for each group
group <- rep(0:1, each = N / 2)
group_0 <- group == 0
group_1 <- group == 1

# Generate a latent continuous variable to discretizing observations based on it's quantile 
latent <- rbeta(N, shape1 = ifelse(group == 0, 1, 3), shape2 = 5) * 100 # Scale to 0–100

# The latent difference (0.208 unscaled) shifts the mean of the beta-distributed latent variable.
# After scaling (×100), the difference is 20.8, but the ordinal bins compress this effect.
# In ordinal regression, a latent mean difference of ~20.8 typically translates to a log-odds difference of 2 (I keep 2 to be similar to the semicont_normal)

# Assigning observations to classes
y <- cut(latent, breaks = quantile(latent, probs = seq(0, 1, length.out = K + 1)), labels = 1:K, include.lowest = TRUE) |> as.numeric()
y <- as.integer(y)

# Threshold calculation
p <- rep(1/K, K) # Equal probabilities
cumprobs <- cumsum(p)[1:(K - 1)]
quantiles <- qlogis(cumprobs)

# Stan data list
stan_data <- list(
  N = N,
  K = K,
  y = y,
  c = quantiles,
  group_0 = as.numeric(group_0),
  group_1 = as.numeric(group_1)
)

# Fit the simplified model
mod_simplified <- cmdstan_model('stan_ordinal_simplified_hoang_230625.stan', cpp_options = list(stan_threads = TRUE))

fit_simplified <- mod_simplified$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 1
)

print(fit_simplified$summary())

# Extract samples and calculate estimate and 95% credible interval limits for the difference of betas
samples_beta <- fit_simplified$draws("beta_1") - fit_simplified$draws("beta_0")
res_bayes_simplified <- c(
  est = median(samples_beta),
  lwr = quantile(samples_beta, 0.025),
  upr = quantile(samples_beta, 0.975)
)


# Extract samples and calculate estimate and 95% credible interval limits for the difference of betas
samples_beta_imp <- fit_improved$draws("beta_1") - fit_improved$draws("beta_0")
res_bayes_imp <- c(
  est = median(samples_beta_imp),
  lwr = quantile(samples_beta_imp, 0.025),
  upr = quantile(samples_beta_imp, 0.975)
)

# Run frequentist ordinal regression model
fit_form <- rms::orm(y ~ group)
est <- fit_form$coefficients['group']
se <- sqrt(diag(vcov(fit_form))['group'])
res_freq <- c(est = est, lwr = est - 1.96 * se, upr = est + 1.96 * se)

# Compare the results
res_bayes_simplified
res_freq

# Fit the improved model
mod_improved <- cmdstan_model('stan_ordinal_improved_hoang_230625.stan', cpp_options = list(stan_threads = TRUE))

fit_improved <- mod_improved$sample(
  data              = stan_data,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 8,
  iter_warmup       = 1000,
  iter_sampling     = 1000,
  seed              = 1
)

print(fit_improved$summary())

# Extract samples and calcullate estimate and 95% credible interval limits for the difference of betas
# Extract draws
draws <- fit_improved$draws(variables = c("beta_1", "beta_0", "lambda"))
draws_mat <- posterior::as_draws_matrix(draws)

# Logistic-scale difference (raw estimate from model)
samples_beta_logistic <- draws_mat[,"beta_1"] - draws_mat[,"beta_0"]

res_bayes_imp <- c(
  est = median(samples_beta_logistic),
  lwr = quantile(samples_beta_logistic, 0.025),
  upr = quantile(samples_beta_logistic, 0.975)
)

# Compare the results
res_bayes_simplified
res_bayes_imp
res_freq

# Using the simplified Bayesian model—whose cut-points are fixed—gives an estimate of 1.82 [1.61, 2.03], slightly low.
# Once the improved Stan model adds a global shift and global stretch for the cut-points, the posterior moves to 2.15, nearly the same with the frequentist fit, the same phenomenon we observed with the semicont_normal example.