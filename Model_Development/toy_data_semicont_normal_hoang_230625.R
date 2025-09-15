# Load the libraries
library(tidyverse)
library(cmdstanr)
library(posterior)
library(rms)

# Generate toy data with semi-continuous variables and there are true differences between the group

# Set seed 
set.seed(1)

# Main parameters
N <- 1000 # Number of observations
prop_zero <- 0.5 # Proportion of observations in the first class (zero)
group_diff <- 2 # True difference in latent means between groups

# Group variable and indicator variable for each group
group <- rep(0:1, each = N / 2)
group_0 <- group == 0
group_1 <- group == 1

# Generate a latent continuous variable to discretizing observations based on it's quantile 
latent <- rnorm(N, mean = ifelse(group == 0, 0, group_diff), sd = 1)

# Define classes
K <- 10 # Number of classes, pre-defined 
n_zeros <- round(N * prop_zero) # Observations in first class
n_nonzero <- N - n_zeros # Observations in other classes

# Assigning observations to classes
y <- rep(NA, N)
zero_indices <- sample(1:N, n_zeros)  
y[zero_indices] <- 1  

nonzero_indices <- setdiff(1:N, zero_indices)
nonzero_latent <- latent[nonzero_indices]

# Quantiles generation
quantiles_latent <- quantile(nonzero_latent, probs = seq(0, 1, length.out = K))
y[nonzero_indices] <- cut(nonzero_latent, breaks = quantiles_latent, labels = 2:K, include.lowest = TRUE) |> as.numeric()
y <- as.integer(y)

# Threshold calculation
p <- numeric(K)
p[1] <- prop_zero
p[2:K] <- (1 - prop_zero) / (K - 1)
cumprobs <- cumsum(p)[1:(K - 1)]
quantiles <- qlogis(cumprobs)

# Stan data list
stan_data <- list(
  N   = N,
  K   = K,
  y   = y,
  c   = quantiles,    
  group_0 = as.numeric(group_0),
  group_1 = as.numeric(group_1)
)

# Fit the simplified model
mod_simplified <- cmdstan_model('stan_ordinal_simplified_hoang_230625.stan', cpp_options = list(stan_threads = TRUE))

fit_simplified <- mod_simplified$sample(
  data              = stan_data,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 4,
  iter_warmup       = 1000,
  iter_sampling     = 1000,
  seed              = 1
)

print(fit_simplified$summary())

# Extract samples and calculate estimate and 95% credible interval limits for
# the difference of betas
samples_beta <- fit_simplified$draws("beta_1") - fit_simplified$draws("beta_0")

res_bayes_simplified <- c(
  est = median(samples_beta),
  lwr = quantile(samples_beta, 0.025) |> as.numeric(),
  upr = quantile(samples_beta, 0.975) |> as.numeric()
)


# Run frequentist ordinal regression model
fit_form <- rms::orm(y ~ group)

est <- fit_form$coefficients['group'] |> as.numeric()
se <- sqrt(diag(vcov(fit_form))['group']) |> as.numeric()

res_freq <- c(
  est = est,
  lwr = est - 1.96 * se,
  upr = est + 1.96 * se)


# Compare the results
res_bayes_simplified
res_freq

# Both models use the logistic link, so they produce logistic-scale results.
# The model’s beta_diff (difference in log-odds, beta_1 - beta_0) is estimated on the logistic scale, so the expected beta_diff is ~1.1 from (2/1.81)
# However, the simplified model can still be improved (we can see that the lwr and upr does not included 1.1), the performance is not as good as the frequentist.
# With the bayesian model’s MCMC sampling is robust (1.0 R_hat, high ESS).
# The frequentist ORM estimates a group effect of 0.97, but the CI included 1.1, which is quite good.
# When we try to reduce the number prop_zero, both performance estimates better, which is expected.

# Fit the improved model
mod_improved <- cmdstan_model('stan_ordinal_improved_hoang_230625.stan', cpp_options = list(stan_threads = TRUE))

fit_improved <- mod_improved$sample(
  data              = stan_data,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 4,
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
res_freq
res_bayes_imp

# The simplified Bayesian model slightly off-track because its rigid cut-points cannot track the spacing of the ordinal bins.
# However, the shift-stretch improved model is enough to bring it into line, and then gave the approximately same answer to the frequentist approaches.