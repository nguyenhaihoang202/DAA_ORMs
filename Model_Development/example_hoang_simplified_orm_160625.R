library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


N <- 10 # Number of observations

# Group variable and indicator variable fo each group
group <- rep(0:1, each = N / 2)
group_0 <- group == 0
group_1 <- group == 1

# A random ordinal variable with distinct values
set.seed(1)
y <- sample(1:N)

# Calculate quantile limits for logistic distribution so that the classes have
# equal probabilities.
k <- length(unique(y)) # Number of classes
probs <- (1:(k - 1)) / k 
quantiles <- qlogis(probs)


data_stan <- list(
  N = N,
  K = k,
  y = y,
  group_0 = group_0,
  group_1 = group_1,
  c = quantiles
)


# A super simplified Bayesian model where the intercepts determining the
# quantiles of the logistic distribution are fixed! Now there is only two
# parameters to be estimated, one for each group.
mod <- stan_model('stan_ordinal_simplified_170625.stan')

fit <- sampling(mod,
                data = data_stan,
                iter = 2000,
                warmup = 1000,
                chains = 4)

print(fit)

# Extract samples and calcullate estimate and 95% credible interval limits for
# the difference of betas
samples_beta <- extract(fit)$beta_1 - extract(fit)$beta_0

res_bayes <- c(
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
res_bayes
res_freq
