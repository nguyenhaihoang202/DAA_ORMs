library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


N <- 10 # Number of observations

# Intercept, group variable and indiator variable fo each group
intercept <- rep(1, N)
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
  intercept = intercept,
  group = group,
  c = quantiles
)



mod <- stan_model('stan_ordinal_betapara_globalprior_070825.stan')

fit <- sampling(mod,
                data = data_stan,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                seed = 1)

print(fit)

# Extract samples and calculate estimate and 95% credible interval for beta
samples_beta <- extract(fit)$beta_diff

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
