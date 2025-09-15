# Preparing the libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(tidyverse)
library(ggpubr)
library(brms)

# Load the microbiome data
load("example_mb_data.rds")

# Take the relative abudances of taxon 2 as an example
ggplot(mb_data, aes(x = taxon_2)) +
  geom_histogram(bins = 500) +
  facet_grid(group ~ .)

# Transformed to categories based on their ordinal value.
taxon_2_ord <- factor(mb_data$taxon_2,
                      levels = sort(unique(mb_data$taxon_2)),
                      labels = 1:length(unique(mb_data$taxon_2))) |> 
  as.integer()

# Check the distribution of the ordinal values
mb_data$taxon_2
taxon_2_ord

# Using rstan to run Bayesian ordinal regression model

# Convert group to binary (0 = control, 1 = case)
group <- ifelse(mb_data$group == "case", 1, 0)

# Prepare data for Stan
stan_data <- list(
  N = length(taxon_2_ord),
  K = max(taxon_2_ord),
  y = taxon_2_ord,
  group = group
)

# Define Stan model as a string
stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1, upper=K> y[N];
  int<lower=0, upper=1> group[N];
}
parameters {
  real beta;
  ordered[K - 1] c;
}
model {
  beta ~ normal(0, 5);
  c ~ normal(0, 5);
  for (n in 1:N)
    y[n] ~ ordered_logistic(beta * group[n], c);
}
"

# Compile and sample from the model
stan_model <- stan_model(model_code = stan_code)

fit_rstan <- sampling(stan_model,
                data = stan_data,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                seed = 1)

# Show results
print(fit_rstan, pars = c("beta", "c"))


# Fit using brms with taxon 2 to compare the results
mbd <- tibble(group = mb_data$group,
              y = taxon_2_ord)

brm2 <- brm(y ~ group,
            data = mbd,
            family = cumulative("logit"),
            iter = 2000,
            warmup = 1000,
            chains = 4,
            seed = 1)

# Print the results
print(brm2)

# There seems to no real difference between the groups (Q2.5 -> Q97.5 not includes 0)
brms::posterior_summary(brm2, variable = "b_groupcase")


