library(tidyverse)
library(cmdstanr)
library(posterior)

# A function to generate Stan code for multiple taxa
gen_stan_code <- function(ordered_c, prior_c, y_model){
  paste("data {
    int<lower=1> N;
    int<lower=1> K;
    array[N, K] int<lower=1> y;
    array[K] int<lower=2> J;
    vector[N] x;
  }

  parameters {
    vector[K] z;
    real <lower=0> tau; // Scale of the prior of beta parameters
    ", ordered_c, 
        "
  }

  transformed parameters {
    vector[K] beta;
    beta = z .* tau;
  }

  model {
    z ~ student_t(2, 0, 1);
    tau ~ normal(0, 1);
    ", prior_c, 
        "
    for (n in 1:N) {", y_model, "}
  }")
}


# A function to transform real numbers into ordinal categories
create_cats <- function(y){
  factor(y,
         levels = sort(unique(y)),
         labels = 1:length(unique(y))) |> 
    as.integer()
}

# Prepare data
load("example_mb_data.rds")

rel_abundances <- mb_data |> select(contains('taxon_'))

rel_cats <- apply(rel_abundances, 2, create_cats)

group <- ifelse(mb_data$group == 'case', 1, 0)

# Number of categories for each taxon
n_cats <- apply(rel_cats, 2, max)

# Data object for Stan
d_stan <- list(N = nrow(rel_abundances), # Number of samples
               K = ncol(rel_abundances), # Number of taxa
               y = rel_cats, # Ordinalized relative abundances
               J = n_cats, # Number of ord categories for each taxon
               x = group) # Group indicator

# Write Stan code
ordered_cs <- character(ncol(rel_cats))
for (k in 1:ncol(rel_cats)) {
  ordered_cs[k] <- paste0('ordered[J[', k, '] - 1] c_', k, ';')
}

ordered_c <- paste(ordered_cs, collapse = ' ')

prior_cs <- character(ncol(rel_cats))
for (k in 1:ncol(rel_cats)) {
  prior_cs[k] <- paste0('c_', k, ' ~ normal(0, 3);')
}

prior_c <- paste(prior_cs, collapse = ' ')

y_models <- character(ncol(rel_cats))
for (k in 1:ncol(rel_cats)) {
  y_models[k] <- paste0('y[n, ', k, '] ~ ordered_logistic(x[n] * beta[', k, '], c_', k, ');')
}

y_model <- paste(y_models, collapse = ' ')

stan_code <- gen_stan_code(ordered_c, prior_c, y_model)

# Track execution time
start_time <- Sys.time() 

# Using cmdstanr instead of rstan for performance and parallelism
# Compile the model, enabling multiple threads
stan_file <- write_stan_file(stan_code)
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE)) # Compiling with multithreading enabled

# Run Stan with parallel chains and multiple threads
n_chains <- 4 # run 4 chains
n_iter <- 200
n_warmup <- 100

# Run the model with cmdstanR
fit <- mod$sample(
  data = d_stan,
  chains = n_chains,
  iter_sampling = n_iter - n_warmup,
  iter_warmup = n_warmup,
  seed = 1,
  parallel_chains = n_chains,
  threads_per_chain = 4
)

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Execution time:", execution_time))

# Save the full model summary
mod_summary <- fit$summary()

# Extract draws for beta
beta_draws <- fit$draws(variables = "beta", format = "draws_df")

# Summarize
res <- beta_draws %>%
  pivot_longer(cols = starts_with("beta["), names_to = "name", values_to = "value") %>%
  mutate(n_taxon = str_extract(name, "\\d+") %>% as.numeric()) %>%
  group_by(n_taxon) %>%
  summarize(
    est = quantile(value, 0.500), # Posterior median
    lwr95 = quantile(value, 0.025),
    lwr90 = quantile(value, 0.050),
    lwr80 = quantile(value, 0.100),
    upr80 = quantile(value, 0.900),
    upr90 = quantile(value, 0.950),
    upr95 = quantile(value, 0.975),
    p_low = mean(value < 0),
    p_high = mean(value > 0)
  ) %>%
  mutate(taxon = colnames(rel_abundances))
res
