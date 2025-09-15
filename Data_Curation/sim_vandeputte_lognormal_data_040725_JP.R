library(tidyverse)

# Means, variances and correlations of absolute abundances based on Vandeputte data
load('mean_vars_vandeputte_250625.rds')
load('cors_vandeputte_250625.rds')

n <- 200 # Sample size
n_taxa <- nrow(mean_vars) # Number of taxa
groups <- rep(0:1, each = n / 2) # Two groups: control (0) and case (1)

# Simulate microbiome data -----------------------------------------------------
set.seed(1)

# First, create latent variables from multivariate standard normal distribution
mvn_latent <- MASS::mvrnorm(n, mu = rep(0, n_taxa), Sigma = cors)

# Create true effects
non_zero_effect <- rbinom(n_taxa, size = 1, prob = .1)
effects <- rnorm(n_taxa, mean = 0, sd = .5) * non_zero_effect
names(effects) <- paste0("taxon_", 1:n_taxa)

# Create effects matrix
effects_matrix <- (cbind(replicate(n, effects)) |> t()) * groups

# Add effects to latent variables
mvn_latent_plus_effect <- mvn_latent + effects_matrix

# Transform normal variates to log-normal variates. Log-normal distributions can
# resemble the distributions of absolute abundances, accounting for skewness.
# Parameters for the log-normal distributions (from Vandeputte data)
abd_means <- mean_vars$mean
abd_vars <- mean_vars$var

# Calculate log-normal parameters: meanlog and sdlog
lognorm_meanlog <- log(abd_means^2 / sqrt(abd_vars + abd_means^2))
lognorm_sdlog <- sqrt(log(abd_vars / abd_means^2 + 1))

m_lognorm_meanlog <- rbind(replicate(n, lognorm_meanlog)) |> t()
m_lognorm_sdlog <- rbind(replicate(n, lognorm_sdlog)) |> t()

# Transform multivariate normal variates to uniform variates and then to
# log-normal variates that are the simulated absolute abundances
mvu_latent <- pnorm(mvn_latent_plus_effect)
true_abundances <- qlnorm(mvu_latent,
                          meanlog = m_lognorm_meanlog,
                          sdlog = m_lognorm_sdlog)

# The measurement process etc introduce taxon-wise biases to the measured abundances
taxonwise_biases <- exp(rnorm(n_taxa, mean = 0, sd = 1))
biased_abundances <- t(taxonwise_biases * t(true_abundances))

# The total number of observed counts for each sample does not reflect the
# total amount of DNA in the sample (absolute abundance). Only relative
# information of abundances is gathered, that is, microbiome data are
# compositional.
relative_biased_abundances <- biased_abundances / rowSums(biased_abundances)

# Library size = the total number of counts in each sample
library_sizes <- 10 ^ rnorm(n, mean = 4.0, sd = .20) |> round()

# Simulate count data
lambdas <- library_sizes * relative_biased_abundances
counts <- rpois(n * n_taxa, lambda = lambdas) |> 
  matrix(nrow = n, ncol = n_taxa)
colnames(counts) <- names(effects)

# Perform differential abundace analysis----------------------------------------

# Transform counts to relative abundances
rel_abs <- counts / rowSums(counts)

# Remove taxa < 3 non_zero counts
rel_abundances <- rel_abs[, colSums(rel_abs > 0) >= 3]
# Function to run frequentist ORM for all taxa
run_orm <- function(abundance, metadata, formula) {
  
  mm <- model.matrix(formula, metadata) |> 
    cbind(abundance) |> 
    as_tibble() |> 
    select(-"(Intercept)")
  
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  score_1 <- fit_1$stats["Score"]
  
  res_freq <- data.frame(
    estimate = fit_1$coefficients[vars],
    se = sqrt(diag(vcov(fit_1))[vars]),
    p_value = NA
  )
  
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


meta <- data.frame(group = factor(groups, labels = c("control", "case")))

res <- rel_abundances |> 
  as.data.frame() |>
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')

# Compare the results to the true effects
true_effects <- tibble(taxon = names(effects),
                       true_effect = effects)

results <- res |> 
  left_join(true_effects, by = "taxon") |>
  mutate(q = p.adjust(p_value, method = "fdr"),
         prevalence = colMeans(rel_abundances > 0),
         significant = q < .05,
         true_positive = sign(estimate) == sign(true_effect) & significant,
         true_negative = true_effect == 0 & !significant,
         false_positive = sign(estimate) != sign(true_effect) & significant,
         false_negative = true_effect != 0 & !significant)

results |> summarize(n_taxa = n(),
                     n_significant = sum(significant),
                     n_true_positive = sum(true_positive),
                     n_true_negative = sum(true_negative),
                     n_false_positive = sum(false_positive),
                     n_false_negative = sum(false_negative),
                     fdr = 1 - n_true_positive / n_significant,
                     power = n_true_positive /
                       (n_true_positive + n_false_negative))

# FDR = False Discovery Rate
