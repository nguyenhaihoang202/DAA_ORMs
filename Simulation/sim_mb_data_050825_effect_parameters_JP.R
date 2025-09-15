library(tidyverse)

# Means, variances and correlastions of absolute abundances based on Vandeputte
# data
load('mean_vars_vandeputte_250625.rds')
load('cors_vandeputte_250625.rds')

n <- 65 + 28 # Sample size
n_taxa <- nrow(mean_vars) # Number of taxa
groups <- c(rep(0, 65), rep(1, 28)) # Two groups: control (0) and case (1)


# Simulate microbiome data -----------------------------------------------------
set.seed(1)

# First, create latent variables from multivariate standard normal distribution
mvn_latent <- MASS::mvrnorm(n, mu = rep(0, n_taxa), Sigma = cors)

# Create true effects based on Vandeputte's (2017! not 2021) study.
# Note that in real life there are no truly exactly zero effects.
effects <- rnorm(n_taxa, mean = -.9, sd = 0.6)

names(effects) <- paste0("taxon_", 1:n_taxa)

# Plot true effects
ggplot(data = tibble(x = effects)) +
  geom_histogram(aes(x = x), bins = 50) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black',
             linewidth = 1) +
  labs(x = 'Effect size', y = 'Count') +
  theme_light()

# Create effects matrix
effects_matrix <- (cbind(replicate(n, effects)) |> t()) * groups

# Add effects to latent variables
mvn_latent_plus_effect <- mvn_latent + effects_matrix

# Transform normal variates to gamma variates. Gamma distributions resemble
# the distributions of absolute abundances better than normal distributions.
# Gamma distributions may, however, not be very good either. You can try e.g.
# generalized gamma distributions or log-normal distributions, too, if you want.

# Parameters for the gamma distributions (from Vandeputte 2021 data)
abd_means <- mean_vars$mean
abd_vars <- mean_vars$var

gamma_shapes <- abd_means ^ 2 / abd_vars
m_gamma_shapes <- rbind(replicate(n, gamma_shapes)) |> t()

gamma_scales <- abd_vars / abd_means
m_gamma_scales <- rbind(replicate(n, gamma_scales)) |> t()

# Transform multivariate normal variates first to uniform variates and then to
# gamma variates that are the simulated absolute abundances
mvu_latent <- pnorm(mvn_latent_plus_effect)
true_abundances <- qgamma(mvu_latent,
                          shape = m_gamma_shapes,
                          scale = m_gamma_scales)

# The measurement process etc introduce taxon-wise biases to the measured
# abundances, i.e., different taxa are detected with different efficiencies.
taxonwise_biases <- exp(rnorm(n_taxa, mean = 0, sd = 1))

biased_abundances <- t(taxonwise_biases * t(true_abundances))


# The total number of observed counts for each sample does not reflect the
# total amount of DNA in the sample (absolute abundance). Only relative
# information of abundances is gathered, that is, microbiome data are
# compositional.
relative_biased_abundances <- biased_abundances / rowSums(biased_abundances)
rowSums(relative_biased_abundances) |> head()

# Library size = the total number of counts in each sample
library_sizes <- 10 ^ rnorm(n, mean = 4.0, sd = .20) |> round()

# Count data can be simulated by sampling from a Poisson distribution
lambdas <- library_sizes * relative_biased_abundances

counts <- rpois(n * n_taxa, lambda = lambdas) |> 
  matrix(nrow = n, ncol = n_taxa)

# Check overall prevalence
mean(counts > 0)

# Taxa with >= 10 non-zero counts. Remove them when analyzing the distribution
# of the absolute sample effects, to get rid of some outliers due to very small
# abundances.
prevalent_taxa <- colSums(counts[groups == 0, ] > 0) >= 5 &
  colSums(counts[groups == 1, ] > 0) >= 5

# Absolute sample effects. These should be similar to those _observed_ in the 
# Vandeputte 2017 data.
absolute_sample_effs <- cbind(group = groups, true_abundances[, prevalent_taxa]) |>
  as.data.frame() |>
  pivot_longer(cols = -group, names_to = 'taxon') |> 
  group_by(taxon, group) |>
  summarize(mean_abs = mean(value)) |>
  ungroup() |> 
  pivot_wider(names_from = group, values_from = mean_abs) |>
  mutate(eff = log(`1` / `0`))

# Plot true sample effects
ggplot(data = absolute_sample_effs, aes(x = eff)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black',
             linewidth = 1) +
  labs(x = 'True effect size (log(case / control))', y = 'Count') +
  theme_light()

absolute_sample_effs |> 
  summarize(m = mean(eff, na.rm = TRUE),
            md = median(eff, na.rm = TRUE),
            sd = sd(eff, na.rm = TRUE))


# When the true population effects are generated as
# effects <- rnorm(n_taxa, mean = -.9, sd = 0.6)
# it seems to lead to absoulute sample effects that have a distribution somewhat
# similar to the observed effects from Vandeputte 2021 data.
# (Note that the summary values depend on the seed.)