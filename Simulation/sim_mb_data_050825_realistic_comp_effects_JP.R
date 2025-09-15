library(tidyverse)

# Means, variances and correlastions of absolute abundances based on Vandeputte
# data
load('mean_vars_vandeputte_250625.rds')
load('cors_vandeputte_250625.rds')

n <- 100 #Sample size
n_taxa <- nrow(mean_vars) # Number of taxa
groups <- rep(0:1, each = n / 2) # Two groups: control (0) and case (1)


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

# Parameters for the gamma distributions (from Vandeputte data)
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

colnames(counts) <- names(effects)




# Perform differential abundace analysis----------------------------------------

# Transform counts to relative abundances
rel_abs <- counts / rowSums(counts)

# Remove taxa < 3 non_zero counts
rel_abundances <- rel_abs[, colSums(rel_abs > 0) >= 5]

counts2 <- counts[, colSums(counts > 0) >= 5]

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

res_orm <- rel_abundances |> 
  as.data.frame() |>
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')


make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  meta_data <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, meta_data)
}

obj_ancombc <- ANCOMBC::ancombc(data = make_physeq(counts2, meta),
                                  formula = "group",
                                  p_adj_method = "fdr",
                                  prv_cut = 0)

obj_linda <- LinDA::linda(
  otu.tab = counts2 |> t(),
  meta = meta,
  formula = "~ group"
)

# Compare the results to the true effects
true_effects <- tibble(taxon = names(effects),
                       true_effect = effects)

res_ancombc <- tibble(taxon = obj_ancombc$res$lfc$taxon,
                      estimate = obj_ancombc$res$lfc$groupcase,
                      se = obj_ancombc$res$se$groupcase,
                      p = obj_ancombc$res$p_val$groupcase,
                      q = obj_ancombc$res$q_val$groupcase) |> 
  left_join(true_effects, by = "taxon") |>
  mutate(significant = q < .05,
         lwr = estimate - 1.96 * se,
         upr = estimate + 1.96 * se,
         true_positive = sign(estimate) == sign(true_effect) & significant,
         true_negative = true_effect == 0 & !significant,
         false_positive = sign(estimate) != sign(true_effect) & significant,
         false_negative = true_effect != 0 & !significant,
         method = 'ancombc')


res_orm2 <- res_orm |> 
  left_join(true_effects, by = "taxon") |>
  mutate(q = p.adjust(p_value, method = "fdr"),
         lwr = estimate - 1.96 * se,
         upr = estimate + 1.96 * se,
         prevalence = colMeans(rel_abundances > 0),
         significant = q < .05,
         true_positive = sign(estimate) == sign(true_effect) & significant,
         true_negative = true_effect == 0 & !significant,
         false_positive = sign(estimate) != sign(true_effect) & significant,
         false_negative = true_effect != 0 & !significant,
         method = 'orm')

res_linda <- obj_linda$output$groupcase |>
  rownames_to_column("taxon") |>
  select(taxon,
         estimate = log2FoldChange,
         se = lfcSE,
         p = pvalue,
         q = padj) |>
  left_join(true_effects, by = "taxon") |>
  mutate(lwr = estimate - 1.96 * se,
         upr = estimate + 1.96 * se,
         prevalence = colMeans(rel_abundances > 0),
         significant = q < .05,
         true_positive = sign(estimate) == sign(true_effect) & significant,
         true_negative = true_effect == 0 & !significant,
         false_positive = sign(estimate) != sign(true_effect) & significant,
         false_negative = true_effect != 0 & !significant,
         method = 'linda')


res_all <- bind_rows(res_orm2, res_ancombc, res_linda) |> 
  mutate(taxon = fct_reorder(taxon, -estimate, .na_rm = T))

rmd <- res_all |>
  group_by(method) |>
  summarize(md = median(estimate, na.rm = TRUE))

(p1 <- ggplot(res_all, aes(estimate, taxon)) +
  geom_point(aes(color = significant), size = 2) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
  geom_point(aes(x = true_effect), color = 'black', alpha = .5) +
  geom_vline(xintercept = 0, linetype = 'solid', color = 'grey30',
             linewidth = 0.5) +
  geom_vline(data = rmd, aes(xintercept = md),
             linetype = 'solid', color = 'blue', linewidth = 0.5) +
  coord_cartesian(xlim = c(-2, 3)) +
  labs(x = 'Effect size', y = 'Taxon') +
  facet_wrap(~ method) +
  theme_light() +
  theme(legend.position = "bottom"))

res_all |> 
  group_by(method) |>
  summarize(n_taxa = n(),
            n_significant = sum(significant),
            n_true_positive = sum(true_positive),
            n_false_positive = sum(false_positive),
            fdr = 1 - n_true_positive / n_significant) |> 
  select(method, n_taxa, n_false_positive, n_true_positive, fdr)

# Now when no true population effects are exactly zero it does not make sense
# to talk about true negatives or false negatives. Insted, the _number_ of true
# positives may be a good metric for power (and FDR a good error rate metric).
# Now: true positive = significant + correct sign
#      false positive = significant + opposite sign

# ORM seems to perform here very well compared to ANCOM-BC and LinDA!