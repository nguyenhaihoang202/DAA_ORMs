library(tidyverse)
library(readxl)
library(MIDASim)

# Extract and preprocess count matrix from DataS1.xlsx
counts <- read_excel("DataS1.xlsx", sheet = "S1-6") %>% 
  rename(sid = '...1') %>% 
  mutate_all(~ replace(., is.na(.) | . == 'NA' | . == '', 0)) %>%  # Handle NA and empty strings
  mutate_all(as.numeric) %>% 
  mutate(reads = rowSums(select(., -sid), na.rm = TRUE))

meta <- read_excel("DataS1.xlsx", sheet = "S1-3") %>% 
  select(sid = '...1', day = Day_Number, id = ID_Number, load = Cell_count_per_gram) %>% 
  mutate_all(~ replace(., is.na(.) | . == 'NA' | . == '', 0)) %>% 
  mutate_all(as.numeric)

abundances <- meta %>% 
  left_join(., counts, 'sid') %>% 
  mutate(sf = reads / load, id = as_factor(id)) %>% 
  filter(!is.na(load) & reads > 10^3 & reads < 10^5 & load > 10^10 & day < 45) %>% 
  select(sf, contains('_')) %>% 
  mutate_at(vars(contains('_')), ~ . / sf) %>% 
  select(-sf) %>% 
  as.matrix()

# Use raw counts matrix (samples as rows, taxa as columns)
# JP: Remove the zero rows, too!
counts_1 <- counts %>% select(contains('_')) %>% as.matrix
counts_matrix <- counts_1[rowSums(counts_1) > 0, ]
  

# Apply MIDASim procedure (non-parametric mode)
midasim_setup_nonpara <- MIDASim.setup(counts_matrix, mode = "nonparametric")
midasim_modified_nonpara <- MIDASim.modify(midasim_setup_nonpara)
simulated.data_nonpara <- MIDASim(midasim_modified_nonpara)
summary(simulated.data_nonpara)

# Apply MIDASim procedure (parametric mode)
midasim_setup_para <- MIDASim.setup(counts_matrix, mode = "parametric")

# Same setup as Juho's approach
n <- dim(counts_matrix)[1] #JP: easer to compare data sets when ns are the same
new_lib_size <- sample(1000:10000, size = n, replace = TRUE) # New library sizes
# JP: Maybe the new library sizes should be similar to the original ones?

# Model group effects with individual relative abundances
groups <- c(rep(0:1, each = n / 2), 1)  # JP: + 1, as n is odd
n_taxa <- ncol(counts_matrix)
beta <- rep(0, n_taxa)  # Effect size per taxon
beta[1:10] <- 0.0  # JP: changed from 0.1 to 0.2
Y <- matrix(groups, ncol = 1)  # Binary covariate (0 = control, 1 = case)
beta_all <- matrix(beta, nrow = 1)

# Initialize individual relative abundances
new_individual_rel_abund <- matrix(midasim_setup_para$mean.rel.abund, nrow = n, ncol = n_taxa, byrow = TRUE)
new_individual_rel_abund <- exp(Y %*% beta_all) * new_individual_rel_abund  # Apply group effect
new_individual_rel_abund <- new_individual_rel_abund / rowSums(new_individual_rel_abund, na.rm = TRUE)  # Normalize

# Modify setup with new library sizes and individual relative abundances
midasim_modified_para <- MIDASim.modify(
  midasim_setup_para,
  lib.size = new_lib_size,
  individual.rel.abund = new_individual_rel_abund
)

simulated.data_para <- MIDASim(midasim_modified_para)
summary(simulated.data_para)

# Perform Differential Abundance Analysis (DAA)----------------------------------------
# Transform counts to relative abundances (use sim_rel directly)
rel_abs <- simulated.data_para$sim_rel

# Remove taxa < 3 non-zero counts
rel_abundances <- rel_abs[, colSums(rel_abs > 0, na.rm = TRUE) >= 3]

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

# Prepare metadata
meta <- data.frame(group = factor(groups, labels = c("control", "case")))

# Run ORM for all taxa
res <- rel_abundances |> 
  as.data.frame() |> 
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')

# Compare the results to true effects
true_effects <- tibble(taxon = colnames(simulated.data_para$sim_rel), true_effect = beta)

results <- res |> 
  left_join(true_effects, by = "taxon") |> 
  mutate(q = p.adjust(p_value, method = "fdr"),
         prevalence = colMeans(rel_abundances > 0, na.rm = TRUE),
         significant = q < .05,
         true_positive = sign(estimate) == sign(true_effect) & significant,
         true_negative = true_effect == 0 & !significant,
         false_positive = sign(estimate) != sign(true_effect) & significant,
         false_negative = true_effect != 0 & !significant)

results |> summarize(n_taxa = n(),
                     n_true_effects = sum(true_effect != 0),
                     n_no_true_effects = sum(true_effect == 0),
                     n_significant = sum(significant),
                     n_true_positive = sum(true_positive),
                     n_true_negative = sum(true_negative),
                     n_false_positive = sum(false_positive),
                     n_false_negative = sum(false_negative),
                     fdr = n_false_positive / n_significant,
                     power = n_true_positive / n_true_effects)

# JP:
# To evaluate the performance of methods, you basically need a metric for error
# rate and for statistical power (how effective the method is at detecting true
# effects).

# The standard metric for error rate in the context of microbiome data analysis
# is False Discovery Rate (FDR, the proportion of false positives among all
# significant results). It is a bit problematic to estimate it when there are no
# significant results, but when you calculate the FDR over many simulated data
# sets, with e.g 1/3, 2/9 and 0/0 false positives/significant results, estimate
# the FDR as (1 + 2 + 0) / (3 + 9 + 0) = 3/12. Then you should avoid the "0/0"
# problem.

# I think the standard metric for statistical power is (statistical power) the
# proportion of true effects that are detected as significant.



# JP: Comparing original and simulated data ------------------------------------
library(vegan)
library(ape)

d_original <- counts_matrix / rowSums(counts_matrix)
# d_simulated <- simulated.data_nonpara$sim_rel
d_simulated <- simulated.data_para$sim_rel
d_both <- rbind(d_original, d_simulated)

# Compare sparsities
mean(d_original == 0)
mean(d_simulated == 0)

# Calculate Bray-Curtis dissimilarities/distances
dist <- vegdist(d_both, method = "bray")
# dist <- vegdist(d_both, method = "jaccard", binary = T) #For presence/absence data

# Calculate and plot PCoA (principal coordinate analysis) values
pcoa <- pcoa(dist, correction = "none")
relative_eigenvalues <- pcoa$values$Relative_eig[1:2] #Var explained by 1st and 2nd PCo

# Extract scores of the 1st and 2nd PCoA axes. You can examine other axes, too,
# if they exlain a large amount of variance
data_pcoa <- pcoa$vectors[, 1:2] |>  # 1st and 2nd PCoA axes
  as.data.frame() |> 
  dplyr::rename(PCoA1 = 1, PCoA2 = 2) |> 
  mutate(data = c(rep('Original', nrow(d_original)),
                  rep('Simulated', nrow(d_simulated))))

ggplot(data_pcoa, aes(x = PCoA1, y = PCoA2, color = data)) +
  geom_point(alpha = .5) +
  labs(x = paste0("PCoA1 (", round(relative_eigenvalues[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(relative_eigenvalues[2] * 100, 2), "%)")) +
  facet_wrap(~ data) +
  theme_light()

# The non-parametric simulations seems to produce very similar data to the
# original, but the data from parametric simulations seem to have less variation.