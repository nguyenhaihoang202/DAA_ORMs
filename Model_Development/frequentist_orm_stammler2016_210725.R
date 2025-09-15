library(tidyverse)
library(TreeSummarizedExperiment)

# Function to run ORM without error
run_orm <- function(abundance, metadata, formula){
  
  # Create the design matrix of the model
  mm <- model.matrix(formula, metadata) |> 
    cbind(abundance) |> 
    tibble::as_tibble() |> 
    dplyr::select(-"(Intercept)")
  
  # Get the indices and names of the variables in the model matrix
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  # Initialize result frame in case of full failure
  res <- data.frame(variable = vars, estimate = NA, se = NA, p_value = NA)
  
  # Fit ordinal regression model.
  # Intercepts are automatically added and the abundance is automatically
  # transformed to ordinal ranks.
  # Fit ordinal regression model with tryCatch and increased iterations
  # Increases 'maxiter = 100' in 'rms::orm' calls to enhance convergence chances, unlike the original's default (20).
  fit_1 <- try(rms::orm(abundance ~ ., data = mm, maxiter = 100), silent = TRUE) 
  
  if (inherits(fit_1, "try-error")) {
    return(res |> tibble::rownames_to_column("variable"))  # Return NAs on failure
  }
  
  # Extract the score statistic that is used to calculate p-values
  score_1 <- fit_1$stats["Score"]
  
  # Extract estimates and SE
  res <- data.frame(estimate = fit_1$coefficients[vars],
                    se = sqrt(diag(vcov(fit_1))[vars]),
                    p_value = NA)
  
  # Calculate p-values with tryCatch for null models
  if (length(inds) > 1) {
    for (i in inds) {
      fit_0 <- try(rms::orm(abundance ~ ., data = mm[, -i], maxiter = 100), silent = TRUE)
      if (inherits(fit_0, "try-error")) {
        res$p_value[i] <- NA
      } else {
        score_0 <- fit_0$stats["Score"]
        res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
      }
    }
  } else {
    # Null model score is 0
    res$p_value <- as.numeric(1 - pchisq(score_1 - 0, df = 1))
  }
  
  return(res |> tibble::rownames_to_column("variable"))
}


# Load data and prepare data----------------------------------------------------
load('mb_datasets_gamboa_tuz.rds')

# The data are in the in the from of a TreeSummarizedExperiment object
# Extract the Stammler 2016 dataset
tse <- data_mbd_raw$Stammler_2016_16S_spikein

# Filter out taxa with prevalence < 0.01
tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)

# Extract counts from the TreeSummarizedExperiment object
counts <- tse |> assay() |> t()

# Relative abundances
rel_abs <- as.data.frame(counts / rowSums(counts))

# Metadata - create appropriate grouping variable
meta <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(
    # Map the study_condition/timepoints into broad early/late categories:
    group = case_when(
      study_condition %in% c("preASCT", "d0", "d7") ~ "uncompleted",
      study_condition %in% c("d14") ~ "completed",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("uncompleted", "completed"))
  ) |> 
  filter(!is.na(group))

# Spike-in bacteria identification
# The Stammler dataset has spike-in bacteria with "XXXX" suffix
  spike_in_ids <- c("AF323500XXXX", "AB247615XXXX", "AB076660XXXX")  # Salinibacter ruber, Rhizobium radiobacter, Alicyclobacillus acidiphilus

# Create ground truth for spike-in bacteria
bvs <- rowData(tse) |> 
  as.data.frame() |> 
  tibble::rownames_to_column('taxon') |> 
  mutate(ground_truth = case_when(
    taxon %in% spike_in_ids ~ 'spike_in',
    TRUE ~ 'endogenous'
  )) |>
  select(taxon, ground_truth)

# Run ORM analysis---------------------------------------------------------------
res <- rel_abs |> 
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')

sl <- 0.10 # Significance level

# Compare results against ground truth
res2 <- res |> 
  left_join(bvs, by = 'taxon') |> 
  mutate(
    significant = p_value < sl
  )

# Spike-in bacteria results
spike_in_details <- res2 |> 
  filter(ground_truth == 'spike_in') |> 
  select(taxon, estimate, p_value, significant)

spike_in_details
# These results are expected: none of the spike-in bacteria are significant,
# indicating proper normalization. Since spike-ins were added at equal concentration to all samples,
# they should not differ between groups; non-significant p-values here confirm that
# the workflow is correctly controlling for technical variation.

# JP: Ok, I think that is a very good idea to test whether the spike-in bacteria
# are significantly different between some samples/groups! Altohugh, in this 
# small dataset, even large differences would not easily become significant.
# There is, however, again the problem that the samples are not independent, but
# there are several samples from the same subjects. But maybe it is still fine
# to use this data. Maybe you can try to compare (all) samples
# from Patients 1 and 2 vs. Patients 3, 4 and 5. Then groups would be more 
# balanced (8 vs. 9 samples), and thus the sample size would effectively be a
# bit higher. That would probably also be more challenging for the DAA methods,
# because there is likely lot of variation between the patients' microbiomes.
# (Maybe you could also try to test between patient 1 (4 samples) vs. patient 2
# (4 samples)?)
# 

# From Patients 1 and 2 vs. Patients 3, 4 and 5

# Metadata - create appropriate grouping variable
meta_group_patients <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(
    # Map the subject_id (patients) into groups: Patients 1 & 2 vs. Patients 3, 4 & 5
    group = case_when(
      subject_id %in% c("Patient1", "Patient2") ~ "groupA",
      subject_id %in% c("Patient3", "Patient4", "Patient5") ~ "groupB",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("groupA", "groupB"))
  ) |> 
  filter(!is.na(group))


# Run ORM analysis---------------------------------------------------------------
res_group_patients <- rel_abs |> 
  map(~ run_orm(., metadata = meta_group_patients, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')

sl <- 0.10 # Significance level

# Compare results against ground truth
res2_group_patients <- res_group_patients |> 
  left_join(bvs, by = 'taxon') |> 
  mutate(
    significant = p_value < sl
  )

# Spike-in bacteria results
spike_in_details_group_patients <- res2_group_patients |> 
  filter(ground_truth == 'spike_in') |> 
  select(taxon, estimate, p_value, significant)

spike_in_details_group_patients



#-------------------------------------------------------------------------------
# JP: Filter samples and taxa in the correct order------------------------------
# Patients 1 vs Patients 2

# The data are in the in the from of a TreeSummarizedExperiment object
# Extract the Stammler 2016 dataset
tse_1vs2 <- data_mbd_raw$Stammler_2016_16S_spikein

## Filter the object to include only healthy and BV samples
tse_1vs2 <- tse_1vs2[, tse_1vs2$subject_id %in% c('Patient1', 'Patient2')]

# Filter out taxa with prevalence < 0.01 and no count
tse_1vs2 <- mia::subsetByPrevalent(tse_1vs2, prevalence = 0.01)


# Extract counts from the TreeSummarizedExperiment object
counts_1vs2 <- tse_1vs2 |> assay() |> t()

# Relative abundances
rel_abs_1vs2 <- as.data.frame(counts_1vs2 / rowSums(counts_1vs2))

summary(colSums(rel_abs_1vs2))
# JP: Now all taxa have at least some non-zeros

# Metadata - create appropriate grouping variable
meta_1vs2 <- tse_1vs2 |>
  colData() |>
  as.data.frame() |>
  mutate(
    # Map the subject_id (patients) into groups: Patients 1 & 2 vs. Patients 3, 4 & 5
    group = case_when(
      subject_id %in% c("Patient1") ~ "groupA",
      subject_id %in% c("Patient2") ~ "groupB",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("groupA", "groupB"))
  ) |>
  filter(!is.na(group))

# Create ground truth for spike-in bacteria
bvs_1vs2 <- rowData(tse_1vs2) |>
  as.data.frame() |>
  tibble::rownames_to_column('taxon') |>
  mutate(ground_truth = case_when(
    taxon %in% spike_in_ids ~ 'spike_in',
    TRUE ~ 'endogenous'
  )) |>
  select(taxon, ground_truth)


# Run ORM analysis---------------------------------------------------------------
res_1vs2 <- rel_abs_1vs2 |>
  map(~ run_orm(., metadata = meta_1vs2, formula = ~ group)) |>
  bind_rows(.id = 'taxon')

# JP: It works now

sl <- 0.05 # Significance level

# Compare results against ground truth
res2_1vs2 <- res_1vs2 |>
  left_join(bvs_1vs2, by = 'taxon') |>
  mutate(,
    significant = p_value < sl
    # JP: As you are only intersted in three taxa, there
    # is on need for multiplicity correction. Define
    # significance thus with the undajusted p-value. 
  )

# Spike-in bacteria results
spike_in_details_1vs2 <- res2_1vs2 |>
  filter(ground_truth == 'spike_in') |>
  select(taxon, estimate, p_value, significant)

spike_in_details_1vs2

# Only 3 taxa, more sensible to use the p-value.

