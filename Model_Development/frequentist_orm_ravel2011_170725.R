library(tidyverse)
library(TreeSummarizedExperiment)

# Function to run ORM
run_orm <- function(abundance, metadata, formula){
  
  # Create the design matrix of the model
  mm <- model.matrix(formula, metadata) |> 
    cbind(abundance) |> 
    tibble::as_tibble() |> 
    dplyr::select(-"(Intercept)")
  
  # Get the indices and names of the variables in the model matrix
  inds <- 1:(ncol(mm) - 1)
  vars <- colnames(mm)[inds]
  
  # Fit ordinal regression model.
  # Intercepts are automatically added and the abundance is automatically
  # transformed to ordinal ranks.
  fit_1 <- rms::orm(abundance ~ ., data = mm)
  
  # Extract the score statistic that is used to calculate p-values
  score_1 <- fit_1$stats["Score"]
  
  # Extract the log odds estimate and its standard error for each variable
  res <- data.frame(estimate = fit_1$coefficients[vars],
                    se = sqrt(diag(vcov(fit_1))[vars]),
                    p_value = NA)
  
  # Calculate the p-value based on score test for each variable
  if(length(inds) > 1){
    
    for(i in inds){
      
      #Fit the "null model" (the model without the variable of interest)
      fit_0 <- rms::orm(abundance ~ ., data = mm[, -i])
      score_0 <- fit_0$stats["Score"]
      
      # p-value is based on the difference of score statistics. Under the null
      # hypothesis it follows the chi squared distribution with 1 degree of
      # freedom. 
      res$p_value[i] <- as.numeric(1 - pchisq(score_1 - score_0, df = 1))
    }
    
  }else{
    
    # If there is only one variable in the model, the "null model" would include
    # only intercepts. The score statistic for such null model would be zero.
    res$p_value <- as.numeric(1 - pchisq(score_1 - 0, df = 1))
  }
  
  return(res |> tibble::rownames_to_column("variable"))
}


# Load data and prepare data----------------------------------------------------
load('mb_datasets_gamboa_tuz.rds')

# The data are in the in the from of a TreeSummarizedExperiment object
# Extract the Ravel 2011 BV dataset
tse <- data_mbd_raw$Ravel_2011_16S_BV

# Filter the object to include only healthy and BV samples
tse <- tse[, tse$study_condition %in% c('healthy', 'bacterial_vaginosis')]

# Filter out taxa with prevalence < 0.0001
# tse <- mia::subsetByPrevalent(tse, prevalence = 0.01)
# JP: Maybe you could use here, too, 20% filtering to have the same filtering
# as in the Gamboa-Tuz paper.
tse <- mia::subsetByPrevalent(tse, prevalence = 0.20)

# In the G-T paper thay say that they summarized at genus level. I'm not sure
# what that means. At least the code below does not work somehow.. it looses
# almost all the counts..
# tse <- mia::agglomerateByRank(tse, rank = 'genus', prevalence = .20)

#Extract counts from the TreeSummarizedExperiment object
counts <- tse |> assay() |> t()

#Relative abundances
rel_abs <- as.data.frame(counts / rowSums(counts))

# Metadata
meta <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(group = factor(study_condition,
                        levels = c('healthy', 'bacterial_vaginosis'),
                        labels = c('healthy', 'bv')))


# BV signatures. These define the "ground truth".
# BV associated bacteria should be more abundant in women with bacterial vaginosis
# HV associated bacteria should be more abundant in healthy women
bvs <- rowData(tse) |> 
  as.data.frame() |> 
  rownames_to_column('taxon') |> 
  mutate(ground_truth = case_when(taxon_annotation == 'bv-associated' ~ 'bv',
                                  taxon_annotation == 'hv-associated' ~ 'healthy',
                        T ~ 'none')) |>
  select(taxon, ground_truth)


# Run DAA-----------------------------------------------------------------------

res <- rel_abs |> 
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')


sl <- .10 # Significance level

# Compare the results against the ground truth
# Positive estimate indicates that the taxon is more abundant in women with BV   
res2 <- res |> 
  left_join(bvs, by = 'taxon') |> 
  mutate(q = p.adjust(p_value, method = 'BH'),
         res = case_when(q < sl & estimate > 0 ~ 'bv',
                         q < sl & estimate < 0 ~ 'healthy',
                         T ~ 'ns'),
         correct = ifelse(res == ground_truth, T, F),
         incorrect = ifelse(res == 'bv' & ground_truth == 'healthy' |
                      res == 'healthy' & ground_truth == 'bv', T, F))

# Summarize the results against the ground truth
res2 |> summarize(correct = sum(correct),
                  incorrect = sum(incorrect),
                  n = sum(ground_truth != 'none'))

# 20 correctly identified taxa among 28 taxa with known ground truth
# No incorrect findings

# JP: When using 20% filtering, the result is 14 vs. 0 (out of 16)
# 14 seems to be more than for Wilcox.TSS (10) in Fig 2a) in G-T paper.
# Maybe it is due to not agglomerating to genus level? But I'm not sure..
