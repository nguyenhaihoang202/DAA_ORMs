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

# Extract HMP_2012_16S_gingival_V13 dataset
tse <- data_mbd_raw$HMP_2012_16S_gingival_V13

# Filter to relevant body subsites
tse <- tse[, tse$body_subsite %in% c('subgingival_plaque', 'supragingival_plaque')]

# Filter taxa with prevalence < 0.05, increased due to bigger dataset
# JP: Ok, that filtering is good.

# tse <- mia::subsetByPrevalent(tse, prevalence = 0.05)
# It may be more sensible to agglomerate to genus level. As far as I know, the
# 16S sequencing is not very reliable belowe genus level (but I'm not sure) 
tse <- mia::agglomerateByPrevalence(tse, rank = 'genus', prevalence = 0.05)

# Extract counts
counts <- tse |> assay() |> t()

# Relative abundances 
rel_abs <- as.data.frame(counts / rowSums(counts))

# Metadata with corrected factor levels 
meta <- tse |> 
  colData() |> 
  as.data.frame() |> 
  mutate(group = factor(body_subsite,
                        levels = c('subgingival_plaque', 'supragingival_plaque'),
                        labels = c('subgingival', 'supragingival')))

# Define ground truth (filter to aerobic/anaerobic annotations)
taxonomy_df <- rowData(tse) |> as.data.frame() |> 
  tibble::rownames_to_column("taxon")
taxonomy_df <- taxonomy_df %>% filter(taxon_annotation %in% c("aerobic", "anaerobic"))

# JP: If you run this table on OTU level data, you see that aerobic/anaerobic
# classification is made only on genus level.
table(taxonomy_df$genus, taxonomy_df$taxon_annotation)

bvs <- taxonomy_df %>%
  mutate(
    ground_truth = case_when(
      taxon_annotation == "aerobic"   ~ "supragingival",
      taxon_annotation == "anaerobic" ~ "subgingival",
      TRUE                            ~ "none"
    )
  ) %>%
  select(taxon, ground_truth)

# Run DAA-----------------------------------------------------------------------
res <- rel_abs |>
  map(~ run_orm(., metadata = meta, formula = ~ group)) |> 
  bind_rows(.id = 'taxon')

sl <- 0.1 # Significance level

# Compare results against ground truth
# Positive estimate: more abundant in supragingival 
res2 <- res |> 
  dplyr::left_join(bvs, by = 'taxon') |> 
  tidyr::drop_na() |>  # Remove all rows with any NA values after the join
                       # JP: I would not make this filtering. At least not 
                       #     before calculating the q-values.
  mutate(
    q = p.adjust(p_value, method = 'BH'),
    res = case_when(
      q < sl & estimate > 0 ~ 'supragingival',
      q < sl & estimate < 0 ~ 'subgingival',
      T ~ 'ns'
    ),
    correct = ifelse(res == ground_truth, T, F),
    incorrect = ifelse(res == 'supragingival' & ground_truth == 'subgingival' |
                         res == 'subgingival' & ground_truth == 'supragingival', T, F)
  )

# Summarize the results against the ground truth
res2 |> 
  summarize(
    correct = sum(correct),
    incorrect = sum(incorrect),
    n = sum(!is.na(ground_truth))  # Count non-NA values instead of non-"none"
  )
# 455 correctly identified taxa, 30 incorrectly identified taxa among 2422 taxa with known ground truth (under significance level).
# Low performance?

# JP: I don't think that is bad performance. It is normal that everything is not
#     found. But difficult to say without comparing it to another method.

# JP: On genus level, the result is 22 vs. 1 (out of 31) 