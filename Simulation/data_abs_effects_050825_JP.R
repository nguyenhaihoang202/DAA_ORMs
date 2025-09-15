library(tidyverse)

# The data includes measured _absolute_ abundances of genera.
# The data was used by Vandeputte 2017 and Vieira-Silva 2019
# (https://www.nature.com/articles/s41564-019-0483-9)
load('data_vieira_silva_050825.rds')

d |> group_by(Diagnosis) |> summarize(n = n())

# Observed absolute abundance effects between healthy subjects and those with
# Crohn's disease. This data was used to demonstrate the effects of
# compositionality on DAA https://www.nature.com/articles/nature24460
# The absolute DA effects in this data set are thus realistic but probably close
# to as strong as there exist in human gut microbiome studies. 

abs_effects <- d |> 
  filter(Diagnosis %in% c('mHC', 'CD')) |> 
  mutate(group = ifelse(Diagnosis == 'mHC', 'healthy', 'CD')) |>
  select(group, Abiotrophia:unclassified_Xanthomonadaceae) |> 
  pivot_longer(cols = -group, names_to = 'genus') |> 
  group_by(genus, group) |> 
  mutate(prevalence = sum(value > 0)) |> 
  ungroup() |> 
  group_by(genus) |> 
  mutate(min_prevalence = min(prevalence)) |> 
  filter(min_prevalence >= 5) |> 
  group_by(genus, group) |>
  summarise(mean_abundance = mean(value)) |> 
  ungroup() |> 
  pivot_wider(names_from = group, values_from = mean_abundance) |> 
  mutate(abs_effect = log(CD /healthy))

ggplot(abs_effects, aes(x = abs_effect)) +
  geom_histogram(bins = 50) +
  labs(x = 'Absolute effect size (log(CD / healthy))',
       y = 'Number of genera')

abs_effects |> summarize(m = mean(abs_effect),
                         md = median(abs_effect),
                         sd = sd(abs_effect))

# The observed effects (log ratio of artithmetic means) are roughly normally
# distributed, with the mean of about -1 and SD about 1.1. Note however that
# these are _observed_ effects (in a sample with N = 65 + 28), not the true
# effects (of a population). Therefore, these observed effects also include
# lots of random variation. The true effects thus likely have smaller SD.