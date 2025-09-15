library(tidyverse)
library(readxl)


#Means and SDs from Vandeputte data---------------------------------------------
counts <- read_excel("DataS1.xlsx", sheet = "S1-6") %>% 
  rename(sid = '...1') %>% 
  mutate_all(~ str_replace(., 'NA', '')) %>% 
  mutate_all(as.numeric) %>% 
  mutate(reads = rowSums(select(., -sid), na.rm = T))

meta <-  read_excel("DataS1.xlsx", sheet = "S1-3") %>% 
  select(sid = '...1',
         day = Day_Number,
         id = ID_Number,
         load = Cell_count_per_gram) %>% 
  mutate_at(vars(sid, day, id, load), ~ as.numeric(str_replace(., 'NA', '')))

abundances <- meta %>% 
  left_join(., counts, 'sid') %>% 
  mutate(sf = reads / load,
         id = as_factor(id)) %>% 
  filter(!is.na(load) &
           reads > 10 ^ 3 & reads < 10 ^ 5 &
             load > 10 ^ 10 & day < 45) %>% 
  select(sf, contains('_')) %>% 
  mutate_at(vars(contains('_')), ~ . / sf) %>% 
  select(-sf) %>% 
  as.matrix()


mean_vars <- cbind(mean = colMeans(abundances),
                   var = Rfast::colVars(abundances)) %>% 
  as.data.frame() |> 
  filter(mean > 0)

taxa <- rownames(mean_vars)

cors0 <- cor(abundances, method = 'spearman')
cors <- cors0[taxa, taxa]


# Save the results--------------------------------------------------------------
save(cors, file = 'cors_vandeputte_250625.rds')
save(mean_vars, file = 'mean_vars_vandeputte_250625.rds')
