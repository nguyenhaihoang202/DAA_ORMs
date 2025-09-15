library(tidyverse)
library(ggpubr)
library(brms)
options(mc.cores = parallel::detectCores())


# Here I show how to create ideal data for ordinal logistic regresssion---------

# Example of an ordinal variable with five classes and the standard logistic
# distribution as the latent variable.
# The thresholds of the classes are at -3.0, -1.0, 0.5 and 1.5

d1 <- tibble(x = seq(-5, 5, .1),
             density = dlogis(x, location = 0, scale = 1))

ggplot(d1, aes(x = x, y = density)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-3, -1, 0.5, 1.5), color = 'red') +
  geom_line() +
  theme_light()

# The probability to get value from class 1 is now
(p1 <- plogis(-3)) # P(x < -3)

# The probability to get value from class 2 is now
(p2 <- plogis(-1) - plogis(-3)) # P(-3 < x < -1)

# The probability to get value from class 3 is now
(p3 <- plogis(0.5) - plogis(-1)) # P(-1.0 < x < 0.5)

# etc...
(p4 <- plogis(1.5) - plogis(0.5))
(p5 <- 1 - plogis(1.5))

# The values of the ordinal variable can now be generated as follows
# (sampe size n = 100)
n <- 100
set.seed(1)
y <- sample(1:5, size = n, replace = TRUE, prob = c(p1, p2, p3, p4, p5))

c(p1, p2, p3, p4, p5) # The true probabilities of the classes
prop.table(table(y)) # The observed proportions of the classes


# Run Bayesian ordinal regresssion model for the data
brm1 <- brm(y ~ 1, data = tibble(y = y),
            family = cumulative("logit"),
            iter = 2000,
            warmup = 1000,
            chains = 4,
            seed = 1)

# Print the estimates from the model. You can see that the estimates of the
# thresholds (Intercepts) are close to the true values of -3, -1, 0.5 and 1.5
print(brm1)


# Two group model---------------------------------------------------------------

# Assume there are two groups A and B. In the group A the probabilities are as
# above, but in group B the probabilities are "shifted" so that smaller classes
# have higher probabillities. The log odds ratio is 0.8.

plot_A <- ggplot(d1, aes(x = x, y = density)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-3, -1, 0.5, 1.5), color = 'red') +
  geom_line() +
  labs(title = 'Group A') +
  theme_light()

plot_B <- ggplot(d1, aes(x = x, y = density)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-3, -1, 0.5, 1.5) + 0.8, color = 'red') +
  geom_line() +
  labs(title = 'Group B') +
  theme_light()

ggarrange(plot_A, plot_B, nrow = 2)

# The true probabilities of the classes in group B are now
(p1_B <- plogis(-3 + 0.8)) # P(x < -3 + 0.8)
(p2_B <- plogis(-1 + 0.8) - plogis(-3 + 0.8)) # P(-3 < x < -1 + 0.8)
(p3_B <- plogis(0.5 + 0.8) - plogis(-1 + 0.8)) # P(-1.0 < x < 0.5 + 0.8)
(p4_B <- plogis(1.5 + 0.8) - plogis(0.5 + 0.8)) # P(0.5 < x < 1.5 + 0.8)
(p5_B <- 1 - plogis(1.5 + 0.8)) # P(x > 1.5 + 0.8)


# We can now create an example data set with two groups A and B
set.seed(1)
y2 <- c(sample(1:5, size = n, replace = TRUE,
               prob = c(p1, p2, p3, p4, p5)),
        sample(1:5, size = n, replace = TRUE,
               prob = c(p1_B, p2_B, p3_B, p4_B, p5_B)))

group <- c(rep('A', n), rep('B', n))


# Run Bayesian ordinal regresssion model for the two group data
d2 <- tibble(y = y2, group = group)
brm2 <- brm(y ~ group, data = d2,
            family = cumulative("logit"),
            iter = 2000,
            warmup = 1000,
            chains = 4,
            seed = 1)

# Print the estimates from the model. You can see that the estimates of the
# thresholds (Intercepts) in the baseline group level A are somwhat close to the
# true values of -3, -1, 0.5 and 1.5. The true log OR between the group (0.8) is 
# also quite well estimated
print(brm2)



# Analyzing microbial relative abundances with ORM------------------------------

# The relative abundances of microbes have typically many zeros (the "semi"
# part) and then values between 0 and 1 (the "continuous" part), so they maybe
# can be considered as a "semi-continuous" variable.

# Here is an example dataset from subjects with colorectal cancer (CRC, group =
# "case") and from healthy controls (group = "control"). It also includes 
# relative abundances of 250 microbial taxa.

load("example_mb_data.rds")

# Take the relative abudances of taxon 3 as an example
ggplot(mb_data, aes(x = taxon_3)) +
  geom_histogram(bins = 500) +
  facet_grid(group ~ .)

# For ordinal regression the relative abundances have to be first transformed to 
# categories (integera) based on their ordinal value.
taxon_3_ord <- factor(mb_data$taxon_3,
                      levels = sort(unique(mb_data$taxon_3)),
                      labels = 1:length(unique(mb_data$taxon_3))) |> 
  as.integer()

# Check the distribution of the ordinal values
mb_data$taxon_3
taxon_3_ord


# Run Bayesian ordinal regression model to compare relative abundances between
# case and control groups
mbd <- tibble(group = mb_data$group,
              y = taxon_3_ord)

brm3 <- brm(y ~ group,
            data = mbd,
            family = cumulative("logit"),
            iter = 2000,
            warmup = 1000,
            chains = 4,
            seed = 1)

# Print the results
print(brm3)

# There seems to no real difference between the groups
brms::posterior_summary(brm3, variable = "b_groupcase")









