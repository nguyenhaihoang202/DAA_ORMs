// ============================================================================
// Bayesian ordinal logistic model for microbiome DAA
// ============================================================================

functions {
  // Asymmetric Laplace log-PDF for beta | tau, nu (mu = 0).
  // tau > 0 (scale), nu in (0,1) (asymmetry; nu=0.5 is symmetric Laplace).
  real asym_laplace_lpdf(real x, real tau, real nu) {
    return log(nu * (1 - nu)) - log(tau)
         - (nu / tau) * fmax(x, 0)
         - ((1 - nu) / tau) * fmax(-x, 0);
  }

  // Parallelizable partial log-likelihood over a slice of observations.
  // Likelihood:
  //   Ordered logistic. Group is centered as g = group - 0.5.
  real partial_log_lik(array[] int y_slice, int start, int end,
                       array[] int taxon_idx, array[] int group,
                       array[] vector c_taxon, vector beta) {
    real lp = 0;
    int L = end - start + 1;
    for (i in 1:L) {
      int n = start + i - 1;                 // row in the full dataset
      int m = taxon_idx[n];                  // taxon for this row
      real g   = group[n] - 0.5;             // centered group: control -0.5, case +0.5
      real eta = beta[m] * g;                // ordinal linear predictor
      vector[rows(c_taxon[m])] c_adj = c_taxon[m];
      lp += ordered_logistic_lpmf(y_slice[i] | eta, c_adj);
    }
    return lp;
  }
}

data {
  int<lower=1> MN;                           // total number of observations (N * M)
  int<lower=1> M;                            // number of taxa
  int<lower=2> K;                            // number of ordinal categories (e.g., 4)
  array[MN] int<lower=1, upper=K> y;         // ordinal outcome per observation
  array[MN] int<lower=0, upper=1> group;     // group indicator: 0 control, 1 case
  array[MN] int<lower=1, upper=M> taxon_idx; // taxon id per observation
}

transformed data {
  // Weak prior center for cutpoints: evenly spaced on logit scale
  vector[K-1] c_prior_mean;
  for (k in 1:(K-1)) {
    real p = (k * 1.0) / K;
    c_prior_mean[k] = logit(p);
  }
}

parameters {
  // Taxon-specific group effects on abundance ranks
  vector[M] beta;

  // Global prior hyperparameters for asymmetric Laplace on beta
  real<lower=0> tau;                 // global scale
  real<lower=0, upper=1> nu;         // global asymmetry (nu=0.5 symmetric)

  // Per-taxon ordered cutpoints (no global base/scale)
  array[M] ordered[K-1] c_taxon;
}

model {
  // ----- Priors -----

  // Global hyperpriors (weakly-informative)
  // Half-Normal on tau (scale on log-odds)
  target += normal_lpdf(tau | 0, 1);
  // Laplace prior on nu centered at 0.5 with tight scale
  target += double_exponential_lpdf(nu | 0.5, 0.05);

  // Asymmetric Laplace prior on beta 
  for (m in 1:M)
    target += asym_laplace_lpdf(beta[m] | tau, nu);

  // Cutpoint priors: taxon-specific, weakly-informative around common template
  // Each taxon has its own ordered[K-1] vector.
  for (m in 1:M)
    c_taxon[m] ~ normal(c_prior_mean, 1.0);

  // ----- Likelihood -----
  target += reduce_sum(
             partial_log_lik, y, 1000,              // grain-size = 1000
             taxon_idx, group, c_taxon, beta
           );
}

generated quantities {
  // Store mean-centered beta for easier cross-taxon comparison
  vector[M] beta_centered_mean;
  real beta_center_mean;
  beta_center_mean = mean(beta);
  for (m in 1:M) beta_centered_mean[m] = beta[m] - beta_center_mean;
}

