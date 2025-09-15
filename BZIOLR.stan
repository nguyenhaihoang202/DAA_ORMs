// ============================================================================
  // Bayesian zero-inflated ordered logistic model - revised 
// Key changes:
  //   1) Per-taxon ordered cutpoints c_taxa[m], hierarchically shrunk to a template.
//   2) Zero inflation only as baseline per-taxon structural-zero prob (no group slope).
//   3) Tighter global shrinkage for beta: beta_mu ~ N(0, 0.3), beta_sigma ~ half-N(0.3).
//   4) Soft compositional regularization: mean(beta) ~ N(0, 0.3).
// ============================================================================
  
  functions {
    // Partial log-lik over a slice for reduce_sum
    real partial_log_lik(array[] int y_slice, int start, int end,
                         array[] int taxon_idx, array[] int group,
                         array[] vector c_taxa, vector beta, vector phi) {
      real lp = 0;
      for (i in 1:(end - start + 1)) {
        int n = start + i - 1;
        int m = taxon_idx[n];
        real g   = group[n] - 0.5;              // 0 -> -0.5, 1 -> +0.5
        real eta = beta[m] * g;
        vector[rows(c_taxa[m])] c = c_taxa[m];
        real phi_m = phi[m];
        
        if (y_slice[i] == 1) {
          // Mixture for the lowest category
          lp += log_sum_exp(
            bernoulli_lpmf(1 | phi_m),                                   // structural zero
            bernoulli_lpmf(0 | phi_m) + ordered_logistic_lpmf(1 | eta, c) // non-structural -> cat 1
          );
        } else {
          // Positive category must be non-structural
          lp += bernoulli_lpmf(0 | phi_m)
          + ordered_logistic_lpmf(y_slice[i] | eta, c);
        }
      }
      return lp;
    }
  }

data {
  int<lower=1> MN;                           // total observations
  int<lower=1> M;                            // taxa
  int<lower=2> K;                            // ordinal categories
  array[MN] int<lower=1, upper=K> y;         // responses
  array[MN] int<lower=0, upper=1> group;     // 0 control, 1 case
  array[MN] int<lower=1, upper=M> taxon_idx; // taxon id
}

transformed data {
  // Evenly spaced cutpoint template on logit scale
  vector[K-1] c_template;
  for (k in 1:(K-1)) {
    real p = (k * 1.0) / K;
    c_template[k] = logit(p);
  }
}

parameters {
  // Hierarchical taxon effects on ordinal part
  real beta_mu;
  real<lower=0> beta_sigma;
  vector[M] beta_raw;
  
  // Zero-inflation intercept per taxon (no group slope)
  vector[M] a_phi;
  
  // Per-taxon ordered cutpoints with hierarchical shrinkage
  array[M] ordered[K-1] c_taxa;
  real<lower=0> sigma_c; // shared cutpoint deviation scale
}

transformed parameters {
  vector[M] beta = beta_mu + beta_sigma * beta_raw;
  vector[M] phi  = inv_logit(a_phi); // structural-zero probabilities
}

model {
  // Priors for beta: tighter pooling as suggested
  beta_mu    ~ normal(0, 0.3);
  beta_sigma ~ normal(0, 0.3);     // half-normal via <lower=0>
    beta_raw   ~ normal(0, 1);
  
  // Soft compositional regularization
  target += normal_lpdf(mean(beta) | 0, 0.3);
  
  // Zero inflation prior: favor modest structural-zero rates, still flexible
  a_phi ~ normal(-2, 1.0);         // logit^-1(-2) ≈ 0.12 baseline
  
  // Cutpoint hierarchy
  sigma_c ~ normal(0, 0.5);        // half-normal
  for (m in 1:M)
    c_taxa[m] ~ normal(c_template, sigma_c);
  
  // Likelihood via parallel reduce_sum
  target += reduce_sum(
    partial_log_lik, y, 1000,
    taxon_idx, group, c_taxa, beta, phi
  );
}

generated quantities {
  // Centered beta for reporting
  vector[M] beta_centered_mean;
  real beta_center_mean = mean(beta);
  for (m in 1:M) beta_centered_mean[m] = beta[m] - beta_center_mean;
}
