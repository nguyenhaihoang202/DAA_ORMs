// ============================================================================
// Bayesian zero-inflated ordered logistic model for microbiome DAA
// - Presence process: Bernoulli for structural zeros with per-taxon intercept a_phi[m]
//   and a single global group slope b_phi that shifts prevalence for all taxa.
// - Abundance process: ordered logistic (proportional odds) on ordinal categories
//   with taxon-specific group effects beta[m]. Global cutpoints c_base are shared,
//   with per-taxon shift delta_c[m]. 
// - Hierarchical prior on beta enables partial pooling (informed by the overall group data) across taxa.
// - Compositional regularization: softly shrink mean(beta) to 0.
// - Computation: reduce_sum parallelizes the likelihood over observations.
// - Identifiability: Anchors on mean(delta_c) and mean(beta) prevent global shifts; monitor divergences in sparse data.
// ============================================================================

functions {
  // Parallelizable partial log-likelihood over a slice of observations.
  // Arguments:
  //   y_slice    : integer outcomes in {1, ..., K}, sliced for reduce_sum
  //   start, end : indices mapping the slice back to full data
  //   taxon_idx  : taxon index per observation (1..M)
  //   group      : group indicator per observation in {0,1}
  //   c_shifted  : per-taxon adjusted cutpoints (vector[K-1] per taxon)
  //   beta       : per-taxon group effects for the ordinal part
  //   phi_ctrl   : per-taxon P(structural zero) for control
  //   phi_case   : per-taxon P(structural zero) for case
  //
  // Likelihood:
  //   If y == 1 (lowest category), combine two paths:
  //     (i) structural zero from Bernoulli(phi_ni)
  //     (ii) not structural zero times ordinal probability of category 1
  //   If y > 1, outcome must be from the ordinal path and not a structural zero.
  real partial_log_lik(array[] int y_slice, int start, int end,
                       array[] int taxon_idx, array[] int group,
                       array[] vector c_shifted, vector beta,
                       vector phi_ctrl, vector phi_case) {
    real lp = 0;
    for (i in 1:(end - start + 1)) {
      int n = start + i - 1;                       // row in the full dataset
      int m = taxon_idx[n];                        // taxon for this row
      real g   = group[n] - 0.5;                   // centered group: control -0.5, case +0.5
      real eta = beta[m] * g;                      // ordinal linear predictor
      real phi_ni = (group[n] == 1) ? phi_case[m]  // precomputed P(structural zero)
                                    : phi_ctrl[m];
      vector[rows(c_shifted[m])] c_adj = c_shifted[m];  // taxon cutpoints

      if (y_slice[i] == 1) {
        // Zero-inflated ordinal: mixture for the lowest category
        lp += log_sum_exp(
               bernoulli_lpmf(1 | phi_ni),                                  // structural zero
               bernoulli_lpmf(0 | phi_ni) + ordered_logistic_lpmf(1 | eta, c_adj) // non-structural zero -> category 1
             );
      } else {
        // Positive category: must be non-structural
        lp += bernoulli_lpmf(0 | phi_ni)
              + ordered_logistic_lpmf(y_slice[i] | eta, c_adj);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> MN;                           // total number of observations
  int<lower=1> M;                            // number of taxa
  int<lower=2> K;                            // number of ordinal categories
  array[MN] int<lower=1, upper=K> y;         // ordinal outcome per observation
  array[MN] int<lower=0, upper=1> group;     // group indicator: 0 control, 1 case
  array[MN] int<lower=1, upper=M> taxon_idx; // taxon id per observation
}

transformed data {
  // Prior mean for global cutpoints: evenly spaced on logit scale
  vector[K-1] c_mean;
  for (k in 1:(K-1)) {
    real p = (k * 1.0) / K;
    c_mean[k] = logit(p); // weak prior location for c_base
  }
}

parameters {
  // Hierarchical prior for taxon effects on abundance ranks
  real beta_mu;                   // global mean effect
  real<lower=0> beta_sigma;       // global SD of effects
  vector[M] beta_raw;             // standardized taxon effects

  // Zero inflation (prevalence) submodel
  vector[M] a_phi;                // per-taxon intercept for structural zeros
  real      b_phi;                // global group slope for structural zeros

  // Ordered logistic cutpoints: global template plus taxon shifts
  ordered[K-1] c_base;           // global ordered cutpoints
  vector[M]     delta_c;         // per-taxon shift in cutpoints
}

transformed parameters {
  // Taxon effects on abundance ranks with pooling
  vector[M] beta = beta_mu + beta_sigma * beta_raw;

  // Build per-taxon cutpoints: location shift only (no per-taxon scale)
  array[M] vector[K-1] c_shifted;
  for (m in 1:M)
    c_shifted[m] = c_base + delta_c[m];

  // Precompute zero-inflation probabilities for speed & stability
  vector[M] phi_ctrl = inv_logit(a_phi);            // group == 0
  vector[M] phi_case = inv_logit(a_phi + b_phi);    // group == 1
}

model {
  // Priors for taxon effects
  beta_mu ~ normal(0, 3);
  beta_sigma ~ exponential(1);
  beta_raw ~ normal(0, 1);

  // Compositional regularization: softly center the mean effect near 0
  // This discourages a global drift of all taxa in one direction.
  target += normal_lpdf(mean(beta) | 0, 0.3);

  // Zero inflation (prevalence) priors
  a_phi ~ normal(0, 1);      // baseline prevalence per taxon
  b_phi ~ normal(0, 1);      // global prevalence shift for group

  // Cutpoint priors: global template plus modest taxon deviations
  c_base  ~ normal(c_mean, 0.5);
  delta_c ~ normal(0, 0.5);

  // Soft anchor to fix the overall location of thresholds across taxa
  target += normal_lpdf(mean(delta_c) | 0, 0.1);

  // Parallelized likelihood over observations using reduce_sum
  target += reduce_sum(
             partial_log_lik, y, 1000,        // grain-size tuned for MN ~ O(1e4-1e5)
             taxon_idx, group, c_shifted, beta, phi_ctrl, phi_case
           );
}

generated quantities {
  // Store mean-centered beta for easier cross-taxon comparison 
  vector[M] beta_centered_mean;
  real beta_center_mean;
  beta_center_mean = mean(beta);
  for (m in 1:M) beta_centered_mean[m] = beta[m] - beta_center_mean;
}
