//-------------------------------------------------------------------
//  To give some kind of freedom for the intercepts , we do Ordinal logistic with an extra:
//
//   λ  = global stretch → expands or compresses spacing
//
//  They fix most bias from using a rigid cut-point template
//  while adding only a scalar parameters.
//-------------------------------------------------------------------
functions {
  real partial_log_lik(array[] int y_slice, int start, int end,
                       array[] real group_0,
                       array[] real group_1,
                       vector c,
                       real beta_0,
                       real beta_1,
                       real lambda) {
    real lp = 0;
    for (i in 1:(end - start + 1)) {
      int n = start + i - 1;
      real eta = beta_0 * group_0[n] + beta_1 * group_1[n];
      vector[rows(c)] c_adj = lambda * c;
      lp += ordered_logistic_lpmf(y_slice[i] | eta, c_adj);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> y;
  array[N] real group_0;
  array[N] real group_1;
  vector[K - 1] c;  // base cutpoints (e.g. qlogis of cumulative probs)
}

parameters {
  real beta_0;
  real beta_1;
  real<lower=0.1, upper=3> lambda;  // Restricts lambda to prevent dilution
}

model {
  beta_0 ~ normal(0, 3);
  beta_1 ~ normal(0, 3);
  lambda ~ lognormal(0, 0.1); // Centered on 1, reasonable spread

  target += reduce_sum(partial_log_lik, y, 50,
                       group_0, group_1, c,
                       beta_0, beta_1, lambda);
}
