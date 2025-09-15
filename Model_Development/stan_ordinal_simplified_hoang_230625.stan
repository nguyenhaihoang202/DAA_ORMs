functions {
   /**
   * partial_log_lik - Computes part of the log-likelihood for a chunk of data.
   * This function is called in parallel using reduce_sum() to enable multi-threading.
   */
  real partial_log_lik(array[] int y_slice, int start, int end,
                       array[] real group_0,
                       array[] real group_1,
                       vector c,
                       real beta_0,
                       real beta_1) {
    real lp = 0;
    // Loop over the slice (1-based index for y_slice)
    for (i in 1:(end - start + 1)) {
      int n = start + i - 1;
      real eta = beta_0 * group_0[n] + beta_1 * group_1[n];
      lp += ordered_logistic_lpmf(y_slice[i] | eta, c);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=2> K;                              // Must be at least 2
  array[N] int<lower=1, upper=K> y;
  array[N] real group_0;
  array[N] real group_1;
  vector[K - 1] c;                             // K-1 cutpoints for K classes
}

parameters {
  real beta_0;
  real beta_1;
}

model {
  beta_0 ~ normal(0, 3);
  beta_1 ~ normal(0, 3);

  // Use reduce_sum to parallelize likelihood computation across multiple threads
  // The second argument '50' is the grainsize (chunk size), which can be adjusted for better performance (small dataset, larger grainsize)
  target += reduce_sum(partial_log_lik, y, 50, group_0, group_1, c, beta_0, beta_1);
}
