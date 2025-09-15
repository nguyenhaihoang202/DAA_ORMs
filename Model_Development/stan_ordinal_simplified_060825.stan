data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1, upper=K> y[N];
  real intercept[N];
  real group[N];
  vector[K - 1] c;
}

parameters {
  real alpha;
  real beta;
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);

  for (n in 1:N)
    y[n] ~ ordered_logistic(alpha * intercept[n] + beta * group[n], c);
}
