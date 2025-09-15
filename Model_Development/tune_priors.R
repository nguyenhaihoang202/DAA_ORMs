# =============================================================================
# Toy data + prior tuning for BALOR
# - Simulate from your ordinal model (K=4)
# - Fit once with current priors (tau ~ half-N(0,1), nu ~ Laplace(0.5, 0.05))
# - Reweight draws to emulate alternative priors for tau and nu
# - Score precision/recall using posterior sign probability (>= 0.95)
# =============================================================================
library(tidyverse)
library(cmdstanr)
library(posterior)
library(stringr)

# ----------------------------- Sim helpers -----------------------------------

# Draw per-taxon ordered cutpoints (K=4) with random location + positive gaps
draw_cutpoints <- function(M, mu1 = 0, sd1 = 1.0, rate_gap = 1.0, K = 4) {
  stopifnot(K == 4)
  c_list <- vector("list", M)
  for (m in 1:M) {
    c1 <- rnorm(1, mu1, sd1)
    d1 <- rexp(1, rate_gap)
    d2 <- rexp(1, rate_gap)
    c_list[[m]] <- c(c1, c1 + d1, c1 + d1 + d2)
  }
  c_list
}

# Simulate sparse, possibly asymmetric taxon effects
simulate_beta <- function(M, prop_signal = 0.10, pos_frac = 0.7, effect_scale = 0.8) {
  beta <- rep(0, M)
  S <- max(1, round(M * prop_signal))
  idx <- sample.int(M, S)
  mags <- rexp(S, rate = 1 / effect_scale)        # heavy-tailed magnitudes
  signs <- ifelse(runif(S) < pos_frac, +1, -1)
  beta[idx] <- signs * mags
  list(beta = beta, signal_idx = sort(idx))
}

# Sample one ordinal outcome from ordered-logistic with K=4
rordlogit1 <- function(eta, c_vec) {
  cdf <- c(plogis(c_vec - eta), 1)
  probs <- c(cdf[1], diff(cdf))
  sample.int(4, 1, prob = probs)
}

# Build toy dataset in the exact Stan layout of BALOR.stan
make_toy <- function(N = 120, M = 150, K = 4,
                     prop_signal = 0.1, pos_frac = 0.6, effect_scale = 0.8,
                     seed = 1) {
  set.seed(seed)
  g_bin <- rbinom(N, 1, 0.5)     # 0 control, 1 case (Stan expects 0/1)
  g     <- g_bin - 0.5           # centered for the linear predictor
  
  bet   <- simulate_beta(M, prop_signal, pos_frac, effect_scale)
  beta_true <- bet$beta
  cut_list  <- draw_cutpoints(M, K = K)
  
  y_vec <- integer(N * M)
  taxon_idx <- rep(1:M, times = N)
  group_vec <- rep(g_bin, each = M)
  k <- 1L
  for (n in 1:N) {
    for (m in 1:M) {
      eta <- beta_true[m] * g[n]
      y_vec[k] <- rordlogit1(eta, cut_list[[m]])
      k <- k + 1L
    }
  }
  
  list(
    stan_data = list(
      MN = N * M, M = M, K = K,
      y = y_vec,
      group = group_vec,
      taxon_idx = taxon_idx
    ),
    truth = list(
      beta_true = beta_true,
      signal_idx = bet$signal_idx,
      cutpoints_true = cut_list,
      g_centered = g
    )
  )
}

# ---------------------- Importance reweighting tools -------------------------

# Log pdfs for priors 
log_half_normal <- function(x, sd) {  # x >= 0
  ifelse(x < 0, -Inf, log(sqrt(2/pi)) - log(sd) - 0.5 * (x/sd)^2)
}
log_laplace <- function(x, mu, b) {   # basic Laplace kernel (truncated to [0,1] in Stan by support)
  -log(2*b) - abs(x - mu)/b
}

# Weighted quantile for summaries 
wtd_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  o <- order(x)
  x <- x[o]; w <- w[o]; w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which.min(abs(cw - p))])
}

# Apply prior reweighting given new tau_sd and nu_scale
reweight_draws <- function(draws_df, tau_sd_new, nu_scale_new,
                           tau_sd_old = 1.0, nu_scale_old = 0.05) {
  tau <- as.numeric(draws_df$tau)
  nu  <- as.numeric(draws_df$nu)
  # compute log-weights as ratio of new to old priors
  lw <- log_half_normal(tau, tau_sd_new) - log_half_normal(tau, tau_sd_old) +
    log_laplace(nu, 0.5, nu_scale_new) - log_laplace(nu, 0.5, nu_scale_old)
  lw <- lw - max(lw)                        # stabilize
  w  <- exp(lw); w <- w / sum(w)
  w
}

# Summarize betas with weights and score FDR/power at sign-prob threshold
summarize_weighted <- function(draws_df, weights, beta_true, thresh = 0.95) {
  beta_cols <- grep("^beta\\[", names(draws_df), value = TRUE)
  res <- map_dfr(seq_along(beta_cols), function(j) {
    b <- as.numeric(draws_df[[beta_cols[j]]])
    probs <- wtd_quantile(b, weights, c(0.025, 0.5, 0.975))
    tibble(
      taxon_id = j,
      est = probs[2],
      lwr95 = probs[1],
      upr95 = probs[3],
      p_gt0 = sum(weights[b > 0]),
      p_lt0 = sum(weights[b < 0])
    )
  })
  df <- res |>
    mutate(true_beta = beta_true[taxon_id],
           pp_sign = pmax(p_gt0, p_lt0),
           call = pp_sign >= thresh,
           correct_sign = sign(est) == sign(true_beta) | true_beta == 0,
           true_signal = true_beta != 0)
  tp <- sum(df$call & df$true_signal & df$correct_sign)
  fp <- sum(df$call & !df$true_signal)
  fn <- sum(!df$call & df$true_signal)
  list(
    table = df,
    metrics = tibble(
      TP = tp, FP = fp, FN = fn,
      precision = ifelse(tp + fp == 0, NA_real_, tp/(tp+fp)),
      recall = ifelse(tp + fn == 0, NA_real_, tp/(tp+fn))
    )
  )
}

# ---------------------------- Fit and tune -----------------------------------

# 1) Compile 
mod <- cmdstan_model("BALOR.stan", cpp_options = list(stan_threads = TRUE))

# 2) Make toy data (adjust N, M, sparsity as you like)
toy <- make_toy(
  N = 120, M = 150, K = 4,
  prop_signal = 0.10,   # 10 percent of taxa are truly DA
  pos_frac = 0.7,       # 70 percent of signals are positive
  effect_scale = 0.8,   # typical nonzero effect magnitude on log-odds
  seed = 42
)

# 3) Fit once with the current priors in the Stan file
fit <- mod$sample(
  data = toy$stan_data,
  chains = 4, parallel_chains = 4, threads_per_chain = 3,
  iter_warmup = 500, iter_sampling = 500,
  seed = 1, refresh = 100
)

# 4) Gather draws of beta, tau, nu
draws <- fit$draws(variables = c("beta", "tau", "nu")) |> as_draws_df()

# 5) Grid of alternative priors to evaluate 
tau_sd_grid <- c(0.5, 1.0, 2.0, 3.0)    # half-N(0, tau_sd)
nu_sc_grid  <- c(0.05, 0.10, 0.20)      # Laplace(0.5, scale = nu_sc)

results <- list()
for (s_tau in tau_sd_grid) {
  for (s_nu in nu_sc_grid) {
    w <- reweight_draws(draws, tau_sd_new = s_tau, nu_scale_new = s_nu)
    summ <- summarize_weighted(draws, w, toy$truth$beta_true, thresh = 0.95)
    results[[length(results) + 1]] <- summ$metrics |>
      mutate(tau_prior_sd = s_tau, nu_prior_scale = s_nu)
  }
}
tuning_table <- bind_rows(results) |>
  arrange(tau_prior_sd, nu_prior_scale)

print(tuning_table, n = nrow(tuning_table))


tuning_table |> mutate(F1 = 2 * precision * recall / (precision + recall)) |> arrange(desc(F1)) |> print(n = Inf)

# 6) Inspect how strong the reweighted inference is for top pick
best <- tuning_table |> mutate(F1 = 2 * precision * recall / (precision + recall)) |> slice_max(F1, n = 1)
best
w_best <- reweight_draws(draws, tau_sd_new = best$tau_prior_sd, nu_scale_new = best$nu_prior_scale)
best_summ <- summarize_weighted(draws, w_best, toy$truth$beta_true, thresh = 0.95)

# Show top 10 taxa by posterior sign probability
best_summ$table |>
  transmute(taxon_id, est, lwr95, upr95, pp_sign = pmax(p_gt0, p_lt0), true_beta) |>
  arrange(desc(pp_sign)) |>
  head(10) |>
  print(n = 10)
