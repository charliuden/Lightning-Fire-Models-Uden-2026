//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
data {
  int<lower=0> N;           // Number of observations
  vector[N] S;              // Lightning strike rates (observations)
  vector[N] CAPE;            // Climate variable 1
}

parameters {
  real a_alpha;             // Intercept for alpha
  real b_alpha;             // Coefficient for CAPE in alpha
  real a_beta;              // Intercept for beta
  real b_beta;              // Coefficient for CAPE in beta
}

model {
  vector[N] alpha;          // Shape parameter for the gamma distribution
  vector[N] beta;           // Rate parameter for the gamma distribution

  // Define alpha and beta as transformed functions of SWR, T, RH, W, and P
  alpha = exp(a_alpha + b_alpha * CAPE);
  beta = exp(a_beta + b_beta * CAPE);

  // Likelihood
  S ~ gamma(alpha, beta);

  // Priors
  a_alpha ~ normal(0, 1);
  b_alpha ~ normal(0, 1);
  a_beta ~ normal(0, 1);
  b_beta ~ normal(0, 1);
}

generated quantities {
  vector[N] log_lik;        // Log-likelihood for each observation

  for (n in 1:N) {
    // Recalculate alpha and beta in the generated quantities block
    real alpha_n = exp(a_alpha + b_alpha * CAPE[n]);
    real beta_n = exp(a_beta + b_beta * CAPE[n]);
    log_lik[n] = gamma_lpdf(S[n] | alpha_n, beta_n);
  }
}
