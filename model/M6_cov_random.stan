data {
 int<lower = 1> n_ind_surv;                             // number of individuals for survival
 vector[n_ind_surv] LOGTIME;                           // time to event (log scale) for survivorship
 int<lower = 1> n_year;                                 // number of years
 int<lower = 1> n_cov;                                  // number of covariates
 matrix[n_ind_surv, n_cov] X;                          // covariate matrix for survivorship
 int<lower = 1, upper = n_year> YEAR[n_ind_surv];      // Year index for survivorship data
 int<lower = 0, upper = 1> SURVIVAL[n_ind_surv];       // Survivorship index (all are dead)
 real<lower = 0> prior_scale_for_intercept;        // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_residual_sd;           // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_slope;            // depends on log TIME scale and TIME scale
 real<lower = 0> prior_scale_for_sigma_frailty;         // depends on log TIME scale and TIME scale
}

parameters {
 real unscaled_intercept;                // Intercept as real for both survivorship and reproduction
 real unscaled_slope[n_cov];             // Slope for covariate as real for both survivorship and reproduction
 real<lower = 0> unscaled_sigma2;             // Scale parameter for the frailty                   
 vector<lower = 0>[3] tau;                    // Parameter for beta2 distribution
 simplex[3] prop;        // variance partitioning
 vector<lower = 0>[n_ind_surv] Z;             // Individual Frailty
 vector[n_year] u_year;                    // random effect parameter
}

transformed parameters {
 real intercept;
 real sigma_year;
 real residual_sd;
 real R_sq;
 real sigma_frailty;
 vector[n_ind_surv] linear_predictor_surv;
 vector[n_ind_surv] residuals_surv;
 real slope[n_cov];
 real beta;
 // scale for frailty
 sigma_frailty = prior_scale_for_sigma_frailty * sqrt(prop[1] * unscaled_sigma2 / tau[1]);
 beta = 1 / sigma_frailty;
 // residual std. dev.
 residual_sd = prior_scale_for_residual_sd * sqrt(prop[2] * unscaled_sigma2 / tau[2]);
 // scale for year effect
 sigma_year = prior_scale_for_residual_sd * sqrt(prop[3] * unscaled_sigma2 / tau[3]);
 // slope
 for (n in 1:n_cov) {
 slope[n] = unscaled_slope[n] * prior_scale_for_slope;
}
 // intercept
 intercept = unscaled_intercept * prior_scale_for_intercept;

 // linear predictors
 for(i in 1:n_ind_surv) {
  // part without the covariates
  linear_predictor_surv[i] = intercept - (Z[i]  * sigma_frailty) + sigma_year * u_year[YEAR[i]];
  for (n in 1:n_cov) {
   // add each covariate
   linear_predictor_surv[i] += X[i, n] * slope[n];
  }
  // now linear_predictor is complete: compute residual
  residuals_surv[i] = LOGTIME[i] - linear_predictor_surv[i];
 }
 
 R_sq = 1 - variance(residuals_surv) / variance(LOGTIME);
}

model {
 // weakly informative priors
 unscaled_intercept ~ normal(0.0, 1.0);
 unscaled_sigma2 ~ gamma(0.5, 1.0);
 for (n in 1:n_cov) {
  unscaled_slope[n] ~ normal(0.0, 2);
 }
 for(i in 1:n_year) {
  u_year[i] ~ normal(0.0, 2);
 }
 tau ~ gamma(1.0, 1.0);
 Z ~ exponential(1.0);
 
 for (i in 1:n_ind_surv) {
  if(SURVIVAL[i] == 0) { 
    target += normal_lpdf(LOGTIME[i] | linear_predictor_surv[i], residual_sd);
  }
  else { // censored
    target += normal_lccdf(LOGTIME[i] | linear_predictor_surv[i], residual_sd);
  }
 }
}

generated quantities {
 vector[n_ind_surv] y_rep_surv;
 vector[n_ind_surv] log_lik_surv;

 for(i in 1:n_ind_surv) {
  y_rep_surv[i] = normal_rng(linear_predictor_surv[i], residual_sd);
  if(SURVIVAL[i] == 0) { 
   log_lik_surv[i] = normal_lpdf(LOGTIME[i] | linear_predictor_surv[i], residual_sd);
  }
  else { // censored
   log_lik_surv[i] = normal_lccdf(LOGTIME[i] | linear_predictor_surv[i], residual_sd);
  }
 } 
}
