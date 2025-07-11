## ---------------------------
##
## Script name: 002_model_selection.R
##
## Purpose of script: Fit all the models.
##
## Author: Etienne Rouby, Floriane Plard, Matthieu Authier
##
## Date Created: 10-03-2023
##
## Copyright (c) Etienne Rouby, 2023
## Email: etiennerouby95@gmail.com
##
## ---------------------------
##
## Notes:  
##   
##  This script is build to fit all the models previously defined and compiled.  
##  This is the longest part of the job. For the moment, it is not possible to choose
##  if you want to fit only some models. All have to run. 
##  Then, a table with wAIC and LOOIC will be produce and save to select the better model. 
##
## ---------------------------

### Load models
load(paste(OutDir, "/output_fits","models_compiled.RData", sep = "/"))

# Parallelisation
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### Fits

### Create Directory for the fit to be saved
dir.create(paste(OutDir, "output_fits/models_fits", sep = "/"))

# Model 1 Both New SR
M1_fit = sampling(object = M1,
                  data = list(n_ind_surv = n_ind_surv,
                              #n_ind_repro = n_ind_repro,
                              #n_ind_diff = n_ind_diff,
                              #n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              #TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              #CENSORING = 2 -  SR_data,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1
                       ), 
                       pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "R_sq", "ppc", "log_lik_surv"),
                       chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                       control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )


### Save fit
save(list = c("M1_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M1_fit.RData", sep = "/")
     )

writeLines("Model 1 fitted. It took:")
print(get_elapsed_time(M1_fit))

# Model 2 
M2_fit = sampling(object = M2,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              n_cov = n_cov,
                              X = X,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1,
                              prior_scale_for_slope = 1
                       ), 
                       pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "slope", "R_sq", "log_lik_surv"),
                       chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                       control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M2_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M2_fit.RData", sep = "/")
     )

writeLines("Model 2 fitted. It took:")
print(get_elapsed_time(M2_fit))

# Model 3 

M3_fit = sampling(object = M3,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              YEAR = YEAR,
                              n_year = n,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1
                              ),
                  pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "sigma_year", "u_year", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M3_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M3_fit.RData", sep = "/")
     )

writeLines("Model 3 fitted. It took:")
print(get_elapsed_time(M3_fit))

# Model 4 

M4_fit = sampling(object = M4,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              YEAR = YEAR,
                              x_trend = x_trend,
                              n_year = n,
                              n_cov = n_cov,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1
                              ),
                  pars = c("intercept", "residual_sd", "beta", "trend", "sigma_frailty", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M4_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M4_fit.RData", sep = "/")
     )

writeLines("Model 4 fitted. It took:")
print(get_elapsed_time(M4_fit))

# Model 5 

M5_fit = sampling(object = M5,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              YEAR = YEAR,
                              x_trend = x_trend,
                              n_year = n,
                              n_cov = n_cov,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1
                  ),
                  pars = c("intercept", "residual_sd", "beta", "trend", "sigma_year", "u_year", "sigma_frailty", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
)

### Save the fit
save(list = c("M5_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M5_fit.RData", sep = "/")
)

writeLines("Model 5 fitted. It took:")
print(get_elapsed_time(M5_fit))

# Model 6 

M6_fit = sampling(object = M6,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              YEAR = YEAR,
                              n_year = n,
                              n_cov = n_cov,
                              X = X,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1,
                              prior_scale_for_slope = 1
                              ),
                  pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "slope", "sigma_year", "u_year", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M6_fit"),
     file = paste(OutDir, "output_fits/models_fits/M6_fit.RData", sep = "/")
     )

writeLines("Model 6 fitted. It took:")
print(get_elapsed_time(M6_fit))

# Model 7

M7_fit = sampling(object = M7,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 -  SR_data,
                              YEAR = YEAR,
                              x_trend = x_trend,
                              n_year = n,
                              n_cov = n_cov,
                              X = X,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1,
                              prior_scale_for_slope = 1
                              ),
                  pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "trend", "slope", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M7_fit"),
     file = paste(OutDir, "output_fits/models_fits/M7_fit.RData", sep = "/")
     )

writeLines("Model 7 fitted. It took:")
print(get_elapsed_time(M7_fit))

# Model 8

M8_fit = sampling(object = M8,
                  data = list(n_ind_surv = n_ind_surv,
                              # n_ind_repro = n_ind_repro,
                              # n_ind_diff = n_ind_diff,
                              # n_ind_total = n_ind_total,
                              LOGTIME = log(age_data),
                              # TIME =  SR_age_data,
                              SURVIVAL = rep(0, n_ind_surv),
                              # CENSORING = 2 - SR_data,
                              YEAR = YEAR,
                              n_year = n,
                              x_trend = x_trend,
                              n_cov = n_cov,
                              X = X,
                              prior_scale_for_intercept = 1.5,
                              prior_scale_for_sigma_frailty = log(2),
                              prior_scale_for_residual_sd = 0.1,
                              prior_scale_for_slope = 1
                              ),
                  pars = c("intercept", "residual_sd", "beta", "sigma_frailty", "slope", "trend", "sigma_year", "u_year", "R_sq", "log_lik_surv"),
                  chains = n_chains, iter = n_iter, warmup = n_warm, thin = n_thin,
                  control = list(adapt_delta = 0.9, max_treedepth = 15)
                  )

### Save the fit
save(list = c("M8_fit"), 
     file = paste(OutDir, "output_fits/models_fits/M8_fit.RData", sep = "/")
     )

writeLines("Model 8 fitted. It took:")
print(get_elapsed_time(M7_fit))

### WAIC SELECTION

# Creation de la table recapitulative des verifications de la qualite du modele
waic_df <- do.call('rbind',
                   lapply(list(M1_fit,
                               M2_fit,
                               M3_fit,
                               M4_fit,
                               M5_fit,
                               M6_fit,
                               M7_fit,
                               M8_fit), 
                          function(stanfit){ as.data.frame(do.call('cbind', 
                                                                   loo::waic(loo::extract_log_lik(stanfit = stanfit, 
                                                                                                  parameter_name = c("log_lik_surv")))[3:8])) }
                   )
)
waic_df$model <- c(1:8)

loo_df <- do.call('rbind',
                  lapply(list(M1_fit,
                              M2_fit,
                              M3_fit,
                              M4_fit,
                              M5_fit,
                              M6_fit,
                              M7_fit,
                              M8_fit), 
                         function(stanfit){ as.data.frame(do.call('cbind', 
                                                                  loo::loo(loo::extract_log_lik(stanfit = stanfit, 
                                                                                                parameter_name = c("log_lik_surv")))[5:10])) }
                  )
)
loo_df$model <- c(1:8)

### SAVE WAIC AND LOO

dir.create(paste(OutDir, "output_tables", sep = "/"))
dir.create(paste(OutDir, "output_tables/model_selection", sep = "/"))

tablewaic <- NULL
tableloo <- NULL

tablewaic <- rbind(tablewaic, waic_df)
tableloo <- rbind(tableloo, loo_df)

write.table(tablewaic, paste(OutDir, "/output_tables/model_selection/table_waic.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

write.table(tableloo, paste(OutDir, "/output_tables/model_selection/table_looic.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

writeLines(paste("All the models are fitted. The best model is N°:", filter(tablewaic, waic == min(tablewaic$waic))$model, sep = ""))

source(paste("scripts", "003_loop.R", sep = "/"))
