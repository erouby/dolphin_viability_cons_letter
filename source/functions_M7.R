
### Survival function
get_surv <- function(stanfit, logt = NULL, alpha = 0.20, y_sim = NULL, n_year = NULL, cov = NULL, what = "survival") {
  
  surv_km <- function(y_sim) {
    with(survival::survfit(survival::Surv(Value, Status) ~ 1, 
                           data = data.frame(Value = y_sim, 
                                             Status = rep(1, length(y_sim))
                           )
    ),
    data.frame(lower = lower,
               mean = surv,
               upper = upper,
               age = time,
               method = "KM"
    )
    )
  }
  
  surv_km_year <- function(y_sim) {
    with(survival::survfit(survival::Surv(Value, Status) ~ 1, 
                           data = data.frame(Value = y_sim, 
                                             Status = rep(1, length(y_sim))
                           )
    ),
    data.frame(lower = lower,
               mean = surv,
               upper = upper,
               age = time,
               method = "KM",
               year = "all"
    )
    )
  }
  
  ### Useful function 
  surv_stan <- function(stanfit, logt = NULL, alpha = 0.20, y_sim = NULL, y = NULL, what = "survival") {
    ### useful functions
    get_ci <- function(x) {
      c(coda::HPDinterval(coda::as.mcmc(x), prob = 1 - alpha)[1], 
        mean(x, na.rm = TRUE), 
        coda::HPDinterval(coda::as.mcmc(x), prob = 1 - alpha)[2]
      )
    }
    
    # sanity checks
    if(is.null(logt)) {
      if(is.null(y_sim)) {
        stop("Must provide either a vector of logt or a dataset y_sim")
      }
      else{
        logt <- log(survival::survfit(survival::Surv(Value, Status) ~ 1, 
                                      data = data.frame(Value = y_sim,
                                                        Status = rep(1, length(y_sim))
                                      ))$time
        )
      }
    }
    
    if (covariable == "cov1_0"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend)) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        
        # annual survival rate fct
        ann_surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt[i+1], mu[i+1], sd[i+1])) - exp(beta[i+1] * (logt[i+1] - mu[i+1]) + 0.5 * (beta[i+1] * sd[i+1])^2) * (1 - pnorm(beta[i+1] * sd  + (logt[i+1]-mu[i+1])/sd[i+1])) / (1 - pnorm(logt[i], mu[i], sd[i])) - exp(beta[i] * (logt[i] - mu[i]) + 0.5 * (beta[i] * sd[i])^2) * (1 - pnorm(beta[i] * sd[i]  + (logt[i]-mu[i])/sd[i]))
        }
        
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard", "ann_surv")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "ann_surv") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 ann_surv(logt = logt,
                                          mu = mu[row_i],
                                          sd = sd[row_i],
                                          beta = beta[row_i],
                                          alpha = alpha[row_i],
                                          sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov2_0"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend)) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov1_1"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept + rstan::extract(stanfit, 'slope')$slope[,1]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend)) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }
    
    if (covariable == "cov2_1"){
      
      ### extract posteriors
      if(any(stanfit@model_pars == "beta")) {
        beta <- drop(rstan::extract(stanfit, 'beta')$beta)
        mu <- drop(rstan::extract(stanfit, 'intercept')$intercept + rstan::extract(stanfit, 'slope')$slope[,2]) + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend)) * (y - mean(x_trend)) / sd(x_trend))
        sd <- drop(rstan::extract(stanfit, 'residual_sd')$residual_sd)
        alpha <- rep(0, length(mu))
        sigma <- rep(0, length(mu))
        
        # survivorship fct
        surv <- function(logt, mu, sd, beta, sigma, alpha) {
          (1 - pnorm(logt, mu, sd)) - exp(beta * (logt - mu) + 0.5 * (beta * sd)^2) * (1 - pnorm(beta * sd  + (logt-mu)/sd))
        }
        # hazard fct
        hazard <- function(logt, mu, sd, beta, sigma, alpha) {
          # pdf of a Reed distribution
          reed_pdf <- function(logt, mu, sd, beta) {
            beta * exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * (1 - pnorm(beta * sd + (logt - mu) / sd)) / exp(logt)
          }
          # survival fct of a Reed distribution
          reed_phi <- function(logt, mu, sd, beta) {
            (1 - pnorm(logt, mu, sd)) - 
              exp(beta * (logt - mu) + 0.5 * (beta * sd) * (beta * sd)) * 
              (1 - pnorm(beta * sd + (logt - mu) / sd))
          }
          # hazard of a reed distribution
          reed_hz <- function(logt, mu, sd, beta) {
            reed_pdf(logt, mu, sd, beta) / reed_phi(logt, mu, sd, beta)
          }
          return(reed_hz(logt, mu, sd, beta))
        }
      }
      
      ### survival
      if(what %in% c("survival", "hazard")) {
        if(what == "survival") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 surv(logt = logt,
                                      mu = mu[row_i],
                                      sd = sd[row_i],
                                      beta = beta[row_i],
                                      alpha = alpha[row_i],
                                      sigma = sigma[row_i]
                                 )
                               }
          )
        }
        if(what == "hazard") {
          y_pred <- purrr::map(1:length(mu), 
                               function(row_i) {
                                 hazard(logt = logt,
                                        mu = mu[row_i],
                                        sd = sd[row_i],
                                        beta = beta[row_i],
                                        alpha = alpha[row_i],
                                        sigma = sigma[row_i]
                                 )
                               }
          )
        }
      }
      else {
        stop("Must choose either 'survival' or 'hazard' for arg. 'what'")
      }
      
      ### formating
      y_pred <- do.call('rbind', y_pred)
      output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
      names(output) <- c("lower", "mean", "upper")
      output$age <- exp(logt)
      output$method <- "Stan"
      output$year <- y
      return(output)
    }  
    
  }
  
    surv_table <- NULL
    hazard_table <- NULL
    
    if (what == "survival" & cov == "cov1_1") {
      
      for (y in 1:n_year) {
      
        ### Surv
        covariable <- "cov1_1"
        surv_df <- rbind(cbind(surv_km(y_sim = filter(data_deldel_final, COV1 == 1)$age_age), year = rep("none")),
                                surv_stan(stanfit = model_fit, 
                                          logt = age_modelise_surv,
                                          what = "survival",
                                          y = y
                                )
        )
        
        surv_df$n_ind <- n_ind_surv
        surv_df$cov <- covariable
        
        surv_table <- rbind(surv_table, surv_df)
        
      }
      
      return(surv_table)
      
    }
      
    if (what == "hazard" & cov == "cov1_1") { 
        for (y in 1:n_year){
        covariable <- "cov1_1"
        ### Hazard
        hazard_df <- surv_stan(stanfit = model_fit, 
                                      logt = age_modelise_surv,
                                      what = "hazard",
                                      y = y
        )
        hazard_df$n_ind <- n_ind_surv
        hazard_df$cov <- covariable
        
        hazard_table <- rbind(hazard_table, hazard_df)
        
        }
        
        return(hazard_table)
        
      }
    
    if (what == "survival" & cov == "cov1_0") {
        for (y in 1:n_year){
        ### Surv
        covariable <- "cov1_0"
        surv_df <- rbind(cbind(surv_km(y_sim = filter(data_deldel_final, COV1 == 0)$age_age), year = rep("none")),
                                surv_stan(stanfit = model_fit, 
                                          logt = age_modelise_surv,
                                          what = "survival",
                                          y = y
                                )
        )
        
        surv_df$n_ind <- n_ind_surv
        surv_df$cov <- covariable
        
        surv_table <- rbind(surv_table, surv_df)
        }
        
        return(surv_table)
      }
       
    if (what == "hazard" & cov == "cov1_0") {
        for (y in 1:n_year){
        covariable <- "cov1_0"
        ### Hazard
        hazard_df <- surv_stan(stanfit = model_fit, 
                                      logt = age_modelise_surv,
                                      what = "hazard",
                                      y = y
        )
        
        hazard_df$n_ind <- n_ind_surv
        hazard_df$cov <- covariable
        
        hazard_table <- rbind(hazard_table, hazard_df)
        
        }
        return(hazard_table)
    }
    
    if (what == "survival" & cov == "cov2_1") {
      
      for (y in 1:n_year){
        
        ### Surv
        covariable <- "cov2_1"
        surv_df <- rbind(cbind(surv_km(y_sim = filter(data_deldel_final, COV2 == 1)$age_age), year = rep("none")),
                         surv_stan(stanfit = model_fit, 
                                   logt = age_modelise_surv,
                                   what = "survival",
                                   y = y
                         )
        )
        
        surv_df$n_ind <- n_ind_surv
        surv_df$cov <- covariable
        
        surv_table <- rbind(surv_table, surv_df)
        
      }
      
      return(surv_table)
      
    }
    
    if (what == "hazard" & cov == "cov2_1") { 
      for (y in 1:n_year){
        covariable <- "cov2_1"
        ### Hazard
        hazard_df <- surv_stan(stanfit = model_fit, 
                               logt = age_modelise_surv,
                               what = "hazard",
                               y = y
        )
        hazard_df$n_ind <- n_ind_surv
        hazard_df$cov <- covariable
        
        hazard_table <- rbind(hazard_table, hazard_df)
        
      }
      
      return(hazard_table)
      
    }
    
    if (what == "survival" & cov == "cov2_0") {
      for (y in 1:n_year){
        ### Surv
        covariable <- "cov2_0"
        surv_df <- rbind(cbind(surv_km(y_sim = filter(data_deldel_final, COV2 == 0)$age_age), year = rep("none")),
                                surv_stan(stanfit = model_fit, 
                                          logt = age_modelise_surv,
                                          what = "survival",
                                          y = y
                                )
        )
        
        surv_df$n_ind <- n_ind_surv
        surv_df$cov <- covariable
        
        surv_table <- rbind(surv_table, surv_df)
      }
      
      return(surv_table)
    }
    
    if (what == "hazard" & cov == "cov2_0") {
      for (y in 1:n_year){
        covariable <- "cov2_0"
        ### Hazard
        hazard_df <- surv_stan(stanfit = model_fit, 
                                      logt = age_modelise_surv,
                                      what = "hazard",
                                      y = y
        )
        
        hazard_df$n_ind <- n_ind_surv
        hazard_df$cov <- covariable
        
        hazard_table <- rbind(hazard_table, hazard_df)
        
      }
      return(hazard_table)
    }
    
}

# ### Proportion of mature function
# get_repro <- function(stanfit, t = NULL, alpha = 0.20, n_year = NULL, cov = NULL, type = "survival") {
# 
#   fec_table <- NULL
#   
#   ### Get useful function
#   repro <- function(t, stanfit, type = "hazard", covariable = NULL, y = NULL) {
#   # for getting confidence interval
#   get_ci <- function(x) {
#     c(coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[1], 
#       mean(x), 
#       coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[2])
#   }
#   
#   ### hazard
#   hz_fec <- function(t, mu, alpha) {
#     (alpha/(exp(-mu/alpha))^alpha) * t^(alpha - 1)
#   }
#   
#   ### Obtention taux fecondite
#   taux_fec <- function(t, mu, TG, alpha) {
#     (1-(exp(-(t/exp(-mu/alpha))^alpha))) * TG
#   }
#   
#   ### survival
#   surv_fec <- function(t, mu, alpha) {
#     1 - (exp(-(t/exp(-mu/alpha))^alpha))
#   }
#   
#   ### probability density function
#   pdf_fec <- function(t, mu, alpha) {
#     (alpha * t^(alpha - 1)/(exp(-mu/alpha))^alpha) * exp (-(t/exp(-mu/alpha))^alpha)
#   }
#   
#   if (covariable == "cov1_0") {
#     
#     ### Asignements of values
#     mu <- rstan::extract(stanfit, 'intercept')$intercept[,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend))
#     alpha <- rstan::extract(stanfit, 'alpha')$alpha
#     
#     # compute hazard at t
#     if(type == "hazard") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              hz_fec(t = t,
#                                     mu = mu[row_i], 
#                                     alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute cumulative survival at t
#     if(type == "survival") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              surv_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute fecondity rate
#     if(type == "fecondite") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              taux_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i],
#                                       TG = rstan::extract(stanfit, 'TG')$TG[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     
#     # compute pdf
#     if(type == "pdf") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept), 
#                            function(row_i) {
#                              pdf_fec(t = t,
#                                      mu = rstan::extract(stanfit, 'intercept')$intercept[row_i,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend)), 
#                                      alpha = rstan::extract(stanfit, 'alpha')$alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     return(output)
#   }
#   
#   if (covariable == "cov2_0") {
#     ### Asignements of values
#     mu <- rstan::extract(stanfit, 'intercept')$intercept[,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend))
#     alpha <- rstan::extract(stanfit, 'alpha')$alpha
#     
#     # compute hazard at t
#     if(type == "hazard") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              hz_fec(t = t,
#                                     mu = mu[row_i], 
#                                     alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute cumulative survival at t
#     if(type == "survival") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              surv_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute fecondity rate
#     if(type == "fecondite") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              taux_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i],
#                                       TG = rstan::extract(stanfit, 'TG')$TG[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     
#     # compute pdf
#     if(type == "pdf") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              pdf_fec(t = t,
#                                      mu = rstan::extract(stanfit, 'intercept')$intercept[row_i,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend)), 
#                                      alpha = rstan::extract(stanfit, 'alpha')$alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     return(output)
#   }
#   
#   if (covariable == "cov1_1") {
#     ### Asignements of values
#     mu <- rstan::extract(stanfit, 'intercept')$intercept[,2] + rstan::extract(stanfit, 'slope')$slope[,1,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend))
#     alpha <- rstan::extract(stanfit, 'alpha')$alpha
#     
#     # compute hazard at t
#     if(type == "hazard") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              hz_fec(t = t,
#                                     mu = mu[row_i], 
#                                     alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute cumulative survival at t
#     if(type == "survival") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              surv_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute fecondity rate
#     if(type == "fecondite") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              taux_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i],
#                                       TG = rstan::extract(stanfit, 'TG')$TG[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     
#     # compute pdf
#     if(type == "pdf") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              pdf_fec(t = t,
#                                      mu = rstan::extract(stanfit, 'intercept')$intercept[row_i,2] + rstan::extract(stanfit, 'slope')$slope[,1,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend)), 
#                                      alpha = rstan::extract(stanfit, 'alpha')$alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     return(output)
#   }
#   
#   if (covariable == "cov2_1") {
#     ### Asignements of values
#     mu <- rstan::extract(stanfit, 'intercept')$intercept[,2] + rstan::extract(stanfit, 'slope')$slope[,2,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend))
#     alpha <- rstan::extract(stanfit, 'alpha')$alpha
#     
#     # compute hazard at t
#     if(type == "hazard") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              hz_fec(t = t,
#                                     mu = mu[row_i], 
#                                     alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute cumulative survival at t
#     if(type == "survival") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              surv_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     # compute fecondity rate
#     if(type == "fecondite") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              taux_fec(t = t,
#                                       mu = mu[row_i], 
#                                       alpha = alpha[row_i],
#                                       TG = rstan::extract(stanfit, 'TG')$TG[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     
#     
#     # compute pdf
#     if(type == "pdf") {
#       y_pred <- purrr::map(1:length(rstan::extract(stanfit, 'intercept')$intercept[,2]), 
#                            function(row_i) {
#                              pdf_fec(t = t,
#                                      mu = rstan::extract(stanfit, 'intercept')$intercept[row_i,2] + rstan::extract(stanfit, 'slope')$slope[,2,2] + 1 * ((drop(rstan::extract(stanfit, 'trend')$trend[,2])) * (y - mean(x_trend)) / sd(x_trend)), 
#                                      alpha = rstan::extract(stanfit, 'alpha')$alpha[row_i]
#                              )
#                            }
#       )
#       y_pred <- do.call('rbind', y_pred)
#       output <- as.data.frame(t(apply(y_pred, 2, get_ci)))
#       names(output) <- c("lower", "mean", "upper")
#       output$age <- t
#       output$year <- y
#     }
#     return(output)
#   }
#   
# }
#     
#   if (cov == "cov1_1") {
#     for (y in 1:n_year) {
#     ### fec
#     covariable <- "cov1_1"
#     fec_df <- repro(stanfit = model_fit, 
#                         t = age_modelise_fec,
#                         type = "survival",
#                         covariable = cov,
#                         y = y
#     )
#     
#     fec_df$n_ind <- n_ind_repro
#     fec_df$age <- age_modelise_fec
#     fec_df$cov <- covariable
#     
#     fec_table <- rbind(fec_table, fec_df)
#     
#     }
#     return(fec_table)
#   }
#   
#   if (cov == "cov1_0") {
#     for (y in 1:n_year) {
#     covariable <- "cov1_0"
#     fec_df <- repro(stanfit = model_fit, 
#                         t = age_modelise_fec,
#                         type = "survival",
#                         covariable = cov,
#                         y = y
#     )
#     
#     fec_df$n_ind <- n_ind_repro
#     fec_df$age <- age_modelise_fec
#     fec_df$cov <- covariable
#     
#     fec_table <- rbind(fec_table, fec_df)
#     
#     }
#     return(fec_table)
#   }
#   
#   if (cov == "cov2_1") {
#     for (y in 1:n_year) {
#     covariable <- "cov2_1"
#     fec_df <- repro(stanfit = model_fit, 
#                         t = age_modelise_fec,
#                         type = "survival",
#                         covariable = cov,
#                         y = y
#     )
#     
#     fec_df$n_ind <- n_ind_repro
#     fec_df$age <- age_modelise_fec
#     fec_df$cov <- covariable
#     
#     fec_table <- rbind(fec_table, fec_df)
#     
#     }
#     return(fec_table)
#   }
#   
#   if (cov == "cov2_0") {
#     for (y in 1:n_year) {
#     covariable <- "cov2_0"
#     fec_df <- repro(stanfit = model_fit, 
#                         t = age_modelise_fec,
#                         type = "survival",
#                         covariable = cov,
#                         y = y
#     )
#     
#     fec_df$n_ind <- n_ind_repro
#     fec_df$age <- age_modelise_fec
#     fec_df$cov <- covariable
#     
#     fec_table <- rbind(fec_table, fec_df)
#     
#     }
#     return(fec_table)
#   }
#     
# }
