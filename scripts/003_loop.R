## ---------------------------
##
## Script name: 003_loop.R
##
## Purpose of script: Obtain survivorship, fecundity and growth rate
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
##  This script is build to obtain the vital rates and the growth rate. 
##  First, it will use the better model fit. Depending on n_iter and n_warm, 
##  it is possible that better models varies between model 4, 6 and 7. 
##  In case the model 7 is the better, the script checks if it is significantly better than model 6 and 4. 
##  If not, the model with the lowest wAIC after model 7 is selected. 
##  We did that to be more parcimonious.  
##  Then, all the vital rates are computed and stored in the appropriate folder.  
##  
## ---------------------------

### Load the best model fits
table_waic <- read_delim("output/output_tables/model_selection/table_waic.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
  mutate(waic_low = waic - se_waic, 
         waic_upp = waic + se_waic)

loadRData <- function(fileName){ # This function must go in the function folder. But for now, keept it here. 
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

writeLines(paste("The best model seems to be the model", filter(table_waic, waic == min(table_waic$waic))$model, sep = " "))
writeLines("In case it is the model 7, the script evaluates if the wAIC difference is significative. If not, the second model wit hthe lowest wAIC is selected")

if (filter(table_waic, model == 7)$waic < filter(table_waic, model == 4)$waic & 
    filter(table_waic, model == 7)$waic < filter(table_waic, model == 6)$waic &
    filter(table_waic, model == 7)$waic > filter(table_waic, model == 4)$waic_low & filter(table_waic, model == 7)$waic < filter(table_waic, model == 4)$waic_upp &
    filter(table_waic, model == 7)$waic > filter(table_waic, model == 6)$waic_low & filter(table_waic, model == 7)$waic < filter(table_waic, model == 6)$waic_upp) {

  table_waic <- table_waic %>%
    filter(model != 7)
  
  model_fit <- loadRData(paste(OutDir, "output_fits/models_fits",
                               paste("M", filter(table_waic, waic == min(table_waic$waic))$model, "_fit.RData", sep = ""), sep = "/"))

} else {

  model_fit <- loadRData(paste(OutDir, "output_fits/models_fits",
                               paste("M", filter(table_waic, waic == min(table_waic$waic))$model, "_fit.RData", sep = ""), sep = "/"))
}


if (model_fit@model_name == "M4_fit") {
  # dir.create(paste(OutDir, "output_tables/vital_rates/M6_fit", sep = "/"))
  source(paste("source", "functions_M4.R", sep = "/"))
}

if (model_fit@model_name == "M6_fit") {
  # dir.create(paste(OutDir, "output_tables/vital_rates/M6_fit", sep = "/"))
  source(paste("source", "functions_M6.R", sep = "/"))
}

if (model_fit@model_name == "M7_cov_trend") {
  # dir.create(paste(OutDir, "output_tables/vital_rates/M7_fit", sep = "/"))
  source(paste("source", "functions_M7.R", sep = "/"))
}

# VITAL RATES -------------------------------------------------------------------------------------------------------------------------

writeLines("Now, its time to estimate the vital rates...")

if (model_fit@model_name == "M4_fit") { 
  
  survival_rates <- NULL
  
  for (what in c("survival", "hazard")) {
      table <- get_surv(stanfit = model_fit, 
                        logt = age_modelise_surv,
                        what = what,
                        n_year = 23
      )
      
      table$vital_rate <- what
      table$model <- model_fit@model_name
      
      survival_rates <- rbind(survival_rates, table)
  }
  
  writeLines("Survival rates have been estimated and saved.")
  
  # ### Fecundity
  # 
  # prop_mature <- NULL
  # 
  #   table <- get_repro(stanfit = model_fit, 
  #                      t = age_modelise_fec, 
  #                      alpha = 0.20, 
  #                      n_year = 23, 
  #                      type = "survival"
  #   )
  #   
  #   table$vital_rate <- "proportion_of_mature"
  #   table$model <- model_fit@model_name
  #   
  #   prop_mature <- rbind(prop_mature, table)
  # 
  # writeLines("Proportion of mature has been estimated and saved.")
  
}

if (model_fit@model_name == "M6_fit" | model_fit@model_name == "M7_cov_trend") {
  
  ### Survival
  
  survival_rates <- NULL
  
  for (what in c("survival", "hazard")) {
    for (cov in c("cov1_0", "cov1_1", "cov2_0", "cov2_1")) {
    table <- get_surv(stanfit = model_fit, 
              logt = age_modelise_surv,
              what = what,
              n_year = 23, 
              cov = cov
              )
    
    table$vital_rate <- what
    table$model <- model_fit@model_name
    
    survival_rates <- rbind(survival_rates, table)
    }
  }
  
  writeLines("Survival rates have been estimated and saved.")
  
  # ### Fecundity
  # 
  # prop_mature <- NULL
  # 
  # for (cov in c("cov1_0", "cov1_1", "cov2_0", "cov2_1")) {
  #     table <- get_repro(stanfit = model_fit, 
  #                      t = age_modelise_fec, 
  #                      alpha = 0.20, 
  #                      n_year = 23, 
  #                      cov = cov, 
  #                      type = "survival"
  #                      )
  #     
  #     table$vital_rate <- "proportion_of_mature"
  #     table$model <- model_fit@model_name
  #     
  #     prop_mature <- rbind(prop_mature, table)
  # }
  # 
  # writeLines("Proportion of mature has been estimated and saved.")
  
}

write.table(survival_rates, paste(OutDir, "/output_tables/vital_rates/", "survival_rates.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

# write.table(prop_mature, paste(OutDir, "/output_tables/vital_rates/", "prop_mature.txt", sep = ""), quote = TRUE, sep = "\t",
#             row.names = FALSE, col.names = TRUE
# )


# GROWTH RATE ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

writeLines("Now lets obtain the lambdas (growth rates)...")

library(demogR)
time_serie <- 100

if (model_fit@model_name == "M6_fit" | model_fit@model_name == "M7_cov_trend") { 

female_surv <- survivorship_profile <- read_delim("output/output_tables/vital_rates/survival_rates.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(vital_rate == "survival", cov == "cov2_0", year != "none") %>%
  mutate(year = as.numeric(year)) # Select only the females

female_fec <- read_csv("output/maturity_proportions.csv") # Select only the females

}

proj_years <- NULL
elas_years <- NULL
contrib_years <- NULL

### Obtain growth_rates using the dedicated function 
source(paste(FuncDir, "get_growth_rates.R", sep = "/"))

table_lambda <- get_growth_rates(data_surv = female_surv, data_fec = female_fec, n_year = 23)

write.table(table_lambda, paste(OutDir, "/output_tables/vital_rates/", "table_lambda.txt", sep = ""), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
)

writeLines("All lambdas have been estimated and saved.")

### Final words 
end_time <- Sys.time()
duration_time <- end_time - start_time

writeLines("Congrats !!! Analyses have been well performed.")
paste("All the process took:", round(duration_time/60,2), "minutes", sep = " ")



