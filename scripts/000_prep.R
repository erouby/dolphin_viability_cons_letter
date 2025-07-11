## ---------------------------
##
## Script name: 000_prep.R
##
## Purpose of script: Prepare the compiling of the models.
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
##  This script is build to prepare the compiling of all the models used. 
##  After the compiling, it will save the models in an appropriate folder. 
##  
## ---------------------------

start_time <- Sys.time()

### Compile models ----------------------------------------------------------------------------------------------

# Model 1 
M1 <- stan_model(file = paste(ModDir, "M1_null.stan", sep = "/"),
                             model_name = "M1_null"
)

writeLines("M1 compiled and stored")

# Model 2  
M2 <- stan_model(file = paste(ModDir, "M2_cov.stan", sep = "/"),
                          model_name = "M2_cov"
)

writeLines("M2 compiled and stored")

# Model 3 
M3 <- stan_model(file = paste(ModDir, "M3_random.stan", sep = "/"),
                            model_name = "M3_random"
)

writeLines("M3 compiled and stored")

# Model 4  
M4 <- stan_model(file = paste(ModDir, "M4_trend.stan", sep = "/"),
                                model_name = "M4_trend"
)

writeLines("M4 compiled and stored")

# Model 5 
M5 <- stan_model(file = paste(ModDir, "M5_random_trend.stan", sep = "/"),
                 model_name = "M5_random_trend"
)

writeLines("M5 compiled and stored")

# Model 6
M6 <- stan_model(file = paste(ModDir, "M6_cov_random.stan", sep = "/"),
                                model_name = "M6_cov_random"
)

writeLines("M6 compiled and stored")

# Model 7
M7 <- stan_model(file = paste(ModDir, "M7_cov_trend.stan", sep = "/"),
                                    model_name = "M7_cov_trend"
)

writeLines("M7 compiled and stored")

# Model 8 
M8 <- stan_model(file = paste(ModDir, "M8_cov_random_trend.stan", sep = "/"),
                                             model_name = "M8_cov_random_trend"
)

writeLines("M8 compiled and stored")

### Save the model fits in the dedicated repository
dir.create(paste(OutDir, "output_fits", sep = "/"))
save.image(paste(OutDir, "output_fits/models_compiled.RData", sep = "/"), safe = TRUE)

writeLines("All the models have compile")

source(paste(ScriptsDir, "001_analysis.R", sep = "/"))
