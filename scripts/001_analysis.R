## ---------------------------
##
## Script name: 001_analysis.R
##
## Purpose of script: Load and reshape the data.
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
##  This script is build to load and reshape the data.  
##  Reshape process is mandatory, since it allows Rstan and Stan to run.
##  
## ---------------------------

writeLines("Load the data ...")

load("data/20250210_commondolphin.RData")

### Individuals with age and gonads
# dd_deldel_SR <- read_csv("data/dd_deldel_SR.csv", show_col_types = FALSE)

data_deldel_SR <- SR_cor2 %>% 
  select(num_collec, age_age, SR, SR_num, code_sexe, code_capt, annee) %>% 
  mutate(num_collec = as.character(num_collec),
         code_capt = as.character(code_capt),
         code_sexe = as.character(code_sexe),
         annee = as.numeric(annee)) %>% 
  na.omit()

data_deldel_SR <- data_deldel_SR[order(data_deldel_SR$age_age),]

### Individuals with age only 
dd_deldel <- anti_join(deldel, data_deldel_SR, by = "num_collec")

dd_deldel_random <- dd_deldel %>% 
  filter(!is.na(age_age) & !is.na(annee)) %>%  
  mutate(age_age = as.numeric(age_age))

data_deldel_surv <- dd_deldel_random %>% 
  select(num_collec, age_age, code_sexe, code_capt, annee) %>% 
  mutate(num_collec = as.character(num_collec),
         code_capt = as.character(code_capt),
         code_sexe = as.character(code_sexe),
         annee = as.numeric(annee)) %>% 
  na.omit()

data_deldel_surv <- data_deldel_surv[order(data_deldel_surv$age_age),]

rm(dd_deldel_random, dd_deldel)

### Join table
data_deldel <- left_join(data_deldel_surv, data_deldel_SR, by = c("num_collec", "age_age", "code_capt", "code_sexe", "annee"))
data_deldel <- data_deldel[order(data_deldel$SR_num),]
data_deldel_diff <- anti_join(data_deldel_SR, data_deldel_surv, by = c("num_collec"))
data_deldel_final <- rbind(data_deldel_diff, data_deldel)

### select now
for (i in (1:nrow(data_deldel_final))){
  if (data_deldel_final$age_age[i] == 0 ) {
    data_deldel_final$age_age[i] <- round(runif(1, 0.01, 0.2), 2)
  }
  
  if (data_deldel_final$code_capt[i] == 30 | data_deldel_final$code_capt[i] == 40 ) {
    data_deldel_final$COV1[i] <- 1
  }
  
  if (data_deldel_final$code_capt[i] == 0 | data_deldel_final$code_capt[i] == 52 ) {
    data_deldel_final$COV1[i] <- 0
  }
  
  if (data_deldel_final$code_sexe[i] == 1 ) {
    data_deldel_final$COV2[i] <- 1
  }
  
  if (data_deldel_final$code_sexe[i] == 2 | data_deldel_final$code_sexe[i] == 3) {
    data_deldel_final$COV2[i] <- 0
  }
}

# Age and SR for estimation
age_data <- data_deldel_final$age_age # Because of the i + n_ind_diff for the likelihood
SR_age_data <- filter(data_deldel_final, !is.na(SR_num))$age_age
SR_data <- filter(data_deldel_final, !is.na(SR_num))$SR_num

# Number of data
dataset <- "mine"
n_ind_surv <- nrow(data_deldel_final)
# n_ind_repro <- length(SR_data)
# n_ind_diff <- nrow(data_deldel_diff)
# n_ind_total <- nrow(data_deldel_final)

# Age for functions
age_modelise_surv <- log(seq(0.001, 31, length = 100)) # log(seq(0.001, 25, length = 100))
age_modelise_fec <- seq(0.001, 31, length = 100) # seq(0.001, 25, length = 100)

# Years and trend 
YEAR <- (data_deldel_final$annee - min(data_deldel_final$annee)) + 1
n <- max(YEAR) ### BECAREFUL IT IS NOT THE GOOD WAY TO SAY max(YEAR)
n_year <- max(YEAR)
x_trend <- seq(1:n_year)

# Covariates
cov <- c("cov1","cov2")
X <- matrix(c(data_deldel_final$COV1, data_deldel_final$COV2), ncol = 2)
n_cov <- ncol(X)

writeLines("Data are loaded. Now fit the model and select the better one.")

source(paste("scripts", "002_model_selection.R", sep = "/"))
