## ---------------------------
##
## Script name: 00_run.R
##
## Purpose of script: Run the demographical analyzes for the common dolphin in the North-East Atlantic ocean.
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
##  This script is the root to run all the analyzes. Just press run and wait. 
##  Be careful of what n_iter and n_warm you choose. Depending on these, it may take some time. 
##  Approximately 1-2 hours with n_iter = 5000 and n_warm = 2500. 
##  If you just want to try the script and see results, use n_iter = 2000 and n_warm = 1000
##
##  WARNING : Before running, you must have Rstan and Stan installed on your computer. 
## ---------------------------

### Prepare the environment

lapply(c("sp", "tidyverse", "sf", "rstan", "loo", "readxl"), library, character.only = TRUE)

rm(list = ls())

WorkDir <- getwd() # Working directory
DataDir <- paste(WorkDir, "data", sep = "/") # Data folder
FuncDir <- paste(WorkDir, "source", sep = "/") # Functions and Adjacency matrix folder
ScriptsDir <- paste(WorkDir, "scripts", sep = "/") # Functions and Adjacency matrix folder
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder

### Stan parameters and model fit
## MCMC settings
n_chains <- 4
n_thin <- 10
# n_iter <- 10000 * n_thin
# n_warm <- 5000 * n_thin

n_iter <- 5000 #* n_thin
n_warm <- 2500 #* n_thin

### Execution safety command

if (dir.exists(paste(OutDir, "output_fits", "models_fits", sep = "/"))) {
  switch(menu(c("Yes", "No"), title = "A path with models fits already exists. Do you want to overwrite this folder and start a new run?"),
         source(paste(ScriptsDir, "000_prep.R", sep = "/")))
} else {
  source(paste(ScriptsDir, "000_prep.R", sep = "/"))
}
