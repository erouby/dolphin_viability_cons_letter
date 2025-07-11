## ---------------------------
##
## Script name: 01_convergence.R
##
## Purpose of script: Obtain convergence parameters and compare wAICs
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
##  This script is build to obtain the convergence parameters of the model used to produce vital rates. 
##  It will make a traceplot of the appropriate parameters, estimate rhat and make a wAICs graph. 
##  Be careful to correctly specify the n_warm before plotting traceplot. 
##  
## ---------------------------


### Prepare the environment ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

lapply(c("sp", "lubridate", "tidyverse", "ggthemes", "rstan", "loo", "viridisLite", "bayesplot"), library, character.only = TRUE)

rm(list = ls())

WorkDir <- getwd() # Working directory
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder

loadRData <- function(fileName){ 
  # loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### Load the best model fit ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_waic <- read_delim("output/output_tables/model_selection/table_waic.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(waic_low = waic - se_waic, 
         waic_upp = waic + se_waic)

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

if (model_fit@model_name == "M4_trend") {
  dir.create(paste(OutDir, "output_graphs/convergence_graphs", sep = "/"))
  pars <- c("intercept", "residual_sd", "sigma_frailty", "trend")
}

if (model_fit@model_name == "M6_cov_random") {
  dir.create(paste(OutDir, "output_graphs/convergence_graphs", sep = "/"))
  pars <- c("intercept", "residual_sd", "sigma_frailty", "trend")
}

if (model_fit@model_name == "M7_cov_trend") {
  dir.create(paste(OutDir, "output_graphs/convergence_graphs", sep = "/"))
  pars <- c("intercept", "residual_sd", "sigma_frailty", "trend", "slope[1]")
}

## Define the same n_warmup as in the 00_run.R script
n_warmup <- 250
  
### Get the chains -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

posterior <- extract(model_fit, inc_warmup = TRUE, permuted = FALSE)

mcmc_trace(posterior, pars = pars, n_warmup = n_warmup,
           facet_args = list(nrow = 2, labeller = label_parsed)) +
  theme_bw(base_size = 24)

ggsave(paste(OutDir, "output_graphs/convergence_graphs", paste(unique(model_fit@model_name), "traceplot.pdf", sep = "_"), sep = "/"), width = 16, height = 8, dpi = 1200)
ggsave(paste(OutDir, "output_graphs/convergence_graphs", paste(unique(model_fit@model_name), "traceplot.png", sep = "_"), sep = "/"), width = 16, height = 8, dpi = 1200)

### Get the Rhats ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rhats <- rhat(model_fit, pars = pars)

color_scheme_set("brightblue") # see help("color_scheme_set")

mcmc_rhat(rhats, size = 5) + 
  theme_bw(base_size = 24) +
  yaxis_text(hjust = 1) +
  theme(legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(fill = "grey25",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold", size = 22),
        legend.key.size = unit(1, 'cm'),
        legend.key = element_rect(fill = "white", colour = "black")) +
  theme(legend.text = element_text(colour="white", size = 20))

ggsave(paste(OutDir, "output_graphs/convergence_graphs", paste(unique(model_fit@model_name), "rhat.pdf", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 1200)
ggsave(paste(OutDir, "output_graphs/convergence_graphs", paste(unique(model_fit@model_name), "rhat.png", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 1200)

### Get the WAIC graph -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_waic <- read_delim("output/output_tables/model_selection/table_waic.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(waic_low = waic - se_waic, 
         waic_upp = waic + se_waic)

ggplot() +
  geom_point(data = table_waic, aes(x = model, y = waic), size = 7, alpha = 1) +
  geom_errorbar(data = table_waic, aes(x = model, y = waic, ymin = waic_low, ymax = waic_upp)) +
  theme_bw(base_size = 24) +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey75",
                                         size = 0.5, linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold")) +
  theme(legend.text = element_text(colour="white")) +
  scale_x_continuous(breaks = seq(1, 7, by = 1), labels = seq(1, 7, by = 1)) +
  #scale_y_continuous(limits = c(0.955, 1.03), breaks = seq(0.955, 1.03, by = 0.01)) +
  ylab("wAIC value") +
  xlab("Model N°")

ggsave(paste(OutDir, "output_graphs/convergence_graphs", "waic.pdf", sep = "/"), width = 12, height = 8, dpi = 1200)
ggsave(paste(OutDir, "output_graphs/convergence_graphs", "waic.png", sep = "/"), width = 12, height = 8, dpi = 1200)

### Get the LOOIC graph -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_looic <- read_delim("output/output_tables/model_selection/table_looic.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(looic_low = looic - se_looic, 
         looic_upp = looic + se_looic)

ggplot() +
  geom_point(data = table_looic, aes(x = model, y = looic), size = 7, alpha = 1) +
  geom_errorbar(data = table_looic, aes(x = model, y = looic, ymin = looic_low, ymax = looic_upp)) +
  theme_bw(base_size = 24) +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey75",
                                         size = 0.5, linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold")) +
  theme(legend.text = element_text(colour="white")) +
  scale_x_continuous(breaks = seq(1, 7, by = 1), labels = seq(1, 7, by = 1)) +
  #scale_y_continuous(limits = c(0.955, 1.03), breaks = seq(0.955, 1.03, by = 0.01)) +
  ylab("looic value") +
  xlab("Model N°")

ggsave(paste(OutDir, "output_graphs/convergence_graphs", "looic.pdf", sep = "/"), width = 12, height = 8, dpi = 1200)
ggsave(paste(OutDir, "output_graphs/convergence_graphs", "looic.png", sep = "/"), width = 12, height = 8, dpi = 1200)
