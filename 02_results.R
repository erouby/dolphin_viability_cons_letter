## ---------------------------
##
## Script name: 02_results.R
##
## Purpose of script: Obtain vital rates graphs
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
##  This script is build to produce the graphs of the vital rates. 
##  It will automatically load the vital rates table. Then produce and save graphs. 
##  
## ---------------------------

### Prepare the environment -----------------------------------------------------------------------------------------------------------------------------------

lapply(c("tidyverse", "ggthemes", "viridisLite"), library, character.only = TRUE)

rm(list = ls())

WorkDir <- getwd() # Working directory
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder

### Get the graphs ---------------------------------------------------------------------------------------------------------------------

# Survival
## With Year Effect

survivorship_profile <- read_delim("output/output_tables/vital_rates/survival_rates.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

if (unique(survivorship_profile$model) == "M4_trend") {
  
  dir.create(paste(OutDir, "output_graphs/results_graphs", sep = "/"))
  
  survivorship_profile <- read_delim("output/output_tables/vital_rates/survival_rates.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    filter(vital_rate == "survival") %>%
    mutate(year = recode(year, 
                         "none" = "none", "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
                         "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
                         "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019"))
  
  }

if (unique(survivorship_profile$model) == "M6_cov_trend" | unique(survivorship_profile$model) == "M7_cov_trend") {
  
  dir.create(paste(OutDir, "output_graphs/results_graphs", sep = "/"))
  
  survivorship_profile <- read_delim("output/output_tables/vital_rates/survival_rates.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    filter(vital_rate == "survival", cov == "cov2_0") %>%
    mutate(year = recode(year, 
                         "none" = "none", "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
                         "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
                         "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019"))
  
}

ggplot() +
 # geom_ribbon(aes(x = age, ymin = lower, ymax = upper, group = year, fill = year), alpha = 0.2) +
  geom_line(data = filter(survivorship_profile, method == "Stan"), aes(x = age, y = mean, group = desc(year), color = year), size = 1, alpha = 0.6) +
  # geom_line(data = filter(survivorship_profile, method == "Stan" & year == "none"), aes(x = age, y = mean), size = 1, linetype = "dashed", color = "black") +
  scale_color_viridis_d(name = "Period effect") +
  theme_bw(base_size = 24) +
  theme(legend.position = c(.98, .98),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(fill = "grey25",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold", size = 16),
        legend.key.size = unit(0.8, 'cm'),
        legend.key = element_rect(fill = "white", colour = "black")) +
  theme(legend.text = element_text(colour="white", size = 14)) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  ylab("Survivorship (lx)") +
  xlab("Age (years)")

ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "survivorship.pdf", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)
ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "survivorship.png", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)

# Fecundity
## With Year Effect

# fecundity_profile <- read_delim("output/output_tables/vital_rates/prop_mature.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

if (unique(survivorship_profile$model) == "M4_trend") {
  
  fecundity_profile <- read_delim("output/output_tables/vital_rates/prop_mature.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # %>%
    # mutate(year = recode(year, 
                         # "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
                         # "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
                         # "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019"))

  }

# if (unique(survivorship_profile$model) == "M6_cov_trend" | unique(survivorship_profile$model) == "M7_cov_random_trend") {

  # fecundity_profile <- read_delim("output/output_tables/vital_rates/prop_mature.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  # filter(cov == "cov2_0") %>%
  # mutate(year = recode(year, 
  #                      "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
  #                      "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
  #                      "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019"))
  # 
#}

# ggplot() +
#   geom_ribbon(data = filter(fecundity_profile), aes(x = age, ymin = lower, ymax = upper), alpha = 0.2) +
#   # geom_line(data = filter(fecundity_profile, year != "none"), aes(x = age, y = mean, group = year, color = year), size = 1, alpha = 0.6) +
#   geom_line(data = filter(fecundity_profile), aes(x = age, y = mean), size = 1, color = "black") +
#   #ggtitle("Proportion of mature females: year effect") +
#   scale_color_viridis_d(name = "Period effect") +
#   theme_bw(base_size = 24) +
#   theme(legend.position = c(.99, .64),
#         legend.justification = c("right", "top"),
#         legend.box.just = "right",
#         legend.margin = margin(8, 8, 8, 8),
#         legend.background = element_rect(fill = "grey25",
#                                          size=0.5, linetype="solid", 
#                                          colour ="black"),
#         legend.title = element_text(colour="white", face="bold", size = 16),
#         legend.key.size = unit(0.8, 'cm'),
#         legend.key = element_rect(fill = "white", colour = "black")) +
#   theme(legend.text = element_text(size = 14, colour="white")) +
#   scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 2)) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), labels = seq(0, 100, by = 10)) +
#   ylab("Proportion of mature females (Pfx)") +
#   xlab("Age (years)")
# 
# ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(fecundity_profile$model), "asm.pdf", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)
# ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(fecundity_profile$model), "asm.png", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)

## Growth rate

table_lambda <- read_delim("output/output_tables/vital_rates/table_lambda.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(year = recode(Cohort, 
                       "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
                       "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
                       "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019"))

ggplot() +
  # geom_ribbon(aes(x = age, ymin = lower, ymax = upper, group = year, fill = year), alpha = 0.2) +
  geom_point(data = table_lambda, aes(x = year, y = Lambda, group = year), size = 4, alpha = 1) +
  geom_errorbar(data = table_lambda, aes(x = year, y = Lambda, ymin= Lambda_low, ymax = Lambda_upp, group = year), width = 0) +
  geom_hline(yintercept = 1, size = 1, color = "red", linetype = "dotted") +
  # geom_smooth(method = lm) +
  # geom_smooth(method = lm, se = FALSE, col = "navy", size = 2) +
  # geom_hline(yintercept = 0.8684294, linetype = "dashed", size = 2) +
  # geom_hline(data = filter(projection_table, Cohort == "mean"), aes(x = Cohort, y = Lambda), size = 1.2, color = "black") +
  scale_color_viridis_d(name = "Period effect") +
  theme_bw(base_size = 24) +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey25",
                                         size = 0.5, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold", size = 16),
        legend.key.size = unit(0.8, 'cm'),
        legend.key = element_rect(fill = "white", colour = "black")) +
  theme(legend.text = element_text(colour="white", size = 14)) +
  scale_x_discrete(breaks = seq(1997, 2019, by = 2), labels = seq(1997, 2019, by = 2)) +
  scale_y_continuous(limits = c(min(table_lambda$Lambda_low), max(table_lambda$Lambda_upp)), 
                     breaks = seq(round(min(table_lambda$Lambda_low), 3), round(max(table_lambda$Lambda_upp), 3), by = 0.01)) +
  ylab(expression(paste("Population growth rate (", lambda, ")"))) +
  xlab("Year (period)")

ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "growth_rate.pdf", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)
ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "growth_rate.png", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)

## Relative Growth rate

table_lambda <- read_delim("output/output_tables/vital_rates/table_lambda.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(year = recode(Cohort, 
                       "1" = "1997", "2" = "1998", "3" = "1999", "4" = "2000", "5" = "2001", "6" = "2002", "7" = "2003",
                       "8" = "2004", "9" = "2005", "10" = "2006", "11" = "2007", "12" = "2008", "13" = "2009", "14" = "2010", "15" = "2011",
                       "16" = "2012", "17" = "2013", "18" = "2014", "19" = "2015", "20" = "2016", "21" = "2017", "22" = "2018", "23" = "2019")) %>%
  # Calculate relative change from 1997 baseline
  mutate(baseline_lambda = first(Lambda),
         relative_change = ((Lambda - baseline_lambda) / baseline_lambda) * 100,
         relative_change_low = ((Lambda_low - baseline_lambda) / baseline_lambda) * 100,
         relative_change_upp = ((Lambda_upp - baseline_lambda) / baseline_lambda) * 100)

ggplot() +
  geom_point(data = table_lambda, aes(x = year, y = relative_change, group = year), size = 4, alpha = 1) +
  geom_errorbar(data = table_lambda, aes(x = year, y = relative_change, ymin = relative_change_low, ymax = relative_change_upp, group = year), width = 0) +
  geom_hline(yintercept = 0, size = 1, color = "red", linetype = "dotted") +
  scale_color_viridis_d(name = "Period effect") +
  theme_bw(base_size = 24) +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey25",
                                         size = 0.5, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.title = element_text(colour="white", face="bold", size = 16),
        legend.key.size = unit(0.8, 'cm'),
        legend.key = element_rect(fill = "white", colour = "black")) +
  theme(legend.text = element_text(colour="white", size = 14)) +
  scale_x_discrete(breaks = seq(1997, 2019, by = 2), labels = seq(1997, 2019, by = 2)) +
  scale_y_continuous(limits = c(min(table_lambda$relative_change_low), max(table_lambda$relative_change_upp))) +
  ylab("Relative change in population growth rate (%)") +
  xlab("Year (period)")

ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "growth_rate_relative.pdf", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)
ggsave(paste(OutDir, "output_graphs/results_graphs", paste(unique(survivorship_profile$model), "growth_rate_relative.png", sep = "_"), sep = "/"), width = 12, height = 8, dpi = 600)
