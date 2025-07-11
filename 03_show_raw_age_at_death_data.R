## ---------------------------
##
## Script name: 03_show_data.R
##
## Purpose of script: Show data used to perform the analyzes
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
##  This script is build to obtain the distribution of ages and reproductive status at deaths. 
##  It is just graphical, to visualize the data used. 
##  
## ---------------------------

### Prepare the environment

lapply(c("sp", "tidyverse", "sf", "readxl", "ggplot2"), library, character.only = TRUE)

rm(list = ls())

WorkDir <- getwd() # Working directory
DataDir <- paste(WorkDir, "data", sep = "/") # Data folder
FuncDir <- paste(WorkDir, "source", sep = "/") # Functions and Adjacency matrix folder
ScriptsDir <- paste(WorkDir, "scripts", sep = "/") # Functions and Adjacency matrix folder
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder

### Load data ------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

### Make graph -----------------------------------------------------------------

### Age at death

dataset <- data_deldel_final %>%
  mutate(age = age_age,
  year = annee) %>%
  select(age, year) %>%
  na.omit() %>% 
  arrange(-desc(age))

graph_relative_dataset <- dataset %>%
  mutate(age = round(age)) %>%
  group_by(year, age) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) * 100,
         year = factor(year),
         age = factor(age, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                      "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25")),
         total = sum(n)) 

totals <- graph_relative_dataset %>%
  group_by(year) %>%
  summarize(total = sum(n))

# Age at death total 
# Round down the ages as requested
dataset_rounded <- dataset %>%
  mutate(age_rounded = floor(age))

# Count the frequency of each rounded age
age_counts <- dataset_rounded %>%
  group_by(age_rounded) %>%
  summarise(count = n()) %>%
  # Calculate the proportion
  mutate(proportion = count / sum(count) * 100)  # Convert to percentage

# Create the histogram with counts on top of bars
ggplot(age_counts, aes(x = age_rounded, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste("n =", count)), vjust = -0.5, size = 4) +
  labs(
    x = "Age (years)",
    y = "Proportion (%)"
  ) +
  theme_bw(base_size = 18) +
  scale_x_continuous(breaks = seq(0, max(age_counts$age_rounded), by = 1)) +
  # Add some extra space at the top for the labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Age at death by year in proportion 
ggplot(data = graph_relative_dataset, aes(x = year, y = freq, group = age, fill = age)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "mako", name = "Age at death") +
  theme_bw(base_size = 22) +
  xlab("Years") +
  ylab("Frequency (%)") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Make the legend span the full width of the plot
    legend.box.margin = margin(t = 20),
    legend.key.width = unit(1.5, "cm"),
    # Specify fewer rows for the legend
    legend.nrow = 2
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
