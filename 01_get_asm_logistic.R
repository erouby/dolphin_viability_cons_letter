##--------------------------------------------------------------------------------------------------------
## SCRIPT : Common dolphin data
##
## Authors : Matthieu Authier & Etienne Rouby
## Last update : 2025-02-10
##
## R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
## Copyright (C) 2022 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

lapply(c("tidyverse", "readxl", "ggplot2"), library, character.only = TRUE)
rm(list = ls())

###
writeLines("Load the data ...")

load("data/20250210_commondolphin.RData")

# Filter only females -------------------------------
SR_females <- SR_cor2 %>% 
  filter(code_sexe == 2 | code_sexe == 3)

# Basic logistic regression -----------------------------
model <- glm(SR_num == 1 ~ age_age, data = SR_females, family = binomial)
summary(model)

# Calculate age at 50% probability of maturity (P50)
P50 <- -coef(model)[1]/coef(model)[2]
print(paste("Age at 50% maturity (P50):", round(P50, 2), "years"))

# Create a plot with confidence intervals
# Create predictions for a range of ages
ages <- seq(min(SR_females$age_age), max(SR_females$age_age), length.out = 100)
new_data <- data.frame(age_age = ages)

# Get predictions with standard errors for confidence intervals
pred <- predict(model, newdata = new_data, type = "link", se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se <- pred$se.fit

# Transform to response scale (probabilities)
new_data$prob <- plogis(new_data$fit)
new_data$lower <- plogis(new_data$fit - 1.96 * new_data$se)
new_data$upper <- plogis(new_data$fit + 1.96 * new_data$se)

# Plot with confidence interval ribbon, percentage y-axis, and P50 line
ggplot() +
  geom_ribbon(data = new_data, aes(x = age_age, ymin = lower*100, ymax = upper*100), 
              fill = "gray80", alpha = 0.5) +
  geom_line(data = new_data, aes(x = age_age, y = prob*100), color = "black", size = 1) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "red") + 
  geom_vline(xintercept = P50, linetype = "dotted", color = "red") +
  annotate("text", x = 15, y = 60, 
           label = paste0("ASM = ", round(P50, 1), " years"), 
           color = "red", size = 5) +
  labs(x = "Age (years)",
       y = "Proportion of mature females (%)") +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(0,22), breaks = seq(0,22,2)) +
  scale_y_continuous(limits = c(0,100))

# Create a new data frame with ages from 0 to 25
age_predictions <- data.frame(age_age = seq(0.001, 31, length = 100))

# Generate predictions
pred <- predict(model, newdata = age_predictions, type = "link", se.fit = TRUE)
age_predictions$mean <- plogis(pred$fit)
age_predictions$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
age_predictions$upper <- plogis(pred$fit + 1.96 * pred$se.fit)

# Round to 3 decimal places for cleaner output
age_predictions$mean <- round(age_predictions$mean, 3)
age_predictions$lower <- round(age_predictions$lower, 3)
age_predictions$upper <- round(age_predictions$upper, 3)

# Save to CSV
write.csv(age_predictions, "output/maturity_proportions.csv", row.names = FALSE)

# Print confirmation and first few rows
cat("Saved maturity proportions for ages 0-25 to maturity_proportions.csv\n")
head(age_predictions)
