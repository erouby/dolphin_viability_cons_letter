get_growth_rates <- function(data_surv = NULL, data_fec = NULL, n_year = NULL) {

for (y in 1:n_year) {
  
  table <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                      "lx" = filter(data_surv, method == "Stan" & year == y)$mean, # Cumulative survival rate
                      "fec_x" = filter(data_fec)$mean, # Proportion matures
                      "mx" = rep(NA),
                      "sx" = rep(NA)) # Age-Specific survival rate
  
  table_lower <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                            "lx" = filter(data_surv, method == "Stan" & year == y)$lower, # Cumulative survival rate
                            "fec_x" = filter(data_fec)$lower, # Proportion matures
                            "mx" = rep(NA),
                            "sx" = rep(NA)) # Age-Specific survival rate
  
  table_upper <- data.frame("age" = filter(data_surv, method == "Stan" & year == y)$age, # Age-class
                            "lx" = filter(data_surv, method == "Stan" & year == y)$upper, # Cumulative survival rate
                            "fec_x" = filter(data_fec)$upper, # Proportion matures
                            "mx" = rep(NA),
                            "sx" = rep(NA)) # Age-Specific survival rate
  
  for (i in 1:nrow(table)) {
    
    table$sx[i] <- table$lx[i+1]/table$lx[i]
    
    if (table$age[i] >= 0 & table$age[i] < 4 ) {
      table$mx[i] <- table$fec_x[i] * (0) * table$sx[i]
    }
    if (table$age[i] >= 4 & table$age[i] < 9 ) {
      table$mx[i] <- table$fec_x[i] * (0.125) * table$sx[i]
    }
    if (table$age[i] >= 9 & table$age[i] < 13 ) {
      table$mx[i] <- table$fec_x[i] * (0.40) * table$sx[i]
    }
    if (table$age[i] >= 13 & table$age[i] < 18 ) {
      table$mx[i] <- table$fec_x[i] * (0.35) * table$sx[i]
    }
    if (table$age[i] >= 18) {
      table$mx[i] <- table$fec_x[i] * (0.30) * table$sx[i]
    }
  }
  
  for (i in 1:nrow(table_lower)) {
    
    table_lower$sx[i] <- table_lower$lx[i+1]/table_lower$lx[i]
    
    if (table_lower$age[i] >= 0 & table_lower$age[i] < 4 ) {
      table_lower$mx[i] <- table_lower$fec_x[i] * (0) * table_lower$sx[i]
    }
    if (table_lower$age[i] >= 4 & table_lower$age[i] < 9 ) {
      table_lower$mx[i] <- table_lower$fec_x[i] * (0.125) * table_lower$sx[i]
    }
    if (table_lower$age[i] >= 9 & table_lower$age[i] < 13 ) {
      table_lower$mx[i] <- table_lower$fec_x[i] * (0.40) * table_lower$sx[i]
    }
    if (table_lower$age[i] >= 13 & table_lower$age[i] < 18 ) {
      table_lower$mx[i] <- table_lower$fec_x[i] * (0.35) * table_lower$sx[i]
    }
    if (table_lower$age[i] >= 18) {
      table_lower$mx[i] <- table_lower$fec_x[i] * (0.30) * table_lower$sx[i]
    }
  }
  
  for (i in 1:nrow(table_upper)) {
    
    table_upper$sx[i] <- table_upper$lx[i+1]/table_upper$lx[i]
    
    if (table_upper$age[i] >= 0 & table_upper$age[i] < 4 ) {
      table_upper$mx[i] <- table_upper$fec_x[i] * (0) * table_upper$sx[i]
    }
    if (table_upper$age[i] >= 4 & table_upper$age[i] < 9 ) {
      table_upper$mx[i] <- table_upper$fec_x[i] * (0.125) * table_upper$sx[i]
    }
    if (table_upper$age[i] >= 9 & table_upper$age[i] < 13 ) {
      table_upper$mx[i] <- table_upper$fec_x[i] * (0.40) * table_upper$sx[i]
    }
    if (table_upper$age[i] >= 13 & table_upper$age[i] < 18 ) {
      table_upper$mx[i] <- table_upper$fec_x[i] * (0.35) * table_upper$sx[i]
    }
    if (table_upper$age[i] >= 18) {
      table_upper$mx[i] <- table_upper$fec_x[i] * (0.30) * table_upper$sx[i]
    }
  }
  
  matrix <- leslie.matrix(lx = table$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                          mx = table$mx, # age-specific fertility rates
                          L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                          peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                          one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                          infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
  )
  
  matrix_lower <- leslie.matrix(lx = table_lower$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                                mx = table_lower$mx, # age-specific fertility rates
                                L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                                peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                                one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                                infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
  )
  
  matrix_upper <- leslie.matrix(lx = table_upper$lx, # vector of either age-specific cumulative survival or person-years lived in the interval
                                mx = table_upper$mx, # age-specific fertility rates
                                L = FALSE, # if ’FALSE’, lx is taken to be cumulative survival to exact age x + n
                                peryear = 1, # Multiplier for fertility. Defaults to peryear = 5
                                one.sex = TRUE, # FALSE, # if ’TRUE’, fertility rates will be divided by (1+SRB)
                                infant.class = TRUE # ’TRUE’ if lx contains a value for the infant age-class
  )
  
  l <- eigen.analysis(matrix)$lambda
  l_low <- eigen.analysis(matrix_lower)$lambda
  l_upp <- eigen.analysis(matrix_upper)$lambda
  rho <- eigen.analysis(matrix)$rho
  
  proj_dd <- data.frame("Percentage" = l^(0:time_serie),
                        "Time" = seq(0,100,1),
                        "Cohort" = rep(y),
                        "Lambda" = l,
                        "Lambda_low" = l_low,
                        "Lambda_upp" = l_upp,
                        "Rho" = rho)
  
  proj_years <- rbind(proj_years, proj_dd)
  
}

table_lambda <- proj_years %>% 
  select(- Percentage, - Time) %>%
  unique()

return(table_lambda)

}
