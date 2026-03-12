# Helper functions for education distribution: ----

# General function to assign quintile weights for a given vector of proportions ----

assign_quintile_weights <- function(p_vec) {
  k <- length(p_vec)
  quintiles <- seq(0, 1, by = 0.2)
  cum_p <- c(0, cumsum(p_vec))
  
  w_mat <- matrix(0, nrow = k, ncol = 5)
  colnames(w_mat) <- paste0("wq", 1:5)
  rownames(w_mat) <- paste0("ed_", 0:(k - 1))  # assuming categories 0 to 4
  
  for (s in 1:k) {
    c_prev <- cum_p[s]
    c_curr <- cum_p[s + 1]
    p_s <- p_vec[s]
    
    for (q in 1:5) {
      lower_q <- quintiles[q]
      upper_q <- quintiles[q + 1]
      
      overlap_lower <- max(c_prev, lower_q)
      overlap_upper <- min(c_curr, upper_q)
      
      if (overlap_upper > overlap_lower) {
        w_mat[s, q] <- (overlap_upper - overlap_lower) / p_s
      }
    }
  }
  
  return(as.data.frame(w_mat))
}

# Get quintile weights for each approach ----

## Cohort approach ----

get_quintile_weights_for_cohort <- function(sex_value) {
  pop_deaths_2017 %>%
    filter(sex == sex_value) %>%
    group_by(age) %>%
    mutate(total_age_n_census = sum(n_census)) %>%
    group_by(age, years_ed) %>%
    mutate(p = n_census / total_age_n_census) %>%
    ungroup() %>%
    group_by(age) %>%
    arrange(age, years_ed) %>%
    summarise(across(c(years_ed, p), list), .groups = "drop") %>%
    mutate(weights = map(p, assign_quintile_weights)) %>%
    unnest(cols = c(years_ed, p, weights)) %>%
    mutate(sex = sex_value) %>% 
    select(-p)
}

## Lagged cohort approach ----

get_quintile_weights_for_lagged_cohort <- function(pop_data, educ_data, sex_value) {
  
  # recategorizar los niveles de educación 
  
  educ_data <- educ_data %>% 
    mutate(years_ed = case_when(
      years_ed == 0 ~ 0,
      years_ed >= 1 & years_ed <= 6 ~ 1,
      years_ed >= 7 & years_ed <= 12 ~ 2,
      years_ed >= 13 & years_ed <= 16 ~ 3,
      years_ed >= 17 & years_ed <= 20 ~ 4,
      TRUE ~ NA_real_
    )) %>%
    group_by(year, years_ed) %>%  # Group by year and new years_ed categories
    summarise(p = sum(p, na.rm = TRUE), .groups = "drop")  # Sum p within each group
  
  # Filtrar población del sexo correspondiente y calcular el año en que tuvo 26 años
  
  pop_data_sex <- pop_data %>%
    filter(sex == sex_value) %>%
    mutate(year_at_26 = 2017 - age + 26)  # Año en que se estima tenía 26 años
  
  # Añadir proporciones desde educ_data
  
  pop_data_joined <- pop_data_sex %>%
    left_join(
      educ_data %>% rename(year_at_26 = year),
      by = c("year_at_26", "years_ed")
    )
  
  # Aplicar función de pesos por edad
  
  weights_by_age <- pop_data_joined %>%
    group_by(age) %>%
    arrange(age, years_ed) %>%
    summarise(across(c(years_ed, p), list), .groups = "drop") %>%
    mutate(weights = map(p, assign_quintile_weights)) %>%
    unnest(cols = c(years_ed, p, weights)) %>%
    mutate(sex = sex_value) %>% 
    select(-p)
  
  return(weights_by_age)
}

## Period approach ----

get_quintile_weights_for_period <- function(pop_data, sex_value){
  
  pop_data_sex <- pop_data %>% 
    filter(sex == sex_value)
  
  #Compute total population by education level
  
  props <- pop_data_sex %>% 
    group_by(years_ed) %>% 
    summarise(p = sum(n_census, na.rm = TRUE)) %>% 
    mutate(p = p/sum(p))
  
  #Join the proportions back to the original data

  pop_data_sex <- pop_data_sex %>%
  left_join(props, by = "years_ed")

  # Apply function by age
  
  weights_by_age <- pop_data_sex %>%
    group_by(age) %>%
    arrange(age, years_ed) %>%
    summarise(across(c(years_ed, p), list), .groups = "drop") %>%
    mutate(weights = map(p, assign_quintile_weights)) %>%
    unnest(cols = c(years_ed, p, weights)) %>% 
    mutate(sex = sex_value) %>%
    select(-p)
  
  return(weights_by_age)
  
}





