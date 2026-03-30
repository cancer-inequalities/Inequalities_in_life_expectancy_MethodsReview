# Life tables and life expectancy estimations by compositional adjustment using simulated mortality rates: ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(haven, dplyr, purrr, rlang, glue, here, tidyr)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Functions ----

# Function to create the life table input by quintile

create_input <- function(data, pop_data, deaths_var, pop_var) {
  
  data <- data %>% 
    left_join(pop_data, by = c("age", "sex", "years_ed")) %>% 
    mutate(
      denomq1 = round({{ pop_var }} * wq1, 0),
      denomq2 = round({{ pop_var }} * wq2, 0),
      denomq3 = round({{ pop_var }} * wq3, 0),
      denomq4 = round({{ pop_var }} * wq4, 0),
      denomq5 = round({{ pop_var }} * wq5, 0),
      numq1   = round({{ deaths_var }} * wq1, 0),
      numq2   = round({{ deaths_var }} * wq2, 0),
      numq3   = round({{ deaths_var }} * wq3, 0),
      numq4   = round({{ deaths_var }} * wq4, 0),
      numq5   = round({{ deaths_var }} * wq5, 0)
    )
  
  counts <- data %>% 
    group_by(age, sex) %>% 
    summarise(
      quintile1 = sum(denomq1, na.rm = TRUE),
      quintile2 = sum(denomq2, na.rm = TRUE),
      quintile3 = sum(denomq3, na.rm = TRUE),
      quintile4 = sum(denomq4, na.rm = TRUE),
      quintile5 = sum(denomq5, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    pivot_longer(
      cols = starts_with("quintile"), 
      names_to = "quintile", 
      values_to = "counts"
    )
  
  deaths <- data %>% 
    group_by(age, sex) %>%
    summarise(
      quintile1 = sum(numq1, na.rm = TRUE),
      quintile2 = sum(numq2, na.rm = TRUE),
      quintile3 = sum(numq3, na.rm = TRUE),
      quintile4 = sum(numq4, na.rm = TRUE),
      quintile5 = sum(numq5, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    pivot_longer(
      cols = starts_with("quintile"), 
      names_to = "quintile", 
      values_to = "total_deaths"
    )
  
  final <- counts %>% 
    left_join(deaths, by = c("sex", "age", "quintile"))
  
  return(final)
}

# Calculate life tables function

CalculateLifeTable <-
  function (df, x, nx = c(rep(1,74),Inf), Dx, Ex) {
    
    require(dplyr)
    
    df %>%
      transmute(
        x = {{x}},
        nx = {{nx}},
        mx = {{Dx}}/{{Ex}},
        px = exp(-mx*{{nx}}),
        qx = 1-px,
        lx = head(cumprod(c(1, px)), -1),
        dx = c(-diff(lx), tail(lx, 1)),
        Lx = ifelse(mx==0, lx*nx, dx/mx),
        Tx = rev(cumsum(rev(Lx))),
        ex = Tx/lx,
        sd0 =  c(sqrt(sum(dx*(x+.5-ex[1])^2)), rep(0,length(x)-1))
      )
  }

# Calculate central estimates of life-expectancy

le_central_estimates <- function(data){
  
  central_estimate <- data %>% 
    arrange(quintile, sex, age) %>% 
    group_by(quintile, sex) %>% 
    group_modify(~{
      CalculateLifeTable(df = .x, x = age, Dx = total_deaths, Ex = counts)
    }) %>% 
    ungroup()
  
  return(central_estimate)
  
}

# Calculate confidence intervals for life expectancy

confidence_intervals <- function(data){
  
  # Simulations to obtain confidence intervals
  
  set.seed(1987)
  
  simulations <- data %>% 
    expand_grid(id_sim = 1:500) %>% # we can modify the number of simulations
    group_by(quintile, sex, age) %>% 
    mutate(death_total_sim = rpois(500, total_deaths)) %>%
    arrange(age,quintile, sex) %>%
    group_by(id_sim, quintile, sex) %>%
    group_modify(~{
      CalculateLifeTable(df =.x, x = age, Dx = death_total_sim, Ex = counts)
    }) %>%
    ungroup()
  
  # Assemble table with ex statistics
  
  confidence_intervals <- simulations %>% 
    select(id_sim, quintile, sex, x, mx, ex, sd0) %>%
    arrange(id_sim, quintile, sex, x) %>%
    group_by(quintile, sex, x) %>%
    summarise(
      ex_q025 = quantile(ex, 0.025, na.rm = TRUE),
      ex_q975 = quantile(ex, 0.975, na.rm = TRUE),
      sd0_q025 = quantile(sd0, 0.025, na.rm = TRUE),
      sd0_q975 = quantile(sd0, 0.975, na.rm = TRUE)
    )
  
  return(confidence_intervals)
}

# A). Asymmetric simulations ----

# Load data ---- 

simulations_asy <- read_dta(glue('{wd}/Out/simulated_mortality_data.dta'))

# Distribution of deaths and population across educational quintiles ----

# Load weight data generated in Script 2. Education_distribution.R

weights <- readRDS(glue('{cnst$path_out}/quintile_weights_lagged_cohort_by_sex.rds'))

# Quintile weights by sex, age, and years of education (years_ed) were
# generated in Script 2. Education_distribution.R. These weights describe
# how the population within each sex–age–education stratum is distributed
# along the educational rank.

# To construct mortality schedules by educational quintile, both the
# population at risk and the number of deaths must be allocated across
# quintiles. Because the available data provide mortality rates rather
# than death counts, we first reconstruct the expected number of deaths
# by combining the observed population counts with the reported rates.

# The resulting population exposures and deaths are then redistributed
# across quintiles using the corresponding quintile weights.

simulations_asy <- simulations_asy %>%
  mutate(
    deaths_loglin = rate_sim_loglin * n_census / 100000,
    deaths_shape = rate_sim_shape * n_census / 100000
  )

# Next, we construct the life table inputs, carefully selecting the relevant
# death counts and population exposures

input_loglin_asy <- create_input(weights, simulations_asy, deaths_var = deaths_loglin, pop_var = n_census)

input_shape_asy <- create_input(weights, simulations_asy, deaths_var = deaths_shape, pop_var = n_census)

# Once the inputs are defined, we apply the functions that sequentially generate
# the life tables

# Obtain central estimates of life-expectancy ----

le_central_estimate_loglin_asy <- le_central_estimates(input_loglin_asy)
le_central_estimate_shape_asy <- le_central_estimates(input_shape_asy)

# Add some statistics ----

ci_loglin_asy <- confidence_intervals(input_loglin_asy)
ci_shape_asy <- confidence_intervals(input_shape_asy)

# Assemble final table ----

# Assemble all the ex statistics in a single table

life_table_loglin_asy <-
  left_join(
    le_central_estimate_loglin_asy,
    ci_loglin_asy
  )

life_table_shape_asy <- 
  left_join(
    le_central_estimate_shape_asy,
    ci_shape_asy
  )

# Save final tables ----

saveRDS(life_table_loglin_asy, file = glue('{cnst$path_out}/life_table_compositional_adjustment_sim_loglin_asymmetric.rds'))
saveRDS(life_table_shape_asy, file = glue('{cnst$path_out}/life_table_compositional_adjustment_sim_shape_asymmetric.rds'))

rm(list = setdiff(ls(), c("wd", "cnst", "create_input", "CalculateLifeTable", "le_central_estimates",
                          "confidence_intervals")))

# B). Symmetric simulations ----

# Load data ---- 

simulations_sym <- read_dta(glue('{wd}/Out/simulated_mortality_data_simmetric.dta'))

# Distribution of deaths and population across educational quintiles ----

# Construct equal weights to impose a balanced quintile distribution (20% per group)

equal_quintile_weights <- simulations_sym %>%
  distinct(age, sex, years_ed) %>%
  mutate(
    wq1 = if_else(years_ed == 0, 1, 0),
    wq2 = if_else(years_ed == 1, 1, 0),
    wq3 = if_else(years_ed == 2, 1, 0),
    wq4 = if_else(years_ed == 3, 1, 0),
    wq5 = if_else(years_ed == 4, 1, 0)
  ) %>%
  arrange(sex, age, years_ed)

# Obtain population at risk ----

population_at_risk <- simulations_sym %>% 
  group_by(age, sex) %>% 
  summarise(pop_at_risk = sum(n_census)/5)

# Join population at risk and simulated data ----

simulations_sym <- simulations_sym %>% 
  left_join(population_at_risk, by = c("age", "sex"))

# To construct mortality schedules by educational quintile, both the
# population at risk and the number of deaths must be allocated across
# quintiles. Because the available data provide mortality rates rather
# than death counts, we first reconstruct the expected number of deaths
# by combining the observed population counts with the reported rates.

# The resulting population exposures and deaths are then redistributed
# across quintiles using the corresponding quintile weights.

simulations_sym <- simulations_sym %>%
  mutate(
    deaths_loglin = rate_sim_loglin * pop_at_risk / 100000,
    deaths_shape = rate_sim_shape * pop_at_risk / 100000
  )

# Next, we construct the life table inputs, carefully selecting the relevant
# death counts and population exposures.

input_loglin_sym <- create_input(equal_quintile_weights, simulations_sym, deaths_var = deaths_loglin, pop_var = pop_at_risk)

input_shape_sym <- create_input(equal_quintile_weights, simulations_sym, deaths_var = deaths_shape, pop_var = pop_at_risk)

# Once the inputs are defined, we apply the functions that sequentially generate
# the life tables

# Obtain central estimates of life-expectancy ----

le_central_estimate_loglin_sym <- le_central_estimates(input_loglin_sym)
le_central_estimate_shape_sym <- le_central_estimates(input_shape_sym)

# Add some statistics ----

ci_loglin_sym <- confidence_intervals(input_loglin_sym)
ci_shape_sym <- confidence_intervals(input_shape_sym)

# Assemble final table ----

# Assemble all the ex statistics in a single table

life_table_loglin_sym <-
  left_join(
    le_central_estimate_loglin_sym,
    ci_loglin_sym
  )

life_table_shape_sym <- 
  left_join(
    le_central_estimate_shape_sym,
    ci_shape_sym
  )

# Save final tables ----

saveRDS(life_table_loglin_sym, file = glue('{cnst$path_out}/life_table_compositional_adjustment_sim_loglin_symmetric.rds'))
saveRDS(life_table_shape_sym, file = glue('{cnst$path_out}/life_table_compositional_adjustment_sim_shape_symmetric.rds'))

rm(list = ls())
