# Life tables and life expectancy estimations by compositional adjustment: ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(here, glue, dplyr, tidyr, purrr, haven)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Load helper functions ----

source("Scripts/4.1. Life_tables_fxs.R")

# We need to load pop_deaths_2017 dataset because we'll need it
# for denominators, using n_census variable

pop_deaths_2017 <- readRDS(glue('{wd}/Out/pop_deaths_2017.rds'))

# This uses census datasets with quintiles weights by age, sex and years of education created in
# Education_distribution.R script

# Period compositional adjustment ----

weights_period <- readRDS(glue('{cnst$path_out}/quintile_weights_period_by_sex.rds'))

# Create the imput to life tables

input_life_table_period <- create_input(weights_period, pop_deaths_2017)

# Obtain central estimates of life-expectancy

le_central_estimate_period <- le_central_estimates(input_life_table_period)

# Add some statistics

ci_period <- confidence_intervals(input_life_table_period)

# Assemble final table

# Assemble all the ex statistics in a single table

life_table_period <-
  left_join(
    le_central_estimate_period,
    ci_period
  )

# Save final table

saveRDS(life_table_period, file = glue('{cnst$path_out}/life_table_compositional_adjustment_period_approach.rds'))

# Cohort compositional adjustment ----

weights_cohort <- readRDS(glue('{cnst$path_out}/quintile_weights_cohort_by_sex.rds'))

# Create the imput to life tables

input_life_table_cohort <- create_input(weights_cohort, pop_deaths_2017)

# Obtain central estimates of life-expectancy

le_central_estimate_cohort <- le_central_estimates(input_life_table_cohort)

# Add some statistics

ci_cohort <- confidence_intervals(input_life_table_cohort)

# Assemble final table

# Assemble all the ex statistics in a single table

life_table_cohort <-
  left_join(
    le_central_estimate_cohort,
    ci_cohort
  )

# Save final table

saveRDS(life_table_cohort, file = glue('{cnst$path_out}/life_table_compositional_adjustment_cohort_approach.rds'))

# Lagged cohort compositional adjustment ----

weights_lagged_cohort <- readRDS(glue('{cnst$path_out}/quintile_weights_lagged_cohort_by_sex.rds'))

# Create the imput to life tables

input_life_table_lagged_cohort <- create_input(weights_lagged_cohort, pop_deaths_2017)

# Obtain central estimates of life-expectancy

le_central_estimate_lagged_cohort <- le_central_estimates(input_life_table_lagged_cohort)

# Add some statistics

ci_lagged_cohort <- confidence_intervals(input_life_table_lagged_cohort)

# Assemble final table

# Assemble all the ex statistics in a single table

life_table_lagged_cohort <-
  left_join(
    le_central_estimate_lagged_cohort,
    ci_lagged_cohort
  )

# Save final table

saveRDS(life_table_lagged_cohort, file = glue('{cnst$path_out}/life_table_compositional_adjustment_lagged_cohort_approach.rds'))

# Level approach ----

input_life_table_level <- pop_deaths_2017 %>%
  mutate(quintile = case_when(
    years_ed == 0 ~ "quintile1",
    years_ed == 1 ~ "quintile2",
    years_ed == 2 ~ "quintile3",
    years_ed == 3 ~ "quintile4",
    years_ed == 4 ~ "quintile5",
    TRUE ~ NA 
  )) %>% 
  rename(total_deaths = n_deaths,
         counts = n_census) %>% 
  select(-years_ed)

le_central_estimate_level <- le_central_estimates(input_life_table_level)

# Add some statistics

ci_level <- confidence_intervals(input_life_table_level)

# Assemble final table

# Assemble all the ex statistics in a single table

life_table_level <-
  left_join(
    le_central_estimate_level,
    ci_level
  )

saveRDS(life_table_level, file = glue('{cnst$path_out}/life_table_compositional_adjustment_level_approach.rds') )

rm(list = ls())
