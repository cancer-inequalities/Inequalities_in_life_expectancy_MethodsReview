# Estimation of quintile weights by approach: ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(here, glue, dplyr, tidyr, purrr)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Load helper functions ----

source("Scripts/2.1. Education_distribution_fxs.R")

# Load census data ---- 

pop_deaths_2017 <- readRDS(glue('{wd}/Out/pop_deaths_2017.rds'))

# Cohort ----

cohort_weights_by_age <- bind_rows(
  get_quintile_weights_for_cohort(1),
  get_quintile_weights_for_cohort(2)
)

saveRDS(cohort_weights_by_age, file = glue('{cnst$path_out}/quintile_weights_cohort_by_sex.rds'))

# Lagged cohort ----

educ_females <- readRDS(glue('{wd}/Data_folder/educ_females.rds'))
educ_males <- readRDS(glue('{wd}/Data_folder/educ_males.rds'))

lagged_cohort_weights_by_age <- bind_rows(
  get_quintile_weights_for_lagged_cohort(pop_deaths_2017, educ_males, sex_value = 1),
  get_quintile_weights_for_lagged_cohort(pop_deaths_2017, educ_females, sex_value = 2)
)

saveRDS(lagged_cohort_weights_by_age, file = glue("{cnst$path_out}/quintile_weights_lagged_cohort_by_sex.rds"))

# Period ----

period_weights_by_age <- bind_rows(
  get_quintile_weights_for_period(pop_deaths_2017, sex_value = 1),
  get_quintile_weights_for_period(pop_deaths_2017, sex_value = 2),
)

saveRDS(period_weights_by_age, file = glue("{cnst$path_out}/quintile_weights_period_by_sex.rds"))

rm(list = ls())
