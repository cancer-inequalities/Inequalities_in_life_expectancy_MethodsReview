# Life tables and life expectancy estimations using regression-based method: ----

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

# Load datasets

educ_females <- readRDS(glue('{wd}/Data_folder/educ_females.rds'))
educ_males <- readRDS(glue('{wd}/Data_folder/educ_males.rds'))

pop_deaths_2017 <- readRDS(glue('{wd}/Out/pop_deaths_2017.rds'))

# Load functions

source("Scripts/5.1. Regression_based_fxs.R")

# Some transformations of datasets

deaths <- pop_deaths_2017 %>% 
  select(sex, age, years_ed, n_deaths)

educ_females <- educ_females %>% 
  mutate(sex = rep(2,length(year)))

educ_males <- educ_males %>% 
  mutate(sex = rep(1,length(year)))

educ <- educ_males %>% 
  bind_rows(educ_females)

educ <- educ %>% 
  mutate(years_ed = case_when(
    years_ed >= 1 & years_ed <= 6 ~ 1,
    years_ed >= 7 & years_ed <= 12 ~ 2,
    years_ed >= 13 & years_ed <= 16 ~ 3,
    years_ed >= 17 & years_ed <= 20 ~ 4,
    TRUE ~ 0
  )) %>% 
  group_by(sex, year, years_ed) %>% 
  summarise(p = sum(p))

# Full combinations dataset

sexes <- 1:2
ages <- 26:100
years_ed <- 0:4

full_combinations <- expand.grid(sex = sexes,
                                       age = ages,
                                       years_ed = years_ed)

# Period regression based: input data ----

# For this approach we will estimate the proportion of individuals by educational level and sex (not adjusting by age)

numerator_period <- pop_deaths_2017 %>% 
  group_by(sex, years_ed) %>% 
  summarize(n_census = sum(n_census, na.rm = TRUE), .groups = "drop")

denominator_period <- pop_deaths_2017 %>% 
  group_by(sex) %>% 
  summarize(pop_census = sum(n_census, na.rm = TRUE), .groups = "drop")

period_proportions <- pop_deaths_2017 %>% 
  select(-c(n_census, n_deaths)) %>% 
  left_join(denominator_period, by = "sex") %>% 
  left_join(numerator_period, by = c("sex", "years_ed")) %>% 
  mutate(p = n_census/pop_census)

# Then, we will estimate the cumulative proportion (cd), by education level (years_educ) and age

period_proportions <- period_proportions %>% 
  arrange(sex, age, years_ed) %>% 
  group_by(sex, age) %>% 
  mutate(cd = cumsum(p))

# Now, we need to estimate the midpoint of cumulative proportion (pc)

period_proportions <- period_proportions %>% 
  group_by(sex, age) %>% 
  mutate(pc = ifelse(row_number() == 1, cd/2, (cd + lag(cd))/2))

# The next sterp is to estimate the mortality ratio (just adjusted by sex). First we need the denominator: 

denominator_rm_period <- deaths %>% 
  group_by(sex) %>% 
  summarize(pop_deaths = sum(n_deaths), .groups = 'drop')

# Join denominator

period_proportions <- period_proportions %>% 
  left_join(denominator_rm_period, by = "sex")

# Obtain n_deaths

deaths_period <- deaths %>% 
  group_by(sex, years_ed) %>% 
  summarize(n_deaths = sum(n_deaths))

period_proportions <- period_proportions %>% 
  left_join(deaths_period, by = c("sex", "years_ed"))

period_proportions <- period_proportions %>% 
  mutate(rm = (n_deaths/n_census)/(pop_deaths/pop_census))

# Cohort regression based: input data ----

# For the cohort method, the proportion of individuals by years of education, age, and sex corresponds to the observed distribution in a given year. 
# To compute this, we first obtain the total number of individuals by age and sex (denominator), 
# and the number of individuals by years of education, age, and sex (numerator).

cohort_proportions <- pop_deaths_2017 %>% 
  group_by(sex, age) %>% 
  summarise(pop_census = sum(n_census, na.rm = TRUE), .groups = "drop")

cohort_proportions <- pop_deaths_2017 %>%
  left_join(cohort_proportions, by = c("sex", "age"))

cohort_proportions <- cohort_proportions %>% 
  mutate(p = n_census/pop_census) 

# Once we obtained the proportion, we'll need to estimate the cumulative proportion (cd), by educational level (and adjusting by age)

cohort_proportions <- cohort_proportions %>%
  arrange(sex, age, years_ed) %>%
  group_by(sex, age) %>%
  mutate(cd = cumsum(p)) %>%
  ungroup()

# Then, we have to estimate the midpoint of cumulative proportion (pc) for each educational level. 

cohort_proportions <- cohort_proportions %>% 
  group_by(sex,age) %>% 
  mutate(pc = ifelse(row_number() == 1, cd/2, (cd + lag(cd)) / 2))

# To estimate the mortality ratio "rm" (adjusted by sex and age). First, we'll obtain the denominator

denominator_rm_both_cohorts <- deaths %>%
  group_by(sex, age) %>%
  summarize(pop_deaths = sum(n_deaths), .groups = 'drop')

# Join denominator

cohort_proportions <- cohort_proportions %>% 
  left_join(denominator_rm_both_cohorts, by = c("sex", "age"))

cohort_proportions <- cohort_proportions %>% 
  mutate(rm = (n_deaths/n_census)/(pop_deaths/pop_census))

# Save for future data visualization

saveRDS(cohort_proportions, file = glue('{cnst$path_out}/cohort_proportions_for_visualization.rds'))

# Lagged cohort regression based: input data ----

# We need to obtain the proportion of individuals by educational level, sex, and age. In the case of the lagged cohort method,
# this information comes from the educ database (which includes educ_females and educ_males). 
# For a person aged x in 2017, the proportion (p) corresponds to the observed proportion by educational level and sex in 
# the year when that person was 26 years old, i.e., 2017 - x + 26.
# For example, if a woman is 80 years old in 2017 and has an educational level of 1, the corresponding proportion is that of a 
# 26-year-old woman in 1963 (2017 - 80 + 26), which is 0.02599576.  

lagged_cohort_proportions <- full_combinations %>%
  mutate(year = 2017 - age + 26)

lagged_cohort_proportions <- lagged_cohort_proportions %>%
  left_join(educ, by = c("sex", "years_ed", "year")) %>% 
  select(-year)

# Once these proportions are obtained, we then estimate the cumulative proportion (cd) by educational level, adjusted for age.

lagged_cohort_proportions <- lagged_cohort_proportions %>%
  arrange(sex, age, years_ed) %>%
  group_by(sex, age) %>%
  mutate(cd = cumsum(p)) %>% 
  ungroup()

# With the cumulative estimate available, we then estimate the midpoint of the cumulative proportion for each educational level 
# (adjusted for age), which we refer to as "pc" (percentile).

lagged_cohort_proportions <- lagged_cohort_proportions %>% 
  group_by(sex,age) %>% 
  mutate(pc = ifelse(row_number() == 1, cd/2, (cd + lag(cd)) / 2))

# We use the rm value calculated using the cohort method, as it remains the same regardless of the method applied.

lagged_cohort_proportions <- lagged_cohort_proportions %>% 
  left_join(cohort_proportions %>% select(sex, age, years_ed, rm), by = c("sex", "age", "years_ed"))

# Save for future data visualization

saveRDS(lagged_cohort_proportions, file = glue('{cnst$path_out}/lagged_cohort_proportions_for_visualization.rds'))

# Some important checks about rm (must be a 0 rows x n columns tibble)

filter(period_proportions, rm == 0 | is.na(rm))
filter(cohort_proportions, rm == 0 | is.na(rm))
filter(lagged_cohort_proportions, rm == 0 | is.na(rm))

# Regressions by approach

period_regressions <- regression(period_proportions)
cohort_regressions <- regression(cohort_proportions)
lagged_cohort_regressions <- regression(lagged_cohort_proportions)

# As the observed mortality rate (mx) is calculated within the life table function by sex and age,
# we first obtain the required inputs: the number of deaths (pop_deaths) and the population at risk (pop_census).

observed_pops <- pop_deaths_2017 %>% 
  group_by(sex, age) %>% 
  mutate(pop_deaths = sum(n_deaths),
         pop_census = sum(n_census)) %>% 
  select(sex, age, pop_deaths, pop_census) %>% 
  distinct()

# This observed mortality rate is then adjusted by multiplying it by the R factor, which is computed as a function of 
# the parameters a, b, alpha, and beta.
# The adjustment is applied for each combination of quintile, age, and sex in the cohort and lagged cohort methods,
# and by quintile and sex in the period method. The R factor is derived using Equation 3 from Hendi et al. (2021).
# The computation of the factor was implemented using both forms of Equation 3: the definite integral form and the closed-form 
# expression. This allowed us to confirm that both approaches yield equivalent results (see Regression_based_fxs.R).

# We define the values of parameters a and b as follows:

quintiles <- list(
  quintile1 = list(a = 0, b = 0.2),
  quintile2 = list(a = 0.2, b = 0.4),
  quintile3 = list(a = 0.4, b = 0.6),
  quintile4 = list(a = 0.6, b = 0.8),
  quintile5 = list(a = 0.8, b = 1)
)

# We apply the function to compute the R factor, iterating over each quintile, and store the results in a separate 
# dataset for each of the methods.

# Period

R_results_period <- bind_rows(
  period_regressions %>% rowwise() %>% do(processing(., "quintile1", quintiles$quintile1$a, quintiles$quintile1$b)),
  period_regressions %>% rowwise() %>% do(processing(., "quintile2", quintiles$quintile2$a, quintiles$quintile2$b)),
  period_regressions %>% rowwise() %>% do(processing(., "quintile3", quintiles$quintile3$a, quintiles$quintile3$b)),
  period_regressions %>% rowwise() %>% do(processing(., "quintile4", quintiles$quintile4$a, quintiles$quintile4$b)),
  period_regressions %>% rowwise() %>% do(processing(., "quintile5", quintiles$quintile5$a, quintiles$quintile5$b))
)

# Cohort
 
R_results_cohort <- bind_rows(
  cohort_regressions %>% rowwise() %>% do(processing(., "quintile1", quintiles$quintile1$a, quintiles$quintile1$b)),
  cohort_regressions %>% rowwise() %>% do(processing(., "quintile2", quintiles$quintile2$a, quintiles$quintile2$b)),
  cohort_regressions %>% rowwise() %>% do(processing(., "quintile3", quintiles$quintile3$a, quintiles$quintile3$b)),
  cohort_regressions %>% rowwise() %>% do(processing(., "quintile4", quintiles$quintile4$a, quintiles$quintile4$b)),
  cohort_regressions %>% rowwise() %>% do(processing(., "quintile5", quintiles$quintile5$a, quintiles$quintile5$b))
)

# Lagged cohort

R_results_lagged_cohort <- bind_rows(
  lagged_cohort_regressions %>% rowwise() %>% do(processing(., "quintile1", quintiles$quintile1$a, quintiles$quintile1$b)),
  lagged_cohort_regressions %>% rowwise() %>% do(processing(., "quintile2", quintiles$quintile2$a, quintiles$quintile2$b)),
  lagged_cohort_regressions %>% rowwise() %>% do(processing(., "quintile3", quintiles$quintile3$a, quintiles$quintile3$b)),
  lagged_cohort_regressions %>% rowwise() %>% do(processing(., "quintile4", quintiles$quintile4$a, quintiles$quintile4$b)),
  lagged_cohort_regressions %>% rowwise() %>% do(processing(., "quintile5", quintiles$quintile5$a, quintiles$quintile5$b))
)
  
# Next, we join datasets to get the final input 

lagged_cohort_mortality_by_quintile <- R_results_lagged_cohort %>%
  left_join(observed_pops, by = c("age", "sex"))

cohort_mortality_by_quintile <- R_results_cohort %>% 
  left_join(observed_pops, by = c("age", "sex")) 

period_mortality_by_quintile <-  R_results_period %>% 
  left_join(observed_pops, by = c("age", "sex"))


# Life tables ----

# Finally, we construct the life tables for each method by applying the obtain_lt function

## Cohort ----

cohort_lt <- obtain_lt(cohort_mortality_by_quintile)

## Lagged cohort ----

lagged_cohort_lt <- obtain_lt(lagged_cohort_mortality_by_quintile)

## Period ----

period_lt <- obtain_lt(period_mortality_by_quintile)

# Save dataset ----

saveRDS(cohort_lt, file = glue('{cnst$path_out}/life_table_regression_based_cohort.rds')) 
saveRDS(lagged_cohort_lt, file = glue('{cnst$path_out}/life_table_regression_based_lagged_cohort.rds'))
saveRDS(period_lt, file = glue('{cnst$path_out}/life_table_regression_based_period.rds'))

rm(list = ls())
