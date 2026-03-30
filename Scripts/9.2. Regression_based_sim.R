# Life tables and life expectancy estimations by regression-based method using simulated mortality rates: ----

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

# Regression function 

regression <- function(data, rm_var) {
  
  # Create a data frame to store the results
  
  rm_var <- as_name(ensym(rm_var))
  
  regressions <- data %>%
    group_by(sex, age) %>%
    reframe({
      modelo <- lm(
        formula = log(.data[[rm_var]]) ~ ridit,
        data = pick(everything())
      )
      
      data.frame(
        alpha = unname(coef(modelo)[1]),
        beta  = unname(coef(modelo)[2])
      )
    })
  
  return(regressions)
}

# We define the values of parameters a and b as follows:

quintiles_df <- tibble(
  quintile = paste0("quintile", 1:5),
  a = c(0, 0.2, 0.4, 0.6, 0.8),
  b = c(0.2, 0.4, 0.6, 0.8, 1)
)

# Compute the average mortality ratio by quintile following the integral-based method in Hendi et al. (2021)

R_factor <- function(alpha, beta, a, b) {
  R_factor <- function(x) {
    exp(alpha + beta * x)
  }
  result <- (1 / (b - a)) * integrate(R_factor, lower = a, upper = b)$value
  return(result)
}

# An alternative but equivalent expression is included as a consistency check.

second_formula <- function(alpha, beta, a, b) {
  result <- (exp(alpha) / (beta * (b - a))) * (exp(beta * b) - exp(beta * a))
  return(result)
}

compute_R_results <- function(regressions_df, quintiles_df) {
  
  regressions_df %>%
    crossing(quintiles_df) %>%
    mutate(
      R_factor = pmap_dbl(
        list(alpha, beta, a, b),
        \(alpha, beta, a, b) R_factor(alpha, beta, a, b)
      ),
      second_formula = pmap_dbl(
        list(alpha, beta, a, b),
        \(alpha, beta, a, b) second_formula(alpha, beta, a, b)
      )
    )
}

CalculateLifeTable <-
  function (df, x, nx = c(rep(1,74),Inf), Dx, Ex) {
    
    require(dplyr)
    
    df %>%
      transmute(
        x = {{x}},
        nx = {{nx}},
        mx = ({{Dx}}/{{Ex}}) * R_factor, 
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

obtain_lt <- function(data){
  
  set.seed(1987)
  
  central_estimates <-
    data %>%
    arrange(quintile, sex, age) %>%
    group_by(quintile, sex) %>%
    group_modify(~ {
      CalculateLifeTable(df = .x, x = age, Dx = pop_deaths, Ex = pop_census)
    }) %>%
    ungroup()
  
  simulations <-
    data %>%
    expand_grid(id_sim = 1:500) %>%
    group_by(quintile, sex, age) %>%
    mutate(death_total_sim = rpois(500, pop_deaths)) %>%
    arrange(age, quintile, sex) %>%
    group_by(id_sim, quintile, sex) %>%
    group_modify(~ {
      CalculateLifeTable(df = .x, x = age, Dx = death_total_sim, Ex = pop_census)
    }) %>%
    ungroup()
  
  statistics <-
    simulations %>%
    select(id_sim, quintile, sex, x, mx, ex, sd0) %>%
    arrange(id_sim, quintile, sex, x) %>%
    group_by(quintile, sex, x) %>%
    summarise(
      ex_q025 = quantile(ex, 0.025, na.rm = TRUE),
      ex_q975 = quantile(ex, 0.975, na.rm = TRUE),
      sd0_q025 = quantile(sd0, 0.025, na.rm = TRUE),
      sd0_q975 = quantile(sd0, 0.975, na.rm = TRUE)
    )
  
  final_set <-
    left_join(
      central_estimates,
      statistics
    )
  
  return(final_set)
}

# A). Asymmetric simulations ----

# Load data ---- 

simulations_asy <- read_dta(glue('{wd}/Out/simulated_mortality_data.dta'))

# Because the available data provide mortality rates rather than death counts, we first 
# reconstruct the expected number of deaths by combining the observed population counts with the reported rates. 

simulations_asy <- simulations_asy %>%
  mutate(
    deaths_loglin = rate_sim_loglin * n_census / 100000,
    deaths_shape = rate_sim_shape * n_census / 100000
  )

# Keep a data set containg information about expected number of deaths

deaths_loglin_asy <- simulations_asy %>% 
  select(sex, age, years_ed, deaths_loglin)

deaths_shape_asy <- simulations_asy %>% 
  select(sex, age, years_ed, deaths_shape)

# Cohort regression based: input data ----

# For the cohort method, the proportion of individuals by years of education, age, and sex corresponds to the observed distribution 
# in a given year. To compute this, we first obtain the total number of individuals by age and sex (denominator), 
# and the number of individuals by years of education, age, and sex (numerator).

# To estimate the mortality ratio "rm" (adjusted by sex and age). First, we'll obtain the denominator

denominator_rm_loglin_asy <- deaths_loglin_asy %>%
  group_by(sex, age) %>%
  summarize(pop_deaths_loglin = sum(deaths_loglin), .groups = 'drop')

denominator_rm_shape_asy <- deaths_shape_asy %>% 
  group_by(sex, age) %>% 
  summarise(pop_deaths_shape = sum(deaths_shape), .groups = 'drop')

# Join denominator

simulations_asy <- simulations_asy %>% 
  left_join(denominator_rm_loglin_asy, by = c("sex", "age")) %>% 
  left_join(denominator_rm_shape_asy, by = c("sex", "age"))

# Mortality ratio ----

simulations_asy <- simulations_asy %>% 
  left_join(simulations_asy %>% 
              group_by(sex, age) %>% 
              summarise(pop_census = sum(n_census, na.rm = TRUE), .groups = "drop")) %>% 
  mutate(rm_loglin = (deaths_loglin/n_census)/(pop_deaths_loglin/pop_census),
         rm_shape = (deaths_shape/n_census)/(pop_deaths_shape/pop_census))

# Some important checks about rm (must be a 0 rows x n columns tibble)

filter(simulations_asy, rm_loglin == 0 | is.na(rm_loglin))
filter(simulations_asy, rm_shape == 0 | is.na(rm_shape))

# Regressions by approach ----

regressions_loglin_asy <- regression(simulations_asy, rm_var = rm_loglin)
regressions_shape_asy <- regression(simulations_asy, rm_var = rm_shape)

# As the observed mortality rate (mx) is calculated within the life table function by sex and age,
# we first obtain the required inputs: the number of deaths (pop_deaths) and the population at risk (pop_census).

observed_pops_loglin_asy <- simulations_asy %>% 
  group_by(sex, age) %>% 
  mutate(pop_deaths = sum(deaths_loglin),
         pop_census = sum(n_census)) %>% 
  select(sex, age, pop_deaths, pop_census) %>% 
  distinct()

observed_pops_shape_asy <- simulations_asy %>% 
  group_by(sex, age) %>% 
  mutate(pop_deaths = sum(deaths_shape),
         pop_census = sum(n_census)) %>% 
  select(sex, age, pop_deaths, pop_census) %>% 
  distinct()

# This observed mortality rate is then adjusted by multiplying it by the R factor (average mortality ratio by quintile), which is 
# computed as a function of the parameters a, b, alpha, and beta.
# The adjustment is applied for each combination of quintile, age, and sex in the cohort and lagged cohort methods,
# and by quintile and sex in the period method. The R factor is derived using Equation 3 from Hendi et al. (2021).
# The computation of the factor was implemented using both forms of Equation 3: the definite integral form and the closed-form 
# expression. This allowed us to confirm that both approaches yield equivalent results (see Regression_based_fxs.R).

# We apply the function to compute the R factor, iterating over each quintile, and store the results in a separate 
# dataset for each of the methods.

R_results_loglin_asy <- compute_R_results(regressions_loglin_asy, quintiles_df)
R_results_shape_asy  <- compute_R_results(regressions_shape_asy, quintiles_df)

# Next, we join datasets to get the final input 

cohort_mortality_by_quintile_loglin_asy <- R_results_loglin_asy %>% 
  left_join(observed_pops_loglin_asy, by = c("age", "sex")) 

cohort_mortality_by_quintile_shape_asy <- R_results_shape_asy %>% 
  left_join(observed_pops_shape_asy, by = c("age", "sex")) 

# Life tables ----

# Finally, we construct the life tables for each method by applying the obtain_lt function

lt_loglin_asy <- obtain_lt(cohort_mortality_by_quintile_loglin_asy)
lt_shape_asy <- obtain_lt(cohort_mortality_by_quintile_shape_asy)

# Save life tables ----

saveRDS(lt_loglin_asy, file = glue('{cnst$path_out}/life_table_regression_based_sim_loglin_asymmetric.rds'))
saveRDS(lt_shape_asy, file = glue('{cnst$path_out}/life_table_regression_based_sim_shape_asymmetric.rds'))

rm(list = setdiff(ls(), c("wd", "cnst", "regression", "quintiles_df", "R_factor", "second_formula",
                          "compute_R_results", "CalculateLifeTable", "obtain_lt")))

# B). Symmetric simulations ----

# Load data ---- 

simulations_sym <- read_dta(glue('{wd}/Out/simulated_mortality_data_simmetric.dta'))

# To obtain population at risk 

population_at_risk <- simulations_sym %>% 
  group_by(age, sex) %>% 
  summarise(pop_at_risk = sum(n_census)/5)

# Join population at risk and simulations

simulations_sym <- simulations_sym %>% 
  left_join(population_at_risk, by = c("age", "sex"))

# Reconstruct the expected number of deaths 
# by combining the population at risk with the reported rates.

simulations_sym <- simulations_sym %>%
  mutate(
    deaths_loglin = rate_sim_loglin * pop_at_risk / 100000,
    deaths_shape = rate_sim_shape * pop_at_risk / 100000
  )

# Keep a data set containg information about expected number of deaths

deaths_loglin_sym <- simulations_sym %>% 
  select(sex, age, years_ed, deaths_loglin)

deaths_shape_sym <- simulations_sym %>% 
  select(sex, age, years_ed, deaths_shape)

# Cohort regression based: input data ----

# For the cohort method, the proportion of individuals by years of education, age, and sex corresponds to the observed distribution 
# in a given year. To compute this, we first obtain the total number of individuals by age and sex (denominator), 
# and the number of individuals by years of education, age, and sex (numerator).

# To estimate the mortality ratio "rm" (adjusted by sex and age). First, we'll obtain the denominator

denominator_rm_loglin_sym <- deaths_loglin_sym %>%
  group_by(sex, age) %>%
  summarize(pop_deaths_loglin = sum(deaths_loglin), .groups = 'drop')

denominator_rm_shape_sym <- deaths_shape_sym %>% 
  group_by(sex, age) %>% 
  summarise(pop_deaths_shape = sum(deaths_shape), .groups = 'drop')

# Join denominator

simulations_sym <- simulations_sym %>% 
  left_join(denominator_rm_loglin_sym, by = c("sex", "age")) %>% 
  left_join(denominator_rm_shape_sym, by = c("sex", "age"))

# Mortality ratio ----

simulations_sym <- simulations_sym %>% 
  left_join(simulations_sym %>% 
              group_by(sex, age) %>% 
              summarise(pop_census = sum(pop_at_risk, na.rm = TRUE), .groups = "drop")) %>% 
  mutate(rm_loglin = (deaths_loglin/pop_at_risk)/(pop_deaths_loglin/pop_census),
         rm_shape = (deaths_shape/pop_at_risk)/(pop_deaths_shape/pop_census))

# Some important checks about rm (must be a 0 rows x n columns tibble)

filter(simulations_sym, rm_loglin == 0 | is.na(rm_loglin))
filter(simulations_sym, rm_shape == 0 | is.na(rm_shape))

# Regressions by approach ----

regressions_loglin_sym <- regression(simulations_sym, rm_var = rm_loglin)
regressions_shape_sym <- regression(simulations_sym, rm_var = rm_shape)

# As the observed mortality rate (mx) is calculated within the life table function by sex and age,
# we first obtain the required inputs: the number of deaths (pop_deaths) and the population at risk (pop_census).

observed_pops_loglin_sym <- simulations_sym %>% 
  group_by(sex, age) %>% 
  mutate(pop_deaths = sum(deaths_loglin),
         pop_census = sum(pop_at_risk)) %>% 
  select(sex, age, pop_deaths, pop_census) %>% 
  distinct()

observed_pops_shape_sym <- simulations_sym %>% 
  group_by(sex, age) %>% 
  mutate(pop_deaths = sum(deaths_shape),
         pop_census = sum(pop_at_risk)) %>% 
  select(sex, age, pop_deaths, pop_census) %>% 
  distinct()

# This observed mortality rate is then adjusted by multiplying it by the R factor, which is computed as a function of 
# the parameters a, b, alpha, and beta.
# The adjustment is applied for each combination of quintile, age, and sex in the cohort and lagged cohort methods,
# and by quintile and sex in the period method. The R factor is derived using Equation 3 from Hendi et al. (2021).
# The computation of the factor was implemented using both forms of Equation 3: the definite integral form and the closed-form 
# expression. This allowed us to confirm that both approaches yield equivalent results (see Regression_based_fxs.R).

# We apply the function to compute the R factor, iterating over each quintile, and store the results in a separate 
# dataset for each of the methods.

R_results_loglin_sym <- compute_R_results(regressions_loglin_sym, quintiles_df)
R_results_shape_sym  <- compute_R_results(regressions_shape_sym, quintiles_df)

# Next, we join datasets to get the final input 

cohort_mortality_by_quintile_loglin_sym <- R_results_loglin_sym %>% 
  left_join(observed_pops_loglin_sym, by = c("age", "sex")) 

cohort_mortality_by_quintile_shape_sym <- R_results_shape_sym %>% 
  left_join(observed_pops_shape_sym, by = c("age", "sex")) 

# Life tables ----

# Finally, we construct the life tables for each method by applying the obtain_lt function

lt_loglin_sym <- obtain_lt(cohort_mortality_by_quintile_loglin_sym)
lt_shape_sym <- obtain_lt(cohort_mortality_by_quintile_shape_sym)

# Save life tables ----

saveRDS(lt_loglin_sym, file = glue('{cnst$path_out}/life_table_regression_based_sim_loglin_symmetric.rds'))
saveRDS(lt_shape_sym, file = glue('{cnst$path_out}/life_table_regression_based_sim_shape_symmetric.rds'))

rm(list = ls())
