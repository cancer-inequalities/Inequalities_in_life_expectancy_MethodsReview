# Life tables and life expectancy estimations by bounded mortality method using simulated mortality rates: ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(dplyr,readxl, writexl, haven, rlang, glue, here)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Functions ----

# Life table function 

CalculateLifeTable <- function(df, x, nx = c(rep(1, 74), Inf), mx) {
  
  df %>%
    transmute(
      x  = {{ x }},
      nx = nx,
      mx = {{ mx }},
      px = exp(-mx * nx),
      qx = 1 - px,
      lx = head(cumprod(c(1, px)), -1),
      dx = c(-diff(lx), tail(lx, 1)),
      Lx = if_else(mx <= 0, lx * nx, dx / mx),
      Tx = rev(cumsum(rev(Lx))),
      ex = Tx / lx,
      sd0 = c(sqrt(sum(dx * (x + 0.5 - ex[1])^2)), rep(0, length(x) - 1))
    )
}

# Build a life table using a selected mortality rate variable 

build_life_table <- function(data, rate_var) {
  
  rate_sym <- sym(rate_var)
  
  data %>%
    arrange(quintile, sex, age) %>%
    group_by(quintile, sex) %>%
    group_modify(~ CalculateLifeTable(.x, x = age, mx = !!rate_sym)) %>%
    ungroup()
}

# Process a single dataset with lower and upper mortality bounds 

process_life_tables_bounds <- function(data, dataset_name, out_dir = cnst$path_out) {
  
  # Build life tables for both mortality rate bounds
  
  life_table_lb <- build_life_table(data, "rate_lb")
  life_table_ub <- build_life_table(data, "rate_ub")
  
  # Consolidated life expectancy dataset
  # rate_ub -> higher mortality -> lower life expectancy
  # rate_lb -> lower mortality  -> upper life expectancy
  
  life_expectancy_summary <- life_table_ub %>%
    select(sex, x, quintile, ex) %>%
    rename(low_ex = ex) %>%
    left_join(
      life_table_lb %>%
        select(sex, x, quintile, ex) %>%
        rename(upper_ex = ex),
      by = c("sex", "x", "quintile")
    )
  
  # Output file names
  
  file_lb <- glue("{out_dir}/{dataset_name}_life_table_lower_bounds.rds")
  file_ub <- glue("{out_dir}/{dataset_name}_life_table_upper_bounds.rds")
  file_ex <- glue("{out_dir}/{dataset_name}_life_expectancy_lower_upper.rds")
  
  # Save outputs
  
  saveRDS(life_table_lb, file = file_lb)
  saveRDS(life_table_ub, file = file_ub)
  saveRDS(life_expectancy_summary, file = file_ex)
  
  # Return outputs in memory
  
  list(
    dataset_name = dataset_name,
    life_table_lb = life_table_lb,
    life_table_ub = life_table_ub,
    life_expectancy_summary = life_expectancy_summary,
    files = c(
      lower_bounds = file_lb,
      upper_bounds = file_ub,
      life_expectancy = file_ex
    )
  )
}

# A) Asymmetric simulations ----

# Load data ---- 

bounds_loglin_asy <- read_dta(glue('{wd}/Out/mortality_bounds_by_age_sex_quintile_cohort_loglin.dta'))
bounds_shape_asy <- read_dta(glue('{wd}/Out/mortality_bounds_by_age_sex_quintile_cohort_simshape.dta'))

# Compute mortality rates as inputs for the life table (mx) ----
# Expected deaths are expressed per 100,000 population and must be rescaled accordingly.

bounds_loglin_asy <- bounds_loglin_asy %>% 
  mutate(rate_lb = mort_lb/100000, # lower bound rate
         rate_ub = mort_ub/100000) # upper bound rate

bounds_shape_asy <- bounds_shape_asy %>% 
  mutate(rate_lb = mort_lb/100000, # lower bound rate
         rate_ub = mort_ub/100000) # upper bound rate

# Run workflow ----

results_life_tables_loglin_asy <- process_life_tables_bounds(
  data = bounds_loglin_asy,
  dataset_name = "bounds_loglin_asy",
  out_dir = cnst$path_out
)

results_life_tables_shape_asy <- process_life_tables_bounds(
  data = bounds_shape_asy,
  dataset_name = "bounds_shape_asy",
  out_dir = cnst$path_out
)

# Save lists ----

saveRDS(results_life_tables_loglin_asy, file = glue('{cnst$path_out}/life_table_bounded_mortality_sim_loglin_asymmetric.rds'))
saveRDS(results_life_tables_shape_asy, file = glue('{cnst$path_out}/life_table_bounded_mortality_sim_shape_asymmetric.rds'))

rm(list = setdiff(ls(), c("wd", "cnst", "CalculateLifeTable", "build_life_table", "process_life_tables_bounds")))

# B). Symmetric simulations ----

# Load data ---- 

bounds_loglin_sym <- read_dta(glue('{wd}/Out/mortality_bounds_by_age_sex_quintile_cohort_loglin_simmetric.dta'))
bounds_shape_sym <- read_dta(glue('{wd}/Out/mortality_bounds_by_age_sex_quintile_cohort_simshape_simulated_mortality_data_simmetric.dta'))

# Compute mortality rates as inputs for the life table (mx) ----
# Expected deaths are expressed per 100,000 population and must be rescaled accordingly.

bounds_loglin_sym <- bounds_loglin_sym %>% 
  mutate(rate_lb = mort_lb/100000, # lower bound rate
         rate_ub = mort_ub/100000) # upper bound rate

bounds_shape_sym <- bounds_shape_sym %>% 
  mutate(rate_lb = mort_lb/100000, # lower bound rate
         rate_ub = mort_ub/100000) # upper bound rate

# Run workflow ----

results_life_tables_loglin_sym <- process_life_tables_bounds(
  data = bounds_loglin_sym,
  dataset_name = "bounds_loglin_sym",
  out_dir = cnst$path_out
)

results_life_tables_shape_sym <- process_life_tables_bounds(
  data = bounds_shape_sym,
  dataset_name = "bounds_shape_sym",
  out_dir = cnst$path_out
)

# Save lists ----

saveRDS(results_life_tables_loglin_sym, file = glue('{cnst$path_out}/life_table_bounded_mortality_sim_loglin_symmetric.rds'))
saveRDS(results_life_tables_shape_sym, file = glue('{cnst$path_out}/life_table_bounded_mortality_sim_shape_symmetric.rds'))

rm(list = ls())
