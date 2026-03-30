# Generation of simulated mortality inputs under log-linear and hump-shape models: ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(haven, dplyr, glue, here, ggplot2)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Load data ----

data <- read_dta(glue('{wd}/Out/mortality_rates_2017_ranked_lagged_cohort.dta'))

# Ridit estimation ----

data <- data %>%
 mutate(ridit = rank / 100)

# Parameter ----

beta <- -3

sim_data <- data %>%
 group_by(age, sex) %>%
 mutate(
  
# Mean mortality ----

mean_rate = weighted.mean(rate, n_census),
  
  # Log-linear model (deterministic)

  log_z_loglin = beta * ridit,
  
  z_loglin = exp(log_z_loglin),
  
  rate_loglin = mean_rate * z_loglin /
   weighted.mean(z_loglin, n_census),
  
  # Extrema (on rate scale)

  m_max = max(rate_loglin),
  m_min = min(rate_loglin),
  
  # Hump on log scale
  
  # Locate closest observed point to ridit = 0.2

  idx_peak = which.min(abs(ridit - 0.35)),
  
  # Distance from peak (in rank order)

  dist = abs(row_number() - idx_peak),
  
  # Symmetric hump centered at peak

  shape_unit = 1 - dist / max(dist),
  
  # Smooth curvature

  shape_unit = shape_unit^1.5,
  
  # Normalize log scale between min and max

  log_min = min(log_z_loglin),
  log_max = max(log_z_loglin),
  
  # Construct log-risk hump

  log_z_hump = log_min + (log_max - log_min) * shape_unit,
  
  z_hump = exp(log_z_hump),
  
  # Scale to match mean exactly
  
  rate_hump_raw = mean_rate * z_hump /
   weighted.mean(z_hump, n_census),
  
  rate_sim_shape = rate_hump_raw
  
 ) %>% 
 ungroup()

# Plot ----

ages_to_plot <- c(30, 40, 50, 60)

plot_data <- sim_data %>%
 filter(age %in% ages_to_plot, sex == 1)

ggplot(plot_data, aes(x = ridit)) +
 
 geom_point(aes(y = rate), size = 2) +
 
 geom_line(aes(y = rate_loglin, color = "Log-linear model"), linewidth = 1) +
 
 geom_line(aes(y = rate_sim_shape, color = "Hump-shape model"), linewidth = 1) +
 
 facet_wrap(~age, scales = "free_y") +
 
 labs(
  x = "Ridit (education rank)",
  y = "Mortality rate",
  color = "Model"
 ) +
 
 theme_minimal()

# Check means ----

check_means <- sim_data %>%
 group_by(age, sex) %>%
 summarise(
  mean_loglin = weighted.mean(rate_loglin, n_census),
  mean_hump   = weighted.mean(rate_sim_shape, n_census),
  diff = mean_hump - mean_loglin,
  .groups = "drop"
 )

print(check_means)

sim_data$rate_sim_loglin<-sim_data$rate_loglin

# Save data ----

write_dta(
 sim_data,
 path = glue('{cnst$path_out}/simulated_mortality_data.dta')
)

rm(list = ls())
