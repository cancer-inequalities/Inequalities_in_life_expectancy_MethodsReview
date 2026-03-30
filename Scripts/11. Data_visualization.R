# Data visualization ----

# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(here, glue, dplyr, tidyr, purrr, ggplot2)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_figures = glue('{wd}/Figures')
  path_tmp = glue('{wd}/tmp')
})

# Mortality ratio vs midpoint cumulative distribution

pop_deaths_2017 <- readRDS(glue('{cnst$path_out}/pop_deaths_all_years_ed_2017.rds'))

# We need to create a quinquennial age distribution variable and "regroup" 
# the years of education 

pop_deaths_2017 <- pop_deaths_2017 %>%
  mutate(
    sex = factor(as.numeric(sex), labels = c("Men", "Women")),
    quin_age = cut(age,
             breaks = seq(26, 101, by = 5),
             labels = c("26-30", "31-35", "36-40", "41-45", 
                        "46-50", "51-55", "56-60", "61-65", "66-70", 
                        "71-75", "76-80", "81-85", "86-90", "91-95", "96-100"),
             right = FALSE),
    ed_cat = factor(case_when(
      years_ed == 0 ~ "0",
      years_ed >= 1 & years_ed < 7 ~ "1 - 6",
      years_ed >= 7 & years_ed < 13 ~ "7 - 12",
      years_ed >= 13 & years_ed < 17 ~ "13 - 16",
      years_ed >= 17 ~ "17 - 20"),
      levels = c("0", "1 - 6", "7 - 12", "13 - 16", "17 - 20")))

# Generate a dataset grouped by sex, age group, and educational category

data2 <- pop_deaths_2017 %>%
  group_by(sex, quin_age, ed_cat) %>%
  summarise(
    n_census = sum(n_census),
    n_deaths = sum(n_deaths),
    .groups = 'drop') 

# Create a table for analyzing the mortality ratio in relation to the midpoint of the cumulative distribution

table <- data2 %>%
  mutate(n_denom = n_census - n_deaths) %>%
  group_by(sex, quin_age) %>%                      
  arrange(ed_cat, .by_group = TRUE) %>%      
  mutate(p = n_denom / sum(n_denom),             
         cd = cumsum(p),                     
         cd_prev = lag(cd, default = 0),     
         pc = (cd_prev + cd)/2,              
         tm = 100000*n_deaths/n_census) %>% 
  ungroup()    

table <- table %>%
  group_by(sex, quin_age) %>% 
  mutate(
    n_deaths_total = sum(n_deaths),     
    n_census_total = sum(n_census), 
    ta = n_deaths_total / n_census_total*100000 
  ) %>%
  ungroup()

table <- table %>%
  mutate(rm = tm/ta)

# Figures ----

## Mortality ratio vs Midpoint cumulative distribution by age group ----

women_table <- table %>%
  filter(sex == "Women")

men_table <- table %>%
  filter(sex == "Men")

# Women

women_plot <- ggplot(women_table, aes(x = pc, y = rm, color = quin_age, group = quin_age)) +
  geom_point(size = 2) +
  labs(
    x = expression(
      "Midpoint of cumulative distribution " * 
        ((C[a[i]*", "*s[j]*", "*y[k]]^{e[l]} + C[a[i]*", "*s[j]*", "*y[k]]^{e[l-1]}) / 2)
    ),
    y = expression("Mortality ratio"~(R[a[i]*", "*s[j]*", "*y[k]*", "*e[l]])),
    color = "Age group"
  ) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  theme_minimal()

ggsave(filename = as.character(glue('{cnst$path_figures}/mortality_ratio_vs_midpoint_cumulative_distribution_women.png')), 
       plot = women_plot, 
       width = 12, 
       height = 8,
       dpi = 300)

# Men

men_plot <- ggplot(men_table, aes(x = pc, y = rm, color = quin_age, group = quin_age)) +
  geom_point(size = 2) +
  labs(
    x = expression(
      "Midpoint of cumulative distribution " * 
        ((C[a[i]*", "*s[j]*", "*y[k]]^{e[l]} + C[a[i]*", "*s[j]*", "*y[k]]^{e[l-1]}) / 2)
    ),
    y = expression("Mortality ratio"~(R[a[i]*", "*s[j]*", "*y[k]*", "*e[l]])),
    color = "Age group"
  ) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  theme_minimal()

ggsave(filename = as.character(glue('{cnst$path_figures}/mortality_ratio_vs_midpoint_cumulative_distribution_men.png')), 
       plot = men_plot, 
       width = 12, 
       height = 8,
       dpi = 300)

# Observed mortality rate by age vs predicted mortality rate according to the regression based approach ----

cohort_proportions <- readRDS(glue('{cnst$path_out}/cohort_proportions_for_visualization.rds'))

lagged_cohort_proportions <- readRDS(glue('{cnst$path_out}/lagged_cohort_proportions_for_visualization.rds'))

proportions <- cohort_proportions %>% 
  select(sex, age, years_ed, n_census, n_deaths, pop_census, pop_deaths, pc, rm) %>% 
  rename(pc_coh = pc) %>% 
  left_join(lagged_cohort_proportions %>% select(sex, age, years_ed, pc) %>% rename(pc_lagco = pc), by = c("sex", "age", "years_ed"))

ages <- c(26, 40, 60, 80)

# Generate plots of mortality ratio versus cumulative proportion ----

# Methods vector

methods <- c("Cohort", "Lagged-cohort")

# Function to convert sex codes (1/2) to text labels ("male"/"female")

sex_labels <- function(value) {
  ifelse(value == 1, "Men", "Women")
}

for (method in methods) {
  for (a in ages) {
    
    var_percentile <- if (method == "Cohort") "pc_coh" else "pc_lagco"
    
    df_filtered <- proportions %>%
      filter(age == a) %>%
      group_by(sex, age) %>%
      mutate(
        alpha = coef(lm(log(rm) ~ .data[[var_percentile]]))[1],
        beta = coef(lm(log(rm) ~ .data[[var_percentile]]))[2],
        log_rm_predicted = alpha + beta * .data[[var_percentile]],
        rm_predicted = exp(log_rm_predicted),
        log_rm = log(rm)
      ) %>%
      ungroup() %>%
      mutate(
        rm_numerator = n_deaths / n_census,
        rm_denominator = pop_deaths / pop_census,
        rm_predicted_numerator = rm_predicted * rm_denominator,
        sex_label = sex_labels(sex)
      )
    
    graph <- ggplot(df_filtered, aes(x = .data[[var_percentile]])) +
      geom_line(aes(y = rm_predicted_numerator, group = "Predicted", colour = "Predicted"), linewidth = 1, linetype = "dashed") +
      geom_line(aes(y = rm_numerator, group = "Observed", colour = "Observed"), linewidth = 1, linetype = "solid") +
      facet_wrap(~ sex_label, ncol = 1) +
      scale_color_manual(values = c("Observed" = "#2C3E50", "Predicted" = "#95A5A6"), name = NULL) +
      labs(
        title = paste0(method, " approach - Age ", a),
        x = bquote("Cumulative proportion of education (" ~ c[ a[i] * "," * s[j] * "," * y[k] ] ^ { e[l] } ~ ")"),
        y = bquote("Mortality rate (" ~ m[ a[i] * "," * s[j] * "," * y[k] * "," * e[l] ] ~ ")")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.title.x = element_text(margin = margin(t = 15)),  
            axis.title.y = element_text(margin = margin(r = 15)),  
            plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
            plot.margin = margin(t = 20, r = 26, b = 10, l = 6))
    
    file_name <- paste0("mortality_rate_cumulative_proportion_", method, "_age_", a, ".png")
    file_path  <- file.path(cnst$path_figures, file_name)
    ggsave(filename = file_path, plot = graph, width = 6, height = 7, dpi = 300)
  }
}

# Mortality rates by age and sex for each quintile for each combination of method and approach ----

# Load life tables [compositional adjustment]

lagged_cohort <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_lagged_cohort_approach.rds')) %>% 
select(quintile, sex, x, mx) %>% 
filter( x %in% ages)

cohort <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_cohort_approach.rds')) %>% 
  select(quintile, sex, x, mx) %>% 
  filter( x %in% ages)

period <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_period_approach.rds')) %>%
  select(quintile, sex, x, mx) %>%
  filter(x %in% ages)

# Load life tables [regression based]

lagged_cohort_rb <- readRDS(glue('{cnst$path_out}/life_table_regression_based_lagged_cohort.rds')) %>%
  select(quintile, sex, x, mx) %>%
  filter(x %in% ages)

cohort_rb <- readRDS(glue('{cnst$path_out}/life_table_regression_based_cohort.rds'))  %>%
  select(quintile, sex, x, mx) %>%
  filter(x %in% ages)

period_rb <- readRDS(glue('{cnst$path_out}/life_table_regression_based_period.rds'))  %>%
  select(quintile, sex, x, mx) %>%
  filter(x %in% ages)

# Integrate datasets according to method using a left join ----

lagged_cohort_join <- lagged_cohort %>% 
  left_join(lagged_cohort_rb, by = c("quintile", "sex", "x"))

cohort_join <- cohort %>% 
  left_join(cohort_rb, by = c("quintile", "sex", "x"))

period_join <- period %>% 
  left_join(period_rb, by = c("quintile","sex","x"))

# Figures by approach ----

# Lagged cohort

lagged_cohort_join <- lagged_cohort_join %>%
  mutate(sex = factor(sex, levels = c(1, 2), labels = c("Men", "Women")),
         quintile = factor(quintile, levels = c("quintile1", "quintile2", "quintile3","quintile4",
                                              "quintile5"), labels = c("Quintile 1", "Quintile 2",
                                                                      "Quintile 3", "Quintile 4",
                                                                      "Quintile 5")))

lagged_cohort_rates <- ggplot(lagged_cohort_join, aes(x = x)) +
  geom_line(aes(y = mx.y, group = "Regression-based", colour = "Regression-based"), linewidth = 1, , linetype = "dashed") +
  geom_line(aes(y = mx.x, group = "Compositional adjustment", colour = "Compositional adjustment"), linewidth = 1, linetype = "solid") +
  facet_grid(sex ~ quintile) + 
  scale_color_manual(values = c("Compositional adjustment" = "blue", "Regression-based" = "red"), , name = NULL) +
  labs(
    x = "Age",
    y = "Mortality rate (Log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),  
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.background = element_rect(fill = "gray58"))  


ggsave(filename = as.character(glue('{cnst$path_figures}/lagged_cohort_mortality_rate_sex_age_quintile.png')), 
       plot = lagged_cohort_rates, 
       width = 8, 
       height = 6,
       dpi = 300)

# Cohort

cohort_join <- cohort_join %>%
  mutate(sex = factor(sex, levels = c(1, 2), labels = c("Men", "Women")),
         quintile = factor(quintile, levels = c("quintile1", "quintile2", "quintile3","quintile4",
                                              "quintile5"), labels = c("Quintile 1", "Quintile 2",
                                                                      "Quintile 3", "Quintile 4",
                                                                      "Quintile 5")))

cohort_rates <- ggplot(cohort_join, aes(x = x)) +
  geom_line(aes(y = mx.y, group = "Regression-based", colour = "Regression-based"), linewidth = 1, , linetype = "dashed") +
  geom_line(aes(y = mx.x, group = "Compositional adjustment", colour = "Compositional adjustment"), linewidth = 1, linetype = "solid") +
  facet_grid(sex ~ quintile) + 
  scale_color_manual(values = c("Compositional adjustment" = "blue", "Regression-based" = "red"), , name = NULL) +
  labs(
    x = "Age",
    y = "Mortality rate (Log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),  
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.background = element_rect(fill = "gray58")) 

ggsave(filename = as.character(glue('{cnst$path_figures}/cohort_mortality_rate_sex_age_quintile.png')), 
       plot = lagged_cohort_rates, 
       width = 8, 
       height = 6,
       dpi = 300)

# Period

period_join <- period_join %>%
  mutate(sex = factor(sex, levels = c(1, 2), labels = c("Men", "Women")),
         quintile = factor(quintile, levels = c("quintile1", "quintile2", "quintile3","quintile4",
                                              "quintile5"), labels = c("Quintile 1", "Quintile 2",
                                                                      "Quintile 3", "Quintile 4",
                                                                      "Quintile 5")))

period_rates <- ggplot(period_join, aes(x = x)) +
  geom_line(aes(y = mx.y, group = "Regression-based", colour = "Regression-based"), linewidth = 1, , linetype = "dashed") +
  geom_line(aes(y = mx.x, group = "Compositional adjustment", colour = "Compositional adjustment"), linewidth = 1, linetype = "solid") +
  facet_grid(sex ~ quintile) + 
  scale_color_manual(values = c("Compositional adjustment" = "blue", "Regression-based" = "red"), , name = NULL) +
  labs(
    x = "Age",
    y = "Mortality rate (Log scale)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),  
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.background = element_rect(fill = "gray58"))  


ggsave(filename = as.character(glue('{cnst$path_figures}/period_mortality_rate_sex_age_quintile.png')), 
       plot = lagged_cohort_rates, 
       width = 8, 
       height = 6,
       dpi = 300)

# Simulations figures ----

# Asymmetric figure

compositional_loglin_asy <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_sim_loglin_asymmetric.rds'))
compositional_shape_asy  <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_sim_shape_asymmetric.rds'))

regression_loglin_asy <- readRDS(glue('{cnst$path_out}/life_table_regression_based_sim_loglin_asymmetric.rds'))
regression_shape_asy  <- readRDS(glue('{cnst$path_out}/life_table_regression_based_sim_shape_asymmetric.rds'))

bounds_loglin_asy <- readRDS(glue('{cnst$path_out}/life_table_bounded_mortality_sim_loglin_asymmetric.rds'))[["life_expectancy_summary"]]
bounds_shape_asy  <- readRDS(glue('{cnst$path_out}/life_table_bounded_mortality_sim_shape_asymmetric.rds'))[["life_expectancy_summary"]]

# Prepare data

compositional_loglin_asy <- compositional_loglin_asy %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(method = "Compositional adjustment", type = "loglin")

compositional_shape_asy <- compositional_shape_asy %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(method = "Compositional adjustment", type = "shape")

regression_loglin_asy <- regression_loglin_asy %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(method = "Regression-based", type = "loglin")

regression_shape_asy <- regression_shape_asy %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(method = "Regression-based", type = "shape")

bounds_loglin_asy <- bounds_loglin_asy %>%
  filter(x == 26) %>%
  mutate(
    method = "Bounded mortality",
    type = "loglin",
    quintile = recode(
      quintile,
      "1" = "quintile1", "2" = "quintile2", "3" = "quintile3",
      "4" = "quintile4", "5" = "quintile5"
    )
  ) %>%
  rename(ex_q025 = low_ex, ex_q975 = upper_ex) %>%
  select(-x) %>%
  mutate(ex = (ex_q025 + ex_q975) / 2)

bounds_shape_asy <- bounds_shape_asy %>%
  filter(x == 26) %>%
  mutate(
    method = "Bounded mortality",
    type = "shape",
    quintile = recode(
      quintile,
      "1" = "quintile1", "2" = "quintile2", "3" = "quintile3",
      "4" = "quintile4", "5" = "quintile5"
    )
  ) %>%
  rename(ex_q025 = low_ex, ex_q975 = upper_ex) %>%
  select(-x) %>%
  mutate(ex = (ex_q025 + ex_q975) / 2)

# Bind data

binded_data_asy <- bind_rows(
  compositional_loglin_asy,
  compositional_shape_asy,
  regression_loglin_asy,
  regression_shape_asy,
  bounds_loglin_asy,
  bounds_shape_asy
)

# Prepare plot data

plot_data_asy <- binded_data_asy %>%
  filter(sex == 2) %>%
  mutate(
    quintile = factor(
      quintile,
      levels = paste0("quintile", 1:5),
      labels = 1:5
    )
  )

# Plot (FINAL CLEAN VERSION)

pos <- position_dodge(width = 0.5)

asymmetric_figure <- ggplot(plot_data_asy, aes(x = method, color = quintile, group = quintile)) +
  
  # Intervals
  geom_linerange(
    aes(
      ymin = ifelse(ex_q025 != ex_q975, ex_q025, NA),
      ymax = ifelse(ex_q025 != ex_q975, ex_q975, NA)
    ),
    position = pos,
    linewidth = 1.0
  ) +
  
  # Points (non-bounded)
  geom_point(
    aes(
      y = ifelse(method != "Bounded mortality", ex, NA)
    ),
    position = pos,
    size = 3.2
  ) +
  
  # Points (bounded, degenerate intervals)
  geom_point(
    aes(
      y = ifelse(method == "Bounded mortality" & ex_q025 == ex_q975, ex_q025, NA)
    ),
    position = pos,
    size = 3.0
  ) +
  
  coord_flip() +
  
  facet_wrap(
    ~type,
    labeller = as_labeller(c(
      loglin = "Log-linear",
      shape  = "Hump-shaped"
    ))
  ) +
  
  scale_x_discrete(
    limits = c(
      "Bounded mortality",
      "Regression-based",
      "Compositional adjustment"
    )
  ) +
  
  scale_color_manual(values = c(
    "#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4"
  )) +
  
  labs(
    x = NULL,
    y = expression(e[26]),
    color = "Quintile"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    
    strip.text = element_text(face = "bold", size = 13),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename = as.character(glue('{cnst$path_figures}/asymmetric_simulations_figure.png')), 
       plot = asymmetric_figure, 
       width = 7.5, 
       height = 5,
       dpi = 300)

rm(list = setdiff(ls(), c("wd", "cnst")))

# Symmetric figure

# Load data

compositional_loglin_sym <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_sim_loglin_symmetric.rds'))
compositional_shape_sym  <- readRDS(glue('{cnst$path_out}/life_table_compositional_adjustment_sim_shape_symmetric.rds'))

regression_loglin_sym <- readRDS(glue('{cnst$path_out}/life_table_regression_based_sim_loglin_symmetric.rds'))
regression_shape_sym  <- readRDS(glue('{cnst$path_out}/life_table_regression_based_sim_shape_symmetric.rds'))

bounds_loglin_sym <- readRDS(glue('{cnst$path_out}/life_table_bounded_mortality_sim_loglin_symmetric.rds'))[["life_expectancy_summary"]]
bounds_shape_sym  <- readRDS(glue('{cnst$path_out}/life_table_bounded_mortality_sim_shape_symmetric.rds'))[["life_expectancy_summary"]]

# Prepare data

compositional_loglin_sym <- compositional_loglin_sym %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(
    method = "Compositional adjustment",
    type = "loglin"
  )

compositional_shape_sym <- compositional_shape_sym %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(
    method = "Compositional adjustment",
    type = "shape"
  )

regression_loglin_sym <- regression_loglin_sym %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(
    method = "Regression-based",
    type = "loglin"
  )

regression_shape_sym <- regression_shape_sym %>%
  filter(x == 26) %>%
  select(sex, quintile, ex_q025, ex, ex_q975) %>%
  mutate(
    method = "Regression-based",
    type = "shape"
  )

bounds_loglin_sym <- bounds_loglin_sym %>%
  filter(x == 26) %>%
  mutate(
    method = "Bounded mortality",
    type = "loglin",
    quintile = recode(
      quintile,
      "1" = "quintile1",
      "2" = "quintile2",
      "3" = "quintile3",
      "4" = "quintile4",
      "5" = "quintile5",
      .default = NA_character_
    )
  ) %>%
  rename(
    ex_q025 = low_ex,
    ex_q975 = upper_ex
  ) %>%
  select(-x) %>%
  mutate(ex = (ex_q025 + ex_q975) / 2)

bounds_shape_sym <- bounds_shape_sym %>%
  filter(x == 26) %>%
  mutate(
    method = "Bounded mortality",
    type = "shape",
    quintile = recode(
      quintile,
      "1" = "quintile1",
      "2" = "quintile2",
      "3" = "quintile3",
      "4" = "quintile4",
      "5" = "quintile5",
      .default = NA_character_
    )
  ) %>%
  rename(
    ex_q025 = low_ex,
    ex_q975 = upper_ex
  ) %>%
  select(-x) %>%
  mutate(ex = (ex_q025 + ex_q975) / 2)

# Bind data

binded_data_sym <- bind_rows(
  compositional_loglin_sym,
  compositional_shape_sym,
  regression_loglin_sym,
  regression_shape_sym, 
  bounds_loglin_sym,
  bounds_shape_sym
)

# Prepare plot data

plot_data_sym <- binded_data_sym %>%
  filter(sex == 2) %>%
  mutate(
    quintile = factor(
      quintile,
      levels = paste0("quintile", 1:5),
      labels = 1:5
    )
  )

data_points_sym <- plot_data_sym %>%
  filter(method != "Bounded mortality")

data_bounds_points_sym <- plot_data_sym %>%
  filter(
    method == "Bounded mortality",
    ex_q025 == ex_q975
  )

# Plot

pos <- position_dodge(width = 0.5)

symmetric_figure <- ggplot(plot_data_sym, aes(x = method, color = quintile, group = quintile)) +
  
  # Intervals
  geom_linerange(
    aes(
      ymin = ifelse(ex_q025 != ex_q975, ex_q025, NA),
      ymax = ifelse(ex_q025 != ex_q975, ex_q975, NA)
    ),
    position = pos,
    linewidth = 1.0
  ) +
  
  # Points (non-bounded)
  geom_point(
    aes(
      y = ifelse(method != "Bounded mortality", ex, NA)
    ),
    position = pos,
    size = 3.2
  ) +
  
  # Points (bounded, degenerate)
  geom_point(
    aes(
      y = ifelse(method == "Bounded mortality" & ex_q025 == ex_q975, ex_q025, NA)
    ),
    position = pos,
    size = 3.0
  ) +
  
  coord_flip() +
  
  facet_wrap(
    ~type,
    labeller = as_labeller(c(
      loglin = "Log-linear",
      shape  = "Hump-shaped"
    ))
  ) +
  
  scale_x_discrete(
    limits = c(
      "Bounded mortality",
      "Regression-based",
      "Compositional adjustment"
    )
  ) +
  
  scale_color_manual(values = c(
    "#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4"
  )) +
  
  labs(
    x = NULL,
    y = expression(e[26]),
    color = "Quintile"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.3, "cm"),
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    
    strip.text = element_text(face = "bold", size = 13),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename = as.character(glue('{cnst$path_figures}/symmetric_simulations_figure.png')), 
       plot = symmetric_figure, 
       width = 7.5, 
       height = 5,
       dpi = 300)

rm(list = ls())