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

# Load census data ---- 

pop_deaths_2017 <- readRDS(glue('{wd}/Out/pop_deaths_2017.rds'))

distribution <- pop_deaths_2017 %>% 
  group_by(sex, age) %>% 
  summarise(pop = sum(n_census),
            total_deaths = sum(n_deaths))

# Life expectancy function

# simple piecewise-exponential life-table

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

# This is not by quintile

global_le_central_estimates <-
  distribution %>%
  arrange(sex, age) %>%
  group_by(sex) %>%
  group_modify(~{
    CalculateLifeTable(df =.x,x = age,Dx = total_deaths,Ex = pop)
  }) %>%
  ungroup()

# Create Poisson life-table replicates 

# create life table replicates by group and sex
# based on repeatedly sampling death counts from a Poisson

simulations <-
  distribution %>%
  expand_grid(id_sim = 1:500) %>%
  group_by(sex, age) %>%
  mutate(death_total_sim = rpois(500, total_deaths)) %>%
  arrange(age,sex) %>%
  group_by(id_sim, sex) %>%
  group_modify(~{
    CalculateLifeTable(df =.x,x = age,Dx = death_total_sim,Ex = pop)
  }) %>%
  ungroup()

# Assemble table with ex statistics 

# central estimates of life-expectancy, annual life-expectancy difference,
# and average annual life-expectancy difference 2015 to 2020

# 95% uncertainty intervals around the central estimates

confidence_intervals <-
  simulations %>%
  select(id_sim, sex, x, mx, ex, sd0) %>%
  arrange(id_sim, sex, x) %>%
  group_by(sex, x) %>%
  summarise(
    ex_q025 = quantile(ex, 0.025, na.rm = TRUE),
    ex_q975 = quantile(ex, 0.975, na.rm = TRUE),
    sd0_q025 = quantile(sd0, 0.025, na.rm = TRUE),
    sd0_q975 = quantile(sd0, 0.975, na.rm = TRUE)
  )

# assemble all the ex statistics in a single table
# for further computation

global_life_table <-
  left_join(
    global_le_central_estimates,
    confidence_intervals
  )

# Save the regrouped life table input data

saveRDS(global_life_table, file = glue('{cnst$path_out}/global_life_table.rds'))

rm(list = ls())
