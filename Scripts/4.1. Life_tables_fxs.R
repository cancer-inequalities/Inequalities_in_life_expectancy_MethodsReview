# Create the input datasets to calculate life tables by adjustment approach

create_input <- function(data, pop_data){
  data <- data %>% 
    left_join(pop_data, by = c("age", "sex", "years_ed")) %>% 
    mutate(
      denomq1 = round(n_census * wq1, 0),
      denomq2 = round(n_census * wq2, 0),
      denomq3 = round(n_census * wq3, 0),
      denomq4 = round(n_census * wq4, 0),
      denomq5 = round(n_census * wq5, 0),
      numq1 = round(n_deaths * wq1, 0),
      numq2 = round(n_deaths * wq2, 0),
      numq3 = round(n_deaths * wq3, 0),
      numq4 = round(n_deaths * wq4, 0),
      numq5 = round(n_deaths * wq5, 0)
    )
  
  counts <- data %>% 
    group_by(age, sex) %>% 
    summarise(
      quintile1 = sum(denomq1, na.rm = TRUE),
      quintile2 = sum(denomq2, na.rm = TRUE),
      quintile3 = sum(denomq3, na.rm = TRUE),
      quintile4 = sum(denomq4, na.rm = TRUE),
      quintile5 = sum(denomq5, na.rm = TRUE)
    ) %>% 
    pivot_longer(cols = starts_with("quintile"), 
                 names_to = "quintile", 
                 values_to = "counts")
  
  deaths <- data %>% 
    group_by(age, sex) %>%
    summarise(
      quintile1 = sum(numq1, na.rm = TRUE),
      quintile2 = sum(numq2, na.rm = TRUE),
      quintile3 = sum(numq3, na.rm = TRUE),
      quintile4 = sum(numq4, na.rm = TRUE),
      quintile5 = sum(numq5, na.rm = TRUE)
    ) %>% 
    pivot_longer(cols = starts_with("quintile"), 
                 names_to = "quintile", 
                 values_to = "total_deaths")
  
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
      CalculateLifeTable(df =.x,x = age,Dx = death_total_sim,Ex = counts)
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










