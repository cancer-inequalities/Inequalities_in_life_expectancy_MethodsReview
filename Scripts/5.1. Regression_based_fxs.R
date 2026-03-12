# regression function

regression <- function(data){
  
  # Create a data frame  to store the results
  
  results <- data.frame(sex = integer(), age = integer(), a = numeric(), b = numeric())
  
  # Group by sex and age and apply a regression 
  
  regressions <- data %>%
    group_by(sex, age) %>%
    do({
      modelo <- lm(log(rm) ~ pc, data = .)  # Adjust the regression to each group 
      data.frame(alpha = coef(modelo)[1], beta = coef(modelo)[2])  # Recover the coefficient alpha y beta
    })
  
  return(regressions)
  
}

# Fórmula con la integral

R_factor <- function(alpha, beta, a, b) {
  R_factor <- function(x) {
    exp(alpha + beta * x)
  }
  result <- (1 / (b - a)) * integrate(R_factor, lower = a, upper = b)$value
  return(result)
}

# última fórmula de la igualdad, just to check the results are the same

second_formula <- function(alpha, beta, a, b) {
  result <- (exp(alpha) / (beta * (b - a))) * (exp(beta * b) - exp(beta * a))
  return(result)
}

# Procesar cada fila y aplicar las fórmulas

processing <- function(row, quintile, a, b) {
  alpha <- as.numeric(row['alpha'])
  beta <- as.numeric(row['beta']) 
  integral_val <- R_factor(alpha, beta, a, b)
  final_val <- second_formula(alpha, beta, a, b)
  return(data.frame(sex = row['sex'], age = row['age'], quintile = quintile,
                    R_factor = integral_val, second_formula = final_val))
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




  


  






