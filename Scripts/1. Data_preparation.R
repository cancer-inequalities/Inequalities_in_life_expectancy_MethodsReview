# Load packages ----

if(!require("pacman")) install.packages("pacman")

pacman::p_load(here, glue, dplyr)

rm(list = ls())

# Constants ----

wd <- here()
cnst <- list()

cnst <- within(cnst, {
  path_out = glue('{wd}/Out')
  path_tmp = glue('{wd}/tmp')
})

# Load data ----

# Chilean Census 2017 Database – Individual Microdata

census_2017 <- readRDS(glue('{wd}/Data_folder/Census_2017.rds'))
                                   
# Chilean Mortality 2017 Database - Individual microdata

deaths_counts_2017 <- readRDS(glue('{wd}/Data_folder/Deaths_counts_2017.rds'))

# Create full combination dataset ----

# Generate all possible combinations of sex, age, and years of education.
# This ensures that no combination is missing in later merges,
# even if some combinations are not present in the input datasets.

sex <- c(1,2)
age <- 26:100
years_ed <- 0:20

full_combinations <- expand.grid(sex = sex,
                                       age = age,
                                       years_ed = years_ed)

# Process data ----

## Census ----

census_2017$highest_grade[census_2017$highest_grade == 98 | census_2017$highest_grade == 99] <- NA

# Create ed_type and years_ed columns

census_2017 <- census_2017 %>%
  mutate(ed_type = case_when(
    highest_level <= 3 ~ 0, # no formal education
    highest_level == 4 | highest_level == 5 ~ 1, # Primary education
    highest_level == 7 | highest_level == 8 ~ 2, # Secondary education
    highest_level == 11 ~ 3, # Technical training center/Professional institute
    highest_level == 12 ~ 4, # University education / Undergraduate degree
    highest_level == 13 ~ 5, # Master's degree
    highest_level == 14 ~ 6, # Doctorate / PhD
    TRUE ~ NA
  ),
  years_ed = case_when(
    currently_attending == 3 ~ 0,
    ed_type == 0 ~ 0,
    ed_type == 1 & highest_grade < 8 ~ highest_grade,
    ed_type == 1 & highest_grade >= 8  ~ 8,
    ed_type == 2 & highest_grade < 4 ~ highest_grade + 8,
    ed_type == 2 & highest_grade >= 4 ~ 12,
    ed_type == 3 & highest_grade < 3 ~ highest_grade + 12,
    ed_type == 3 & highest_grade >= 3 ~ 15,
    ed_type == 4 & highest_grade < 5 ~ highest_grade + 12,
    ed_type == 4 & highest_grade >= 5 ~ 17,
    ed_type == 5 & highest_grade < 2 ~ highest_grade + 17,
    ed_type == 5 & highest_grade >= 2 ~ 19,
    ed_type == 6 & highest_grade == 0 ~ 19,
    ed_type == 6 & highest_grade >= 1 ~ 20,
    # Teacher training / Normal school education
    highest_level == 10 & highest_grade < 6 ~ highest_grade + 6,
    highest_level == 10 & highest_grade >= 6 ~ 12,
    # Preparatory primary education
    highest_level == 6 & highest_grade < 6 ~ highest_grade,
    highest_level == 6 & highest_grade >= 6 ~ 6,
    # Secondary education – humanities track
    highest_level == 9 & highest_grade < 6 ~ highest_grade + 6,
    highest_level == 9 & highest_grade >= 6 ~ 12,
    TRUE ~ NA
  ))

# Group by sex, age and years of education and get the counts

grouped_census <- census_2017 %>%
  filter(!is.na(years_ed)) %>%
  group_by(sex, age, years_ed) %>%
  summarise(n_census = n())

# Join to complete full combinations

grouped_census <- full_combinations %>%   
  left_join(grouped_census, by = c("sex","age","years_ed")) %>% 
  mutate(n_census = ifelse(is.na(n_census), 0, n_census))

## Deaths counts 2017 ----

# Manage sex 

sort(unique(deaths_counts_2017$sex)) # values = 1, 2 and 9

# Drop "9" values

deaths_counts_2017 <- deaths_counts_2017 %>%
  filter(sex != 9) 

# Manage age

# Transform age_type and age_value to get age column (age in years). 

sort(unique(deaths_counts_2017$age_type)) # values: 1-4 and 9

# Drop "9" values

deaths_counts_2017 <- deaths_counts_2017 %>%
  filter(age_type != 9)

# Convert age to years only when age_type == 1 (i.e., age is already in years).
# For all other cases (months, days, or hours), set age to 0.

deaths_counts_2017 <- deaths_counts_2017 %>%
  mutate(age = if_else(
    condition = age_type == 1,
    true = age_value,
    false = 0))

# Explore the new column age

summary(deaths_counts_2017$age)   # min: 0, max = 113
sort(unique(deaths_counts_2017$age)) 

# Transform all ages above 100 years old, then filter ages 

deaths_counts_2017 <- deaths_counts_2017 %>%
  mutate(age = if_else(
    condition = age > 100,
    true = 100,
    false = age)) %>%
  filter(age >= 26)

# We will work with the columns education_level and education_grade
# These are numeric variables, and we will explore their possible values

sort(unique(deaths_counts_2017$education_level)) # values: 1-5, 9
sort(unique(deaths_counts_2017$education_grade)) # values: 0-9

table(deaths_counts_2017$education_level)
table(deaths_counts_2017$education_grade)

# A new variable called "years_ed" needs to be created, based on education_level, education_grade, or a combination of both

deaths_counts_2017 <- deaths_counts_2017 %>%
  mutate(education_level = if_else(
    condition = education_level == 9,
    true = NA,
    false = education_level),
    education_grade = if_else(
      condition = education_grade == 9,
      true = NA,
      false = education_grade),
    years_ed = case_when(
      education_level == 5 ~ 0,
      education_level == 4 ~ education_grade,
      education_level == 3 ~ education_grade + 6,
      education_level == 2 ~ education_grade + 8,
      education_level == 1 ~ education_grade + 12)) %>%
  filter(!is.na(years_ed))

### Group data to compute the number of deaths (n_deaths)

deaths_counts_2017 <- deaths_counts_2017 %>%
  group_by(sex, age, years_ed) %>%
  summarise(n_deaths = n(), .groups = "drop")

deaths_counts_2017 <- full_combinations %>% 
  left_join(deaths_counts_2017, by = c("sex", "age", "years_ed")) %>% 
  mutate(n_deaths = ifelse(is.na(n_deaths), 0, n_deaths))

## Census - deaths counts join ----

pop_deaths_2017 <- grouped_census %>% 
  left_join(deaths_counts_2017, by = c("sex","age","years_ed"))

pop_deaths_2017 <- pop_deaths_2017 %>%
  arrange(sex, age, years_ed) 

## n_deaths and n_census correction functions ----

# Correct n_deaths = 0

adjust_n_deaths <- function(df) {
  for (i in 2:nrow(df)) {
    if (df$n_deaths[i] == 0) {
      # Average with the previous row
            df$n_deaths[i] <- (df$n_deaths[i - 1] + df$n_deaths[i]) / 2
      df$n_deaths[i - 1] <- df$n_deaths[i]  # Replace the previous row with the computed average
    }
  }
  return(df)
}

# Apply function

pop_deaths_2017 <- pop_deaths_2017 %>%
  group_by(sex, age) %>%
  do(adjust_n_deaths(.)) %>%
  ungroup()

# Correct n_deaths > n_census

adjust_n_census <- function(df, max_iter = 5, tol = 1e-5) {
  df$n_census_original <- df$n_census  # Store the original value of n_census

  
  iteration <- 1
  
  # Iterative process to adjust the values of n_census
    while (iteration <= max_iter) {
    changes <- FALSE
    
    # Sequentially evaluate all rows to apply cumulative adjustments
        for (i in 2:nrow(df)) {
      if (!is.na(df$n_census[i]) && df$n_census[i] < df$n_deaths[i]) {
        # Calculate the average between the current value and the adjusted value from the previous row
                new_average <- (df$n_census[i - 1] + df$n_census[i]) / 2
        
        # Perform the adjustment only if the new average is sufficiently different from the current value
        if (abs(new_average - df$n_census[i]) > tol) {
          df$n_census[i] <- new_average
          df$n_census[i - 1] <- new_average  # Also adjust the value in the previous row
          changes <- TRUE  # Indicate that a modification has been made
        }
      }
    }
    
    # Exit the loop if no changes were made
        if (!changes) {
      break
    }
    
    iteration <- iteration + 1
  }
  
  # Apply the final set of corrections
  
  df <- df %>%
    mutate(
      # Prevent restoring n_census_original values unless n_deaths actually exceeds n_census
      n_census = ifelse(!is.na(lead(n_deaths)) & lead(n_deaths) > lead(n_census), n_census_original, n_census),
      # Update n_census according to n_deaths if n_deaths remains higher
      n_census = ifelse(!is.na(n_census) & n_census < n_deaths, n_deaths, n_census),
      # Prevent replacing NA with zero; keep the original n_census values instead
      n_census = ifelse(is.na(n_census), n_census_original, n_census)
    )
  
  return(df)
}

# Apply function

pop_deaths_2017 <- pop_deaths_2017 %>%
  group_by(sex, age) %>%
  do(adjust_n_census(., max_iter = 5)) %>%
  ungroup() %>% 
  select(-n_census_original)

# Save dataset for data visualization

saveRDS(pop_deaths_2017, file = glue('{cnst$path_out}/pop_deaths_all_years_ed_2017.rds'))

# Group by 5 categories of education

pop_deaths_2017 <- pop_deaths_2017 %>% 
  mutate(years_ed = case_when(
    years_ed == 0 ~ 0,
    years_ed >= 1 & years_ed <= 6 ~ 1,
    years_ed >= 7 & years_ed <= 12 ~ 2,
    years_ed >= 13 & years_ed <= 16 ~ 3,
    years_ed >= 17 & years_ed <= 20 ~ 4,
    TRUE ~ NA
  )) %>% 
  group_by(age, sex, years_ed) %>% 
  summarise(n_census = sum(n_census),
            n_deaths = sum(n_deaths))

## Save final datasets ---- 

saveRDS(pop_deaths_2017, file = glue('{cnst$path_out}/pop_deaths_2017.rds'))

# Save census by sex

# Male census distribution

male_census <- pop_deaths_2017 %>%
  filter(sex == 1) %>% 
  select(-n_deaths)

saveRDS(male_census, file = glue('{cnst$path_out}/male_pop_deaths_2017.rds'))

# Female census distribution

female_census <- pop_deaths_2017 %>%
  filter(sex == 2)%>% 
  select(-n_deaths)

saveRDS(female_census, file = glue('{cnst$path_out}/female_pop_deaths_2017.rds'))

# Save deaths counts by sex 

# Male deaths counts distribution

male_deaths <- pop_deaths_2017 %>%
  filter(sex == 1) %>% 
  select(-n_census)

saveRDS(male_deaths, file = glue('{cnst$path_out}/male_deaths_2017.rds'))

# Female deaths counts distribution

female_deaths <- pop_deaths_2017 %>%
  filter(sex == 2)%>% 
  select(-n_census)

saveRDS(female_deaths, file = glue('{cnst$path_out}/female_deaths_2017.rds'))

rm(list = ls())
