# Inequalities in life expectancy by educational quantile groups: A review of methods and approaches 
This repository accompanies **“Inequalities in life expectancy by educational quantile groups: A review of methods and approaches.”** 
An open-access link to the paper will be added here once available. 
The replication files are primarily implemented in the R statistical programming language, with part of the workflow executed in Stata.

## Overview
We review methodological and conceptual approaches to measuring educational inequalities in life expectancy using percentiles rather than coarse attainment categories. 
A Chilean case study demonstrates how different choices shape the resulting estimates and their comparability over time and across populations.
The repository contains the scripts required to reproduce the analytical workflow used in the study.

## Workflow
The analytical pipeline combines **R** and **Stata** scripts and should be run in the following order:

### 1. Data preparation and education distribution
These scripts prepare the analytical datasets and estimate the distribution of education:

- `1. Data_preparation.R`
- `2. Education_distribution.R`
- `2.1 Education_distribution_fxs.R`

### 2. Life table estimation based on alternative approaches
These scripts estimate life tables using the main analytical approaches implemented in the study:

- `3. Global_life_table.R`
- `4. Life_tables_compositional_adjustment.R`
- `4.1 Life_tables_fxs.R`
- `5. Life_tables_regression_based.R`
- `5.1 Regression_based_fxs.R`

### 3. Preparation of inputs for bounded mortality estimation
This script generates the input datasets required for the bounded mortality workflow in Stata:

- `6. Bounded_mortality_input.R`

### 4. Bounded mortality estimation in Stata
These Stata programs use the outputs from the previous R step to estimate mortality bounds by age, sex, and educational quintile:

- `7. Programs.do`
- `7.1 Cohort_mortality_q.do`
- `7.2 Lagged_cohort_mortality_q.do`
- `7.3 Period_mortality_q.do`

### 5. Life table estimation using bounded mortality inputs
This R script uses the mortality bounds generated in Stata to construct bounded life tables and corresponding life expectancy estimates:

- `8. Life_tables_bounded_mortality.R`

### 6. Simulations
These R scripts and Stata .do files reproduce the workflows required to simulate mortality rates under the log-linear and hump-shape modeling frameworks:

- `9. Simulation.R`
- `9.1 Compostional_adjustment_sim.R`
- `9.2 Regression_bases_sim.R`
- `10. Sim_q.do`
- `10.1. Bounded_mortality_sim.R`

### 7. Data visualization
This script produces the figures and visual summaries used in the empirical application:

- `11. Data_visualization.R`

## Prerequisites
To run the code locally you will need:
A working installation of R
Required R packages installed
Stata (for the bounded mortality estimation step)

## Running the Code
To run this code, do something like:
git clone https://github.com/cancer-inequalities/inequalities_in_life_expectancy_MethodsReview.git
Then execute the scripts sequentially according to the workflow described above.

## Repository structure
Data_folder/
    Input source data

Scripts/
    R and Stata scripts used for data preparation and analysis

Out/
    Output datasets generated during the workflow

Figures/
    Figures produced for the empirical application

## Data availability
The empirical application uses Chilean mortality and education data. 
The repository provides scripted pipelines and functions to reproduce the analyses with authorized inputs.

## Versioning
This repository hosts the pre-publication version of the code (v0.1.0). 
Please report bugs, request features, or share suggestions by opening an issue, or contact the corresponding author(s) via email.

## How to cite
Authors (Year). Inequalities in life expectancy by educational percentiles: a review of methods and approaches. Journal (forthcoming). DOI: TBD.

## License
This work is free. The code comes without any warranty, to the extent permitted by applicable law.


