# ** Paper title ** 

Supporting code and data for **paper title **

The code is written in R version 4.3.2 and has been run using RStudio 2023.12.0.

## Repository contents
- `raw_data/`: a directory containing the raw particle concentration data as `.csv` files for
  - two scenarios: S1 vs S2
    - at particle counter locations: A and B
      - at 3 ventilation rates: 1.5, 3, 6 air changes per hour (ACH)
        - each with 3 replicates (a, b and c)
        
-  `particle_conc_plot.R`: R script that reads in the raw particle concentration data and plots a time series of the mean concentration across the three replicates. The cubicle scenario is plotted as a purple line and the no cubicle scenario is plotted as a green line. Standard errors are plotted as a shaded region. The particle sizes are split into further sub plots.
  
-  `dose_response_monte_carlo.R`: R script that reads in the raw particle concentration data and uses a stochastic Monte Carlo approach to quantify exposure to SARS-CoV-2 and norovirus. Then uses a dose-response model to calculate the risk of exposure as a fraction.

-  `processed_dose_response_output.zip`: `.zip` file containing the output data (`.csv` files from running `dose_response_monte_carlo.R`. If you do not have the time/wish to run `dose_response_monte_carlo.R` unzip the file and use this data to run the `infection_risk_plot.R` script.

-  `infection_risk_plot.R`: R script that reads in the output data from running `dose_response_monte_carlo.R`, located in the `dose_response_output/` directory and plots the risk of infection (%) for SARS-CoV-2 and norovirus at particle counter location A and B.


*Note: `particle_conc_plot.R` can be run independent of other scripts `dose_response_monte_carlo.R` must be run prior to running `infection_risk_plot.R`.*

## **`particle_conc_plot.R`**

Uncomment the relevant counter_loc in lines 29 - 30 to give two plots. If you do not have the relevant packages installed, uncomment lines 11-15 before running.

## `dose_response_monte_carlo.R`

Lines 91 - 101 contain variables e.g. male/female, scenario 1/2, counter location a/b. This should be run for each combination so that 16 `.csv` are written out to the `dose_response_output/` directory. If you do not have the relevant packages installed, uncomment lines 15 - 20 before running.

## `infection_risk_plot.R`

After running `dose_response_monte_carlo.R` this script is run to plot the infection risk using the data in the `dose_response_output/` directory. Lines 82 - 86 contain variables, counter location a/b and SARS-CoV-2/norovirus. This should be run for each combination so that 4 plots are created. If you do not have the relevant packages installed, uncomment lines 82 - 86 before running.
