##############################################################################
### Dose-response model for controlled chamber 5% NaCl toilet experiments  ###
### Reads in particle counts and estimates an inhaled dose                 ###
### From this calculates a probability of infection                        ###
###                                                                        ###
### Ciara A. Higham March 2023 - JHI paper                                 ###
##############################################################################

###############################
###        Pre-lims         ###
###############################

### install packages ###

#install.packages("EnvStats")
#install.packages("ggplot2")
#install.packages("truncnorm")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("reshape2")

### load packages ###
library(ggplot2)
library(EnvStats)
library(truncnorm)
library(tidyr)
library(dplyr)
library(reshape2)

### set seed ###
set.seed(1234)

###############################
###       Functions         ###
###############################

# Hockey-stick distribution for viral load (gc/L)
# detailed in McBride (2005)
hockeyStickDistribution <- function(min, max, median, percentile, N) {
  p = percentile/100 # P-th percentile
  q <- 1 - p
  s <- p - 0.5
  h1 <- 1 / (median - min)
  x_p <- 0.5 * (median + max + 1 / h1 - sqrt((max - median)^2 + 1 / (h1^2) + (median * (2 - 8 * q) + max * (2 - 8 * s)) / h1))
  h2 <- (2 * q) / (max - x_p)
  
  # Linear interpolation for probability densities in between x_0, x_50, x_p and x_100
  # N*1000 values to choose from for N samples
  linear_interpolation <- function(x_values, probability_density_values, N) {
    result_x <- runif(N*1000, min(x_values), max(x_values))
    result_prob_density <- numeric(N*1000)
    
    for (i in 1:(length(x_values) - 1)) { # gradient and intercept between points
      m <- (probability_density_values[i + 1] - probability_density_values[i]) / (x_values[i + 1] - x_values[i])
      b <- probability_density_values[i] - m * x_values[i]
      
      # corresponding probability densities
      result_prob_density[which(result_x >= x_values[i] & result_x <= x_values[i + 1])] <- m * result_x[which(result_x >= x_values[i] & result_x <= x_values[i + 1])] + b
    }
    return(data.frame(x = result_x, prob_density = result_prob_density))
  }
  
  x_values <- c(min, median, x_p, max)
  pd <- c(0, h1, h2, 0)
  
  # F(x) i.e. probability distribution for x_values
  pf_function <- linear_interpolation(x_values, pd, N)
  
  # sample from the x values that have probability densities prob_density
  path_conc_samples <- sample(pf_function$x, size = N, prob = pf_function$prob_density, replace = TRUE)
  
  return(path_conc_samples)
}

# size distribution for initial size of measured particles
# see supplementary info for details
initialSizeDist <- function(size, s, RH, v) {
  t=(-v+sqrt(v^2-4*(-s)*9.81*0.5))/(2*9.81*0.5) # using suvat kinematics equations and quadratic formula

  size_dist <- sqrt((size*1e-6)^2 + 4 * t * 1.1e-9 * (1 - RH))/1e-6 # maximum initial size, Bozic 2021 paper
  
  return(size_dist)
}

###############################
###     Input parameters    ###
###############################

# comment out Male/Female, S2/S1, counter_a/counter_b, SARS/Norovirus

#gender = "Male"
gender = "Female"

arrangement = "S2"
#arrangement = "S1"

counter_loc = "counter_a"
#counter_loc = "counter_b"

#pathogen = "SARS-CoV-2"
pathogen = "Norovirus"

N = 333 # number of simulations - 333*3 = 999 points for each t_0 and ACH

# Generate a pdf for toilet S2 occupancy times
if (gender == "Male"){
  # exposure time for males
  # data from Henning 1977
  male_time_periods <- seq(30, 900, by=30)
  male_frequencies <- c(3, 19, 10, 11, 7, 7, 9, 5, 7, 5, 4, 1, 0, 1, 4, 0, 0, 2, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
  
  t_male <- sample(male_time_periods, size=N, replace = TRUE, prob = male_frequencies)
  
  t_dur = t_male
  
} else if(gender == "Female"){
  # exposure time for female
  # data from Henning 1977
  female_time_periods <- seq(30, 270, by = 30)
  female_frequencies <- c(1, 27, 32, 26, 8, 1, 2, 2, 3)
  
  t_female <- sample(female_time_periods, size=N, replace = TRUE, prob = female_frequencies)
  t_dur = t_female
}

current_wd <- dirname(rstudioapi::getActiveDocumentContext()$path) # find the current wd

if(arrangement == "S2"){ # set the wd to where the data is located
  wd <- paste0(current_wd,"/raw_data/S2")
} else if(arrangement == "S1"){
  wd <- paste0(current_wd,"/raw_data/S1")
}

if(counter_loc == "counter_a"){
  wd <- paste0(wd, "/counter_a")
} else if(counter_loc == "counter_b"){
  wd <- paste0(wd, "/counter_b")
}

if(pathogen == "SARS-CoV-2"){
  
  percentile = 95
  sars_min_dens = 5e3 * 10^{3}  # minimum value in gc/L
  sars_max_dens = 10^{7.6} * 10^{3} # max value in gc/L
  sars_median_dens = 10^{5.1} * 10^{3} # median value in gc/L
  
  # create hockey stick distribution
  sars_2_conc <- hockeyStickDistribution(sars_min_dens, sars_max_dens, sars_median_dens, percentile, N*(25*60 - max(t_dur)+10)*3*3/(10))
  
  rho = sars_2_conc # viral load
  
  k = 1e-2/18.54 # k in exponential model
  
  doseResponse <- function(dose){ # exponential dose reponse model
    response <- 1 - exp(-dose * k)
    return(response)
  }
  
  # set wd for saving data
  save_wd <- paste0(current_wd,"/dose_response_output/sars_cov_2/")
  
} else if (pathogen == "Norovirus"){
  
  percentile = 95
  density_faeces = 1.06 * 10^3 # density of faeces
  norovirus_max_dens = 10^{11.76} * density_faeces # gc/g *g/L
  norovirus_median_dens = 10^{8.93} * density_faeces # gc/g *g/L
  norovirus_min_dens = 10^{6.34} * density_faeces # gc/g *g/L
  
  # create hockey stick distribution
  norovirus_conc <- hockeyStickDistribution(norovirus_min_dens, norovirus_max_dens, norovirus_median_dens, percentile, N*(25*60 - max(t_dur)+10)*3*3/(10)) 
 
  rho = norovirus_conc
  
  # constants for fractional poisson
  alpha = 0.104
  beta = 32.3
  
  doseResponse <- function(dose){ # fractional poisson dose-resposne
    response <- 1 - (1+dose/beta)^(-alpha)
    return(response)
  }
  
  # set wd to save
  save_wd <- paste0(current_wd,"/dose_response_output/norovirus/")
}

setwd(wd) # where the raw_data is stored

#######################################
### Read in the experimental data ##### 
#######################################

# define file names and ACH values
file_names <- c("1.5ach-a.txt", "1.5ach-b.txt", "1.5ach-c.txt",
                "3ach-a.txt", "3ach-b.txt", "3ach-c.txt",
                "6ach-a.txt", "6ach-b.txt", "6ach-c.txt")

# ach values
ach_values <- c(1.5, 3, 6) 
file_data <- lapply(file_names, read.table, header = TRUE)

conc_dist_df <- data.frame() # data frame for c_{i,j,t}

# breathing rate
breath <- rtruncnorm(N*(25*60 - max(t_dur)+10)*3*3/(10), a = 0, b = Inf, mean = 1.2e-2, sd = 2.5e-3)

# f - factor for gc to PFU
f = runif(N*(25*60 - max(t_dur)+10)*3*3/(10), 100, 1000)

# volume of faeces in L
vol_faeces <- rtruncnorm(N*(25*60 - max(t_dur)+10)*3*3/(10), a = 0, b = Inf, mean = 1.82e-1, sd = 1.51e-1)

for (ach in ach_values) { # iterate over 1.5, 3 and 6 ACH
  
  temp_post_flush <- data.frame() # temp list to store the post-flush data frames
  
  # read in data file
  ach_files <- file_names[grep(paste0(ach, "ach-"), file_names)] # filter file to the corresponding ach in loop
  
  for (i in 1:3) { # loop through the replicates for each ach
    
    data <- file_data[[which(file_names == ach_files[i])]]
    
    pre_flush_rows <- data[892:921, ] # data 5 minutes pre-flush
    post_flush_rows <- data[922:nrow(data), ] # data post-flush
    
    colnames(post_flush_rows) <- c("0.3", "0.5", "1", "3", "5", "10") # column names to particle sizes
    
    bckg_conc_df <- apply(pre_flush_rows, 2, median) # assign median to data frame
    
    for (j in 1:6) { # remove background for each size
      post_flush_rows[, j] <- post_flush_rows[, j] - bckg_conc_df[j]
    }
    
    post_flush_rows[post_flush_rows < 0] <- 0 # if < 0 set to 0
    post_flush_rows$Time <- seq(1, nrow(post_flush_rows)) * 10 # time post-flush column
    
    max_time <- max(post_flush_rows$Time) # max time of post-flush data i.e. when counter stops running
    
    if (max_time > 1500) { # if max time is greater than 25 mins
      valid_rows <- post_flush_rows$Time <= 1500 # just keep data up to 25 mins
      post_flush_rows <- post_flush_rows[valid_rows, ]
    } else if (max_time < 1500) { # if the counter was stopped before 25 mins
      new_times <- seq(max_time + 10, 1500, by = 10) # 
      new_rows <- data.frame(Time = new_times)
      new_rows[, c("0.3", "0.5", "1", "3", "5", "10")] <- 0 # set the values after the max time but before 25 min to zero
      post_flush_rows <- rbind(post_flush_rows, new_rows)
    }
    
    post_flush_rows$ACH <- rep(ach, nrow(post_flush_rows)) # add ACH column
    post_flush_rows$Replicate <- rep(i, nrow(post_flush_rows)) # add a column indicating replicate number
    # melt so that there is a column for size
    post_flush_rows <- melt(post_flush_rows, id.vars = c("Replicate", "Time", "ACH"))
    
    colnames(post_flush_rows)[4] <- "Size" # name size col
    colnames(post_flush_rows)[5] <- "Conc" # name conc col
    
    # change size values to initial size (different due to evaporation)
    temp_post_flush <- bind_rows(temp_post_flush, post_flush_rows) # combine data for each replicate
  }
  conc_dist_df <- bind_rows(conc_dist_df, temp_post_flush) # combine each ACH frame
}

#######################################
### Dose-response using Monte Carlo ###
#######################################

final_data_df <- data.frame()  # final response data with all parameters

for (t_0_iter in seq(0, 25*60 - max(t_dur), by = 10)) { # iterate over the possible times for entering the bathroom post flush i.e. between t=10s and t=25mins
  t_exit = t_0_iter + t_dur # exit time
  drop_vol_ratio_df <- data.frame() # df for the total droplet volume inhaled between t_exit and t_0
  
  for (exit_time in t_exit){
    
    # subset of the data for t_0:t_exit
    temp_t_exit_df<- data.frame()
    temp_t_exit_df <- conc_dist_df[conc_dist_df$Time > t_0_iter & conc_dist_df$Time <= exit_time, ]
    
    temp_t_exit_df <- temp_t_exit_df %>%
      mutate(Size = as.numeric(as.character(Size)),) %>%  # convert size to numeric for use in initial size calc
      group_by(ACH, Replicate) %>%
      mutate(RH = runif(1, 0.4, 0.5), # set RH and v forgiven experiment
             v = runif(1, 1, 2)) %>% # but need to be the same for sizes within the experiment
      ungroup()
    
    # set the initial size of particles using the RH and initial size
    temp_t_exit_df$Initial_Size <- initialSizeDist(temp_t_exit_df$Size, s=0.2, RH = temp_t_exit_df$RH, v =temp_t_exit_df$v)
    
    # sum the concentrations across the occupancy times for a given size
    temp_t_exit_df <- temp_t_exit_df %>%
      group_by(ACH, Replicate, Size) %>%
      mutate(Conc_Sum = sum(Conc)) %>%
      ungroup()
    temp_t_exit_df <- unique(temp_t_exit_df[, c("ACH", "Replicate", "Size", "Conc_Sum", "RH", "v", "Initial_Size")])
    
    # multiply by cube of diameter
    temp_t_exit_df <- temp_t_exit_df %>%
      mutate(Conc_Droplet = Conc_Sum  * (as.numeric(Initial_Size) * 1e-6)^3) %>% 
      
      # and then sum concentrations across size for a given replicate
      group_by(ACH, Replicate) %>%
      mutate(Drop_Vol_Ratio = sum(Conc_Droplet))%>%
    ungroup()
    temp_t_exit_df <- unique(temp_t_exit_df[, c("ACH", "Replicate", "RH", "v", "Drop_Vol_Ratio")])
    
    # column for the time spent in toilet
    temp_t_exit_df$t_dur <- exit_time - t_0_iter
    
    # store the data for each ach
    drop_vol_ratio_df <- bind_rows(drop_vol_ratio_df, temp_t_exit_df)
  }
  
  # set the parameter values for each t_0 iteration
  temp_t_0_df <- drop_vol_ratio_df %>%
    mutate(
      RH = RH,
      v = v,
      Drop_Vol_Ratio = Drop_Vol_Ratio,
      Arrangement = arrangement,
      Gender = gender,
      Counter_Location = counter_loc,
      t_enter= t_0_iter
    ) %>%
    ungroup()
  
  # combination of all possible t_0
  final_data_df <- bind_rows(final_data_df, temp_t_0_df)
}

# calculate the dose and the response
# and store the parameter values
dose_response_df <- final_data_df %>%
  mutate(
    Breath_Rate = breath,
    Vol_Faeces = vol_faeces,
    rho = rho,
    f = f,
    Dose = Drop_Vol_Ratio * breath * rho * vol_faeces * (1000) * pi/ (6 * 6 * f * (1 + vol_faeces)),
    Response = doseResponse(Dose), # either exponential or beta-poisson depending on norovirus or sars
  )

setwd(save_wd) # change wd to where the output csv will be saved

# write to file
write.table(dose_response_df, paste0(tolower(gender),"_",tolower(arrangement),"_",tolower(counter_loc),".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
