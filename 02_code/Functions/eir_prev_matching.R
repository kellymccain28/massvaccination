#####################################################
#####################################################
# Have moved the matching to the end of the create_scenarios.R script
#####################################################
#####################################################

# Trying out pfpr and eir matching 
source("02_code/packages_data.R")
# source("02_code/Functions/generate_params.R")

## Model set up 
# year <- 365
# population <- 1000
# warmup <- 9 * year
# sim_length <- 9 * year
# human_population <- 10000


# calibration to average annual pfpr value for last 2 years of simulation 
annual_pfpr_summary <- function(x){
  x <- x[x$timestep > x$warmup,]
  x$year <- ceiling(x$timestep / 365)
  x <- x %>%
    filter(year == max(year) | year == max(year)-1)
  pfpr <- x$n_detect_730_3650 / x$n_730_3650
  year <- x$year
  tapply(pfpr, year, mean)
}

# test summary function 
# warmup <- 365
# sim_length <- 1460
# out <- run_simulation(warmup+sim_length, params)
# annual_pfpr_summary(out)

pr_match <- function(x, y){
  
  data <- readRDS(paste0(path, "03_output/baseline_parameters.rds"))[x,]
  params <- unlist(data$params, recursive = FALSE)
  params$timesteps <- data$sim_length + data$warmup
  
  # defining target as pfpr value in last 2 years of simulation 
  target <- rep(data$pfpr, 2)
  
  set.seed(1234)
  out <- calibrate(parameters = params,
                   target = target,
                   summary_function = annual_pfpr_summary,
                   tolerance = 0.001,
                   low = 0.1,
                   high = 1500)
  
  # store init_EIR results as an .rds file to be read in later
  PR <- data.frame(scenarioID = x, drawID = y)
  PR$starting_EIR <- out
  PR$ID <- data$ID
  
  print(paste0('Finished scenario ',x))
  saveRDS(PR, paste0('03_output/PrEIR/PRmatch_draws_', data$ID, '.rds'))
}



