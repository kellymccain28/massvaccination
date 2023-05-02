# Function to run simulation 
runsim <- function(x){ # x = scenario #
  
  year <- 365
  month <- year / 12
  
  # read in selected scenario
  data <- readRDS("03_output/parameters_torun.rds")[x,]
  # read in collated EIR estimates 
  match <- readRDS("03_output/PrEIR/EIRestimates.rds") |> select(-scenarioID)
 
  # EIR / prev match from "PfPR_EIR_match.R"
  data <- data |> left_join(match, by = c("drawID", "ID"))
  
  # Set EIR equilibrium ----------
  params <- set_equilibrium(unlist(data$params, recursive = F), as.numeric(data$starting_EIR))
  
  # run simulation ----------
  set.seed(123)
  
  output <- run_simulation(
    timesteps = data$warmup + data$sim_length,
    # correlations = correlations,
    parameters = params) |>
    
    # add vars to output
    mutate(ID = data$ID,
           scenario = data$scenarioID,
           drawID = data$drawID,
           EIR = data$starting_EIR,
           warmup = data$warmup,
           sim_length = data$sim_length,
           population = data$population,
           pfpr = data$pfpr,
           # timestep = timestep - data$warmup,
           seasonality = data$seas_name,
           speciesprop = paste(data$speciesprop, sep = ",", collapse = ""),
           ITN = data$ITN,
           ITNuse = data$ITNuse,
           ITNboost = data$ITNboost,
           resistance = data$resistance,
           IRS = data$IRS,
           treatment = data$treatment,
           SMC = data$SMC,
           RTSS = data$RTSS,
           RTSScov = data$RTSScov,
           RTSSage = data$RTSSage,
           RTSSrounds = data$RTSSrounds,
           fifth = data$fifth) |>
    ungroup() |>
    # filter(timestep > 0) |> # remove warmup period
    
    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) |>
    
    # keep only necessary variables
    dplyr::select(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop,
                  ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, timestep,
                  starts_with("n_inc_severe"), starts_with("p_inc_severe"),
                  starts_with("n_pev"),
                  starts_with("n_inc"), starts_with("p_inc"),
                  starts_with("n_detect"), starts_with("p_detect"),
                  starts_with("n_"), -n_bitten, n_treated, n_infections) |>
    
    # take means of populations and sums of cases by month
    group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop,
             ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth) |>
    
    mutate_at(vars(n_0_1825:n_7300_7665, n_730_3650,
                   n_detect_730_3650, p_detect_730_3650), mean, na.rm = TRUE) |>
    mutate_at(vars(n_inc_severe_0_1825:p_inc_clinical_7300_7665,
                   n_treated, n_infections), sum, na.rm = TRUE) |>
    
    dplyr::select(timestep, n_0_1825:n_7300_7665,
                  n_inc_severe_0_1825:p_inc_clinical_7300_7665,
                  n_detect_730_3650, p_detect_730_3650,
                  n_730_3650,
                  n_treated, n_infections) |>
    distinct()
  
  
  # save output ----------
  # saveRDS(output, paste0(path, "./03_output/raw_modelrun_",x,".rds"))
  saveRDS(output, paste0(HPCpath, "HPC/raw_modelrun_", x, ".rds"))
  
}
