HPC_summ <- function(x # index of run
){ 
  # Processing raw model runs
  #' input: index of model run 
  #' process: reads in the raw model run, aggregates over entire simulation per age group 
  #' output: data frame with one row per age group
  #'
  # read in rds file
  df <- readRDS(paste0(HPCpath, 'HPC/raw_modelrun_', format(x, scientific = FALSE),'.rds'))
  
  # get type of vaccine scenario and add doses
  RTSS <- df$RTSS[1]
  
  if(RTSS=='none'){
    df <- df |> rowwise() |>
      mutate(dose1 = 0,
             dose2 = 0,
             dose3 = 0,
             dose4 = 0) |> ungroup()
  }
  
  if(RTSS=='EPI' | RTSS == 'hybrid'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_epi_dose_1,
             dose2 = n_pev_epi_dose_2,
             dose3 = n_pev_epi_dose_3,
             dose4 = n_pev_epi_booster_1,) |> ungroup()
  }
  
  if(RTSS=='SV'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_mass_dose_1,
             dose2 = n_pev_mass_dose_2,
             dose3 = n_pev_mass_dose_3,
             dose4 = n_pev_mass_booster_1,) |> ungroup()
  }
  
  if(RTSS=='SVmass+EPI'| RTSS == 'SVmass+hybrid' | RTSS == 'mass+EPI'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_epi_dose_1 + n_pev_mass_dose_1,
             dose2 = n_pev_epi_dose_2 + n_pev_mass_dose_2,
             dose3 = n_pev_epi_dose_3 + n_pev_mass_dose_3,
             dose4 = n_pev_epi_booster_1 + n_pev_mass_booster_1) |> ungroup()
  }
  
  # Summarize the data 
  output <- df |> 
    # filter to get outcomes over the first 15 years  (where vaccination begins at 1 year)
    filter(year < 16) |>
    # combine 0-0.25 and 0.25-1 age groups 
    mutate(n_0_365 = n_0_91.25 + n_91.25_365,
           n_inc_clinical_0_365 = n_inc_clinical_0_91.25 + n_inc_clinical_91.25_365,
           n_inc_severe_0_365 = n_inc_severe_0_91.25 + n_inc_severe_91.25_365) |>
    group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, seasonality, speciesprop,
             ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    # mean of n in each age group
    mutate_at(vars(n_0_365, n_365_730:n_34675_36500, n_730_3650,
                   n_detect_730_3650), mean, na.rm = TRUE) |>
    # sum of cases and vaccine doses
    mutate_at(vars(n_inc_severe_0_365, n_inc_severe_365_730:n_inc_severe_34675_36500,
                   n_inc_clinical_0_365, n_inc_clinical_365_730:n_inc_clinical_34675_36500,
                   n_treated, n_infections,
                   starts_with('dose')), sum, na.rm = TRUE) |>
    # select relevant variables 
    dplyr::select(n_0_365, n_365_730:n_7300_7665, # 1 year age groups
                  # n_0_1825:n_34675_36500, # 5 year age groups
                  n_inc_clinical_0_365, n_inc_clinical_365_730:n_inc_clinical_7300_7665,# 1 year age groups
                  # n_inc_clinical_0_1825:n_inc_clinical_34675_36500, # 5 year age groups
                  n_inc_severe_0_365, n_inc_severe_365_730:n_inc_severe_7300_7665,# 1 year age groups
                  # n_inc_severe_0_1825:n_inc_severe_34675_36500, # 5 year age groups 
                  n_detect_730_3650,
                  n_730_3650, starts_with('dose'),
                  n_treated, n_infections) |>
    distinct()|>
    ungroup() |>
    # to long age groups -- here I'm pivoting only the 1 year age groups 
    pivot_longer(cols = c(n_0_365:n_7300_7665,
                          n_inc_clinical_0_365:n_inc_clinical_7300_7665,
                          n_inc_severe_0_365:n_inc_severe_7300_7665),
                 names_to = c('age'), values_to = c('value')) |>
    mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
           inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
           inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
           age = gsub('n_inc_clinical_', '', age),                         # combining age vars
           age = gsub('n_inc_severe_', '', age),
           age = gsub('n_', '', age),
           age = gsub('_', '-', age)) |> 
    group_by(age) |> #, year, month
    select(-c(value)) |>
    mutate_at(vars(n:inc_severe), sum, na.rm = TRUE) |> # consolidate
    distinct() |> ungroup() |>
    # cleaning up age group labels
    separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) |>
    mutate(age_lower = as.numeric(age_lower)/365,
           age_upper = as.numeric(age_upper)/365,
           age_grp = paste(age_lower, age_upper, sep = '-'),
           age_grp = factor(age_grp, levels =c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                               '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100')),
           age_cuts = ifelse(age_grp %in% c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'), '1y split','5y split'),
           inc = inc_clinical / n, #* 1000, # incidence
           sev = inc_severe / n, #* 1000, # severe incidence
           # get incidence per 1000 people fully vaccinated (i.e. got 3 doses)
           inc_per1000vax = ifelse(dose3 != 0, inc/dose3*1000, 0),
           sev_per1000vax = ifelse(dose3 != 0, sev/dose3*1000, 0),
           cases = inc_clinical,
           sev_cases = inc_severe)
  
  output$file <- paste0("./HPC/raw_modelrun_", x, ".rds")
  
  return(output)
}
# ggplot(clean_summ) + geom_col(aes(x = age_grp, y = sev_cases))+ labs(y = 'incidence per 1000 population and 1000 doses')

HPC_summ_5y <- function(x # index of run
){ 
  # Processing raw model runs - with 5 year age groups 
  #' input: index of model run 
  #' process: reads in the raw model run, aggregates over entire simulation per age group 
  #' output: data frame with one row per age group (by 5 years)
  #'
  # read in rds file
  df <- readRDS(paste0(HPCpath, 'HPC/raw_modelrun_',format(x, scientific = FALSE),'.rds'))
  
  # get type of vaccine scenario and add doses
  RTSS <- df$RTSS[1]
  
  if(RTSS=='none'){
    df <- df |> rowwise() |>
      mutate(dose1 = 0,
             dose2 = 0,
             dose3 = 0,
             dose4 = 0) |> ungroup()
  }
  
  if(RTSS=='EPI' | RTSS == 'hybrid'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_epi_dose_1,
             dose2 = n_pev_epi_dose_2,
             dose3 = n_pev_epi_dose_3,
             dose4 = n_pev_epi_booster_1,) |> ungroup()
  }
  
  if(RTSS=='SV'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_mass_dose_1,
             dose2 = n_pev_mass_dose_2,
             dose3 = n_pev_mass_dose_3,
             dose4 = n_pev_mass_booster_1,) |> ungroup()
  }
  
  if(RTSS=='SVmass+EPI'| RTSS == 'SVmass+hybrid' | RTSS == 'mass+EPI'){
    df <- df |> rowwise() |>
      mutate(dose1 = n_pev_epi_dose_1 + n_pev_mass_dose_1,
             dose2 = n_pev_epi_dose_2 + n_pev_mass_dose_2,
             dose3 = n_pev_epi_dose_3 + n_pev_mass_dose_3,
             dose4 = n_pev_epi_booster_1 + n_pev_mass_booster_1) |> ungroup()
  }
  
  # Summarize the data 
  output <- df |> 
    # filter to get outcomes over the first 15 years  (where vaccination begins at 1 year)
    filter(year < 16) |>
    # combine 0-0.25 and 0.25-1 age groups 
    mutate(n_0_365 = n_0_91.25 + n_91.25_365,
           n_inc_clinical_0_365 = n_inc_clinical_0_91.25 + n_inc_clinical_91.25_365,
           n_inc_severe_0_365 = n_inc_severe_0_91.25 + n_inc_severe_91.25_365) |>
    group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, seasonality, speciesprop,
             ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    # mean of n in each age group
    mutate_at(vars(n_0_365, n_365_730:n_34675_36500, n_730_3650,
                   n_detect_730_3650), mean, na.rm = TRUE) |>
    # sum of cases and vaccine doses
    mutate_at(vars(n_inc_severe_0_365, n_inc_severe_365_730:n_inc_severe_34675_36500,
                   n_inc_clinical_0_365, n_inc_clinical_365_730:n_inc_clinical_34675_36500,
                   n_treated, n_infections,
                   starts_with('dose')), sum, na.rm = TRUE) |>
    # select relevant variables 
    dplyr::select(#n_0_365, n_365_730:n_7300_7665, # 1 year age groups
      n_0_1825:n_34675_36500, # 5 year age groups
      # n_inc_clinical_0_365, n_inc_clinical_365_730:n_inc_clinical_7300_7665,# 1 year age groups
      n_inc_clinical_0_1825:n_inc_clinical_34675_36500, # 5 year age groups
      # n_inc_severe_0_365, n_inc_severe_365_730:n_inc_severe_7300_7665,# 1 year age groups
      n_inc_severe_0_1825:n_inc_severe_34675_36500, # 5 year age groups
      n_detect_730_3650,
      n_730_3650, starts_with('dose'),
      n_treated, n_infections) |>
    distinct()|>
    ungroup() |>
    # to long age groups -- here I'm pivoting only the 5 year age groups 
    pivot_longer(cols = c(n_0_1825:n_34675_36500,
                          n_inc_clinical_0_1825:n_inc_clinical_34675_36500,
                          n_inc_severe_0_1825:n_inc_severe_34675_36500),
                 names_to = c('age'), values_to = c('value')) |>
    mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
           inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
           inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
           age = gsub('n_inc_clinical_', '', age),                         # combining age vars
           age = gsub('n_inc_severe_', '', age),
           age = gsub('n_', '', age),
           age = gsub('_', '-', age)) |> 
    group_by(age) |> #, year, month
    select(-c(value)) |>
    mutate_at(vars(n:inc_severe), sum, na.rm = TRUE) |> # consolidate
    distinct() |> ungroup() |>
    # cleaning up age group labels
    separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) |>
    mutate(age_lower = as.numeric(age_lower)/365,
           age_upper = as.numeric(age_upper)/365,
           age_grp = paste(age_lower, age_upper, sep = '-'),
           age_grp = factor(age_grp, levels =c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                               '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100')),
           age_cuts = ifelse(age_grp %in% c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'), '1y split','5y split'),
           inc = inc_clinical / n, #* 1000, # incidence
           sev = inc_severe / n, #* 1000, # severe incidence
           # get incidence per 1000 people fully vaccinated (i.e. got 3 doses)
           inc_per1000vax = ifelse(dose3 != 0, inc/dose3*1000, 0),
           sev_per1000vax = ifelse(dose3 != 0, sev/dose3*1000, 0),
           cases = inc_clinical,
           sev_cases = inc_severe)
  
  output$file <- paste0("./HPC/raw_modelrun_5y", x, ".rds")
  
  return(output)
}


#' # HPC_annual <- function(x){
#'   # Processing raw model runs
#'   #' input: index of model run 
#'   #' process: reads in the raw model run, aggregates by year 
#'   #' output: data frame with one row per year 
#'   HPCpath <- "M:/Kelly/massvacc/"
#'   # read in rds file
#'   df <- readRDS(paste0(HPCpath, 'HPC/raw_modelrun_', format(x, scientific = FALSE), '.rds'))
#'   
#'   output <- df |>
#'     select(-c(n_treated, n_infections, n_use_net, timestep, month)) |>
#'     rowwise() |>
#'     # sum to get values for age group that was vaccinated with mass vax (everyone, school-aged)
#'     mutate(n_1825_36500 = sum(across(n_1825_3650:n_34675_36500)),
#'            n_inc_clinical_1825_36500 = sum(across(n_inc_clinical_1825_3650:n_inc_clinical_34675_36500)),
#'            n_inc_severe_1825_36500 = sum(across(n_inc_severe_1825_3650:n_inc_severe_34675_36500)),
#'            # school-aged
#'            n_1825_5475 = sum(across(n_1825_3650:n_5110_5475)),
#'            n_inc_clinical_1825_5475 = sum(across(n_inc_clinical_1825_3650:n_inc_clinical_5110_5475)),
#'            n_inc_severe_1825_5475 = sum(across(n_inc_severe_1825_3650:n_inc_severe_5110_5475))) |>
#'     ungroup() |>
#'     # grouping by year to get values over time 
#'     group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, year, seasonality, speciesprop,
#'              ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
#'     # get n, inc clincial, inc severe, and n_detect for each year 
#'     mutate_at(vars(n_730_3650, n_detect_730_3650,
#'                    n_inc_clinical_0_91.25, n_inc_clinical_0_1825, n_inc_clinical_1825_36500, n_inc_severe_1825_36500,
#'                    n_inc_clinical_1825_5475, n_inc_severe_1825_5475,
#'                    n_0_1825, n_0_91.25, n_1825_36500, n_1825_5475), mean, na.rm = TRUE) |>
#'     distinct() |>
#'     # calculate prevalence per year
#'     mutate(prev_byyear = n_detect_730_3650 / n_730_3650) |>
#'     # calculate incidence per year 
#'     mutate(inc_0_1825 = n_inc_clinical_0_1825 / n_0_1825 * 1000, # inc per 1000 pop
#'            inc_0_91.25 = n_inc_clinical_0_91.25 / n_0_91.25 * 1000, # inc per 1000 pop
#'            # everyone
#'            inc_1825_365000 = n_inc_clinical_1825_36500 / n_1825_36500 * 1000,
#'            sev_1825_365000 = n_inc_severe_1825_36500 / n_1825_36500 * 1000,
#'            #school-aged 
#'            inc_1825_5475 = n_inc_clinical_1825_5475 / n_1825_5475 * 1000,
#'            sev_1825_5475 = n_inc_severe_1825_5475 / n_1825_5475 * 1000
#'            )
#'     
#'   output$file <- paste0("./HPC/raw_modelrun_", x, ".rds")
#'   return(output)
#' }

# Function to run HPC_annual and then clean/save 
process_runs_byyr <- function(x){
  # Processing raw model runs
  #' input: index of model run 
  #' process: reads in the raw model run, aggregates by year 
  #' output: data frame with one row per year 
  HPCpath <- "M:/Kelly/massvacc/"
  # read in rds file
  df <- readRDS(paste0(HPCpath, 'HPC/raw_modelrun_', format(x, scientific = FALSE), '.rds'))
  
  output <- df |>
    select(-c(n_treated, n_infections, n_use_net, timestep, month)) |>
    rowwise() |>
    # sum to get values for age group that was vaccinated with mass vax (everyone, school-aged)
    mutate(n_1825_36500 = sum(across(n_1825_3650:n_34675_36500)),
           n_inc_clinical_1825_36500 = sum(across(n_inc_clinical_1825_3650:n_inc_clinical_34675_36500)),
           n_inc_severe_1825_36500 = sum(across(n_inc_severe_1825_3650:n_inc_severe_34675_36500)),
           # school-aged
           n_1825_5475 = sum(across(n_1825_3650:n_5110_5475)),
           n_inc_clinical_1825_5475 = sum(across(n_inc_clinical_1825_3650:n_inc_clinical_5110_5475)),
           n_inc_severe_1825_5475 = sum(across(n_inc_severe_1825_3650:n_inc_severe_5110_5475))) |>
    ungroup() |>
    # grouping by year to get values over time 
    group_by(ID, scenario, drawID, EIR, warmup, sim_length, population, pfpr, year, seasonality, speciesprop,
             ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    # get n, inc clincial, inc severe, and n_detect for each year 
    mutate_at(vars(n_730_3650, n_detect_730_3650,
                   n_inc_clinical_0_91.25, n_inc_clinical_0_1825, n_inc_clinical_1825_36500, n_inc_severe_1825_36500,
                   n_inc_clinical_1825_5475, n_inc_severe_1825_5475,
                   n_0_1825, n_0_91.25, n_1825_36500, n_1825_5475), mean, na.rm = TRUE) |>
    distinct() |>
    # calculate prevalence per year
    mutate(prev_byyear = n_detect_730_3650 / n_730_3650) |>
    # calculate incidence per year 
    mutate(inc_0_1825 = n_inc_clinical_0_1825 / n_0_1825 * 1000, # inc per 1000 pop
           inc_0_91.25 = n_inc_clinical_0_91.25 / n_0_91.25 * 1000, # inc per 1000 pop
           # everyone
           inc_1825_365000 = n_inc_clinical_1825_36500 / n_1825_36500 * 1000,
           sev_1825_365000 = n_inc_severe_1825_36500 / n_1825_36500 * 1000,
           #school-aged 
           inc_1825_5475 = n_inc_clinical_1825_5475 / n_1825_5475 * 1000,
           sev_1825_5475 = n_inc_severe_1825_5475 / n_1825_5475 * 1000
    )
  
  output$file <- paste0("./HPC/raw_modelrun_", x, ".rds")
  
  print(paste0('Saving run ', x))
  saveRDS(output, paste0(HPCpath, 'HPC_summbyyr/run_summ_', x, '_', output$int_ID[1], '.rds'))
}

get_annu_cr_int <- function(x){
  # input: scenario number from annual output (processed with HPC_annual)
  # process: combining all files for each scenario and getting median, and 95% CrI
  # output: aggregated data frame with 1 row with intervention ID and median/95% CrI
  
  files <- files |>
    # filter by scenario group 
    filter(group == x) 
  
  dat_list <- lapply(files$file, function (x) readRDS(x))
  
  # combine all draws for the specific scenario
  scen_group <- bind_rows(dat_list) |>
    group_by(warmup, sim_length, population, pfpr, year, seasonality, speciesprop, ITN, ITNuse, ITNboost,
             resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    summarise(across(c(prev_byyear, inc_0_1825, inc_0_91.25, inc_1825_365000, sev_1825_365000, inc_1825_5475, sev_1825_5475), 
                     list(lower = ~quantile(.x, 0.025), median = ~quantile(.x, 0.5), upper = ~quantile(.x, 0.975)),
                     .names = "{.col}_{.fn}")) 
  
  print(paste0('Run ',x))
  return(scen_group)
}


get_cr_int <- function(x, df){
  # input: scenario number from overall output (processed with HPC_summ)
  # process: get median and 95% CrI for each scenario/age group
  # output: aggregated data frame with 1 row for each intervention ID/scenario per age group and median/95% CrI columns for each variable 
  
  dalyout <- df %>%
    filter(group == x) |>
    group_by(age_grp, age_cuts, warmup, sim_length, population, pfpr, seasonality, speciesprop, ITN, ITNuse, ITNboost,
             resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    summarise(across(c(inc_clinical, inc_severe, inc, sev, inc_per1000vax, sev_per1000vax, #  incidence rate variables
                       cases, sev_cases,
                       mortality_rate, deaths, cases_lower, cases_upper, deaths_lower, deaths_upper, yll, yll_lower, yll_upper,
                       yld, yld_lower, yld_upper, daly, daly_upper, daly_lower, n_0_1825, n_91.25_1825, u5_cases, u5_severe,
                       u5_dalys, daly_baseline, cases_baseline, severe_baseline, deaths_baseline, u5_dalys_baseline, u5_cases_baseline,
                       u5_severe_baseline, dalys_averted, u5_dalys_averted, cases_averted, u5_cases_averted, deaths_averted, 
                       severe_averted, u5_severe_averted, dalys_avertedper1000vax, u5_dalys_avertedper1000vax, cases_avertedper1000vax, 
                       u5_cases_avertedper1000vax, deaths_avertedper1000vax, severe_avertedper1000vax, u5_severe_avertedper1000vax), 
                     list(lower = ~quantile(.x, 0.025), median = ~quantile(.x, 0.5), upper = ~quantile(.x, 0.975)),
                     .names = "{.col}_{.fn}") ) 
  
  return(dalyout)
  
}


# process by month for dose plot ------------------------------------------
# summarize individual model runs
summ_fordose <- function(x){
  raw <- readRDS(paste0(HPCpath, 'HPC_5355/HPC/raw_modelrun_', format(x, scientific = FALSE), '.rds'))
  
  if(raw$RTSS[1]!='none'){
    df <- raw |> 
      filter(timestep > 0 & timestep < 365*3) |>
      # mutate(int_ID = paste(pfpr, seasonality, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, sep = '_')) |>
      select(month, starts_with('n_pev'), n_detect_730_3650, n_730_3650, seasonality, RTSS, RTSScov, RTSSage, 
             RTSSrounds, fifth) |> #, MDAcov, MDAtiming, int_ID
      pivot_longer(starts_with('n_pev'), names_to = "dose", values_to = "count") |>
      group_by(month, dose) |>
      summarize(dosecount = sum(count),
                n_detect = mean(n_detect_730_3650),
                n = mean(n_730_3650),
                pfpr = n_detect / n,
                seasonality = seasonality, 
                RTSS = RTSS, 
                RTSSage = RTSSage, 
                RTSScov = RTSScov, 
                fifth = fifth,
                RTSSrounds = RTSSrounds,
                int_ID = int_ID) |> distinct()
    
    print(paste0('Saving run ', x))
    saveRDS(df, paste0(HPCpath,'HPC_5355/HPC_bydose/run_bydose_', x, '_', df$int_ID[1], '.rds'))
  }
}


get_cr_int_doses <- function(x){
  # input: scenario number from monthly output for doses (processed with process_for_doses)
  # process: combining all files for each scenario and getting median, and 95% CrI
  # output: aggregated data frame with 1 row with intervention ID and median/95% CrI
  doses_cr <- out %>%
    filter(group == x) |>
    group_by(pfpr, seasonality, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
    summarise(across(c(dosecount, n_detect, n), 
                     list(lower = ~quantile(.x, 0.025), median = ~quantile(.x, 0.5), upper = ~quantile(.x, 0.975)),
                     .names = "{.col}_{.fn}") ) 
  
  return(doses_cr)
}