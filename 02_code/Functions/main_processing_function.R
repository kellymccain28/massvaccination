# Function to process all of the different ways at once 
process_all <- function(pathtoscenarios, # path pointing to the scenarios .rds file to know how many there are 
                        HPCfolder # folder name 
                        ){
  scenarios <- readRDS(pathtoscenarios)
  index <- c(1:nrow(scenarios))
  ## Processing by month (saved in HPC_summbymonth) ----
  # walk(index, 
  #         get_rates_process, 
  #         time_div = 30.41667, 
  #         age_agg = TRUE, 
  #         HPCfolder = HPCfolder,
  #         sevoutcomes = FALSE) 
  # 
  # files <- list.files(path = paste0(HPCpath, HPCfolder, "/HPC_summbymonth/"), pattern = "run_summ_*", full.names = TRUE) 
  # summbymonth <- lapply(files, function (x) readRDS(x)) %>% bind_rows()
  # 
  # saveRDS(summbymonth, paste0(HPCpath, HPCfolder, '/summbymonth_draws.rds'))
  # 
  # # Credible intervals 
  # summbymonth_out <- summbymonth %>%
  #   group_by(t, int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, 
  # EPIbooster, EPIextra, massbooster_rep, MDA) |> 
  #   mutate(across(c(clinical:n, contains('dose'),
  #                   prevalence_2_10, prevalence_0_100), 
  #                 list(lower = ~quantile(.x, 0.025, na.rm = TRUE), 
  #                      median = ~quantile(.x, 0.5, na.rm = TRUE), 
  #                      upper = ~quantile(.x, 0.975, na.rm = TRUE)),
  #                 .names = "{.col}_{.fn}") ) %>% 
  #   distinct(across(-c(scenario, clinical:n, ID, drawID, EIR, dose1:file)))
  # 
  # 
  # saveRDS(summbymonth_out, paste0(HPCpath, HPCfolder, '/summbymonth.rds'))
  # 
  # print('Finished processing by month')
  
  
  ## Processing annually (saved in HPC_summbyyr) ----
  walk(index, get_rates_process, 
          time_div = 365, 
          age_agg = TRUE, 
          HPCfolder = HPCfolder,
          sevoutcomes = FALSE) 
  
  files <- list.files(path = paste0(HPCpath, HPCfolder, "/HPC_summbyyr/"), pattern = "run_summ_*", full.names = TRUE) 
  summbyyr <- lapply(files, function (x) readRDS(x)) %>% bind_rows()
  
  saveRDS(summbyyr, paste0(HPCpath, HPCfolder, '/summbyyr_draws.rds'))
  
  
  summbyyr_out <- summbyyr %>%
    group_by(t, int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, 
             EPIbooster, EPIextra, massbooster_rep, MDA) |> 
    mutate(across(c(clinical:n, contains('dose'),
                    prevalence_2_10, prevalence_0_100), 
                  list(lower = ~quantile(.x, 0.025, na.rm = TRUE), 
                       median = ~quantile(.x, 0.5, na.rm = TRUE), 
                       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                  .names = "{.col}_{.fn}") ) %>% 
    distinct(across(-c(scenario, clinical:n, ID, drawID, EIR, dose1:file)))
  
  
  saveRDS(summbyyr_out, paste0(HPCpath, HPCfolder, '/summbyyr.rds'))
  
  print('Finished processing by year')
  
  ## Process annually by age group (saved in HPC_summbyyrbyage) ----
  walk(index, get_rates_process, 
          time_div = 365, 
          age_agg = FALSE, 
          HPCfolder = HPCfolder, 
          sevoutcomes = TRUE) 
  
  files <- list.files(path = paste0(HPCpath, HPCfolder, "/HPC_summbyyrbyage/"), pattern = "run_summ_*", full.names = TRUE) 
  summbyyrbyage <- lapply(files, function (x) readRDS(x)) %>% bind_rows()
  
  saveRDS(summbyyrbyage, paste0(HPCpath, HPCfolder, '/summbyyrbyage_draws.rds'))
  
  
  summbyyrbyage_out <- summbyyrbyage %>%
    outcomes_averted() %>%
    group_by(age_lower, age_upper, t, int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, 
             EPIbooster, EPIextra, massbooster_rep, MDA) |> 
    mutate(across(c(clinical:n, contains('dose'),
                    prevalence_2_10, prevalence_0_100,
                    mortality_rate:sevper1000vax,
                    daly_baseline:severe_avertedper1000vax), 
                  list(lower = ~quantile(.x, 0.025, na.rm = TRUE), 
                       median = ~quantile(.x, 0.5, na.rm = TRUE), 
                       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                  .names = "{.col}_{.fn}") ) %>% 
    distinct(across(-c(scenario, clinical:prop_n, ID, drawID, EIR, dose1:file)))
  
  
  saveRDS(summbyyrbyage_out, paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))
  
  print('Finished processing annual by age group')
  
  ## Processing over entire simulation by age group (saved in HPC_summarized) ----
  walk(index, get_rates_process, 
          time_div = 7665, 
          age_agg = FALSE, 
          HPCfolder = HPCfolder,
          sevoutcomes = TRUE)
  
  files <- list.files(path = paste0(HPCpath, HPCfolder, "/HPC_summarized/"), pattern = "run_summ_*", full.names = TRUE) 
  summarized <- lapply(files, function (x) readRDS(x)) %>% bind_rows() %>%
    mutate(age_grp = paste(age_lower, age_upper, sep = '-'),
           age_cuts = ifelse(age_grp %in% c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'), '1y split',
                             ifelse(age_grp %in% '0-100', 'all', 
                                    ifelse(age_grp %in% '5-15', 'school-aged', '5y split')))
    )
  
  saveRDS(summarized, paste0(HPCpath, HPCfolder, '/summarized_draws.rds'))
  
  
  summarized_out <- summarized %>%
    outcomes_averted() %>%
    group_by(age_lower, age_upper, int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, 
             EPIbooster, EPIextra, massbooster_rep, MDA) |> 
    summarize(across(c(clinical:prop_n, contains('dose'),
                    prevalence_2_10, prevalence_0_100,
                    mortality_rate:sevper1000vax,
                    daly_baseline:severe_avertedper1000vax), 
                  list(lower = ~quantile(.x, 0.025, na.rm = TRUE), 
                       median = ~quantile(.x, 0.5, na.rm = TRUE), 
                       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                  .names = "{.col}_{.fn}") ) %>% 
    mutate(age_grp = paste(age_lower, age_upper, sep = '-'),
           age_cuts = ifelse(age_grp %in% c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'), '1y split',
                             ifelse(age_grp %in% '0-100', 'all', 
                                    ifelse(age_grp %in% '5-15', 'school-aged', '5y split'))))
  
  saveRDS(summarized_out, paste0(HPCpath, HPCfolder, '/summarized.rds'))
  
  print('Finished processing over simulation by age group')
  
}
