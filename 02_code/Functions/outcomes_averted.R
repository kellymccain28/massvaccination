# input: data frame that has been processed with HPC_summ and combined -- dalyoutput in Processing_model_runs.R
# process: separates out the baseline scenarios and calculates clinical, severe cases, deaths, and DALYs averted
# output: data frame with 1 row per age group per intervention scenario per drawID with new columns for outcomes averted
# unique IDs: ID, int_ID, age_grp

outcomes_averted <- function(df){
  
  # df <- readRDS(paste0(HPCpath, "03_output/dalyoutput_draws.rds"))
  # df <- readRDS(paste0(HPCpath,'/HPC_summarized/aggregatedoversim_1_36.rds'))
  # create vars for childhood cases
  
  output <- df |> 
    ungroup()|>
    # select(-inc, -sev, -mortality_rate) |> # get rid of rate vars
    group_by(age_grp, ID, int_ID) |>
    # create vars for childhood cases
    mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '0-1', '91.25-1825'), n, 0), # u5 denominator
           n_91.25_1825 = ifelse(age == '91.25-1825', n, 0),
           u5_cases = ifelse(age %in% c('0-91.25','0-1', '91.25-1825'), cases, 0),
           u5_severe = ifelse(age %in% c('0-91.25','0-1', '91.25-1825'), sev_cases, 0),
           u5_dalys = ifelse(age %in% c('0-91.25', '0-1','91.25-1825'), daly, 0)) |>
    mutate_at(vars(n, n_0_1825, n_91.25_1825,
                   u5_cases, u5_severe, u5_dalys,
                   inc_clinical, inc_severe,
                   inc_per1000vax, sev_per1000vax,
                   cases, cases_lower, cases_upper, sev_cases,
                   deaths, deaths_lower, deaths_upper,
                   yll:daly_lower),
              sum, na.rm = T) |>  
    select(-age, -age_upper, -age_lower) #, -int_ID, -group

  # separate out baseline scenarios
  none <- output |>
    ungroup() |>
    filter(RTSS == 'none') |>
    rename(daly_baseline = daly,
           cases_baseline = cases,
           severe_baseline = sev_cases,
           deaths_baseline = deaths,
           
           u5_dalys_baseline = u5_dalys,
           u5_cases_baseline = u5_cases,
           u5_severe_baseline = u5_severe) |>
    
    select(file, ID, scenario, drawID, age_grp, daly_baseline, cases_baseline,
           severe_baseline, deaths_baseline, u5_dalys_baseline,
           u5_cases_baseline, u5_severe_baseline, dose3) |> distinct()
  
  # separate out non baseline scenarios and merge
  base_IDs <- unique(none$file)
  
  averted <- output |> filter(!(file %in% base_IDs)) |>
    left_join(none |> select(-file, -scenario, -dose3), by = c('ID', 'drawID', 'age_grp')) |>
    
    # calculate outcomes averted
    mutate(dalys_averted = daly_baseline - daly,
           u5_dalys_averted = u5_dalys_baseline - u5_dalys,
           cases_averted = cases_baseline - cases,
           u5_cases_averted = u5_cases_baseline - u5_cases,
           deaths_averted = deaths_baseline - deaths, 
           severe_averted = severe_baseline - sev_cases, 
           u5_severe_averted = u5_severe_baseline - u5_severe) |>
    
    # calculate outcomes averted per 1000 fully vaccinated people
    mutate(dalys_avertedper1000vax = dalys_averted/dose3 *1000,
           u5_dalys_avertedper1000vax = u5_dalys_averted/dose3*1000,
           cases_avertedper1000vax = cases_averted/dose3*1000,
           u5_cases_avertedper1000vax = u5_cases_averted/dose3*1000,
           deaths_avertedper1000vax = deaths_averted/dose3*1000,
           severe_avertedper1000vax = severe_averted/dose3*1000,
           u5_severe_avertedper1000vax = u5_severe_averted/dose3*1000)
  
  return(averted)
}
