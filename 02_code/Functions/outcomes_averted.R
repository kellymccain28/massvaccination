# input: data frame that has been processed with HPC_summ and combined -- dalyoutput in Processing_model_runs.R
# process: separates out the baseline scenarios and calculates clinical, severe cases, deaths, and DALYs averted
# output: data frame with 1 row per age group per intervention scenario per drawID with new columns for outcomes averted
# unique IDs: ID, int_ID, age_grp and year if the dataset is by year

outcomes_averted <- function(df, byyear = FALSE){
  
  # df <- readRDS(paste0(HPCpath, "03_output/dalyoutput_draws.rds"))
  # df <- readRDS(paste0(HPCpath,'/HPC_summarized/aggregatedoversim_1_36.rds'))
  # create vars for childhood cases
  if(byyear == FALSE){
    df <- df |>
      ungroup() |>
      group_by(age_grp, ID, int_ID)
    
    joinvars <- c('ID', 'drawID', 'age_grp')
    
    baseline_vars <- c('file', 'ID', 'scenario', 'drawID', 'age_grp', 'daly_baseline', 'cases_baseline',
                       'severe_baseline', 'deaths_baseline', 'u5_dalys_baseline',
                       'u5_cases_baseline', 'u5_severe_baseline', 'dose3')
  } else if(byyear == TRUE){
    df <- df |> 
      ungroup()|>
      group_by(age_grp, year, ID, int_ID)
    
    joinvars <- c('ID', 'drawID', 'age_grp', 'year')
    
    baseline_vars <- c('file', 'ID', 'scenario', 'drawID', 'age_grp', 'daly_baseline', 'cases_baseline',
                       'severe_baseline', 'deaths_baseline', 'u5_dalys_baseline',
                       'u5_cases_baseline', 'u5_severe_baseline', 'dose3', 'year')
  }
  # create vars for childhood cases
  output <- df %>%
    mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '0-1', '91.25-1825', '0-1825'), n, 0), # u5 denominator
         n_91.25_1825 = ifelse(age == '91.25-1825', n, 0),
         u5_cases = ifelse(age %in% c('0-91.25','0-1', '91.25-1825', '0-1825'), cases, 0),
         u5_severe = ifelse(age %in% c('0-91.25','0-1', '91.25-1825', '0-1825'), sev_cases, 0),
         u5_dalys = ifelse(age %in% c('0-91.25', '0-1','91.25-1825', '0-1825'), daly, 0)) |>
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
    filter(RTSS == 'none'& booster_rep !='annual') |># 
    rename(daly_baseline = daly,
           cases_baseline = cases,
           severe_baseline = sev_cases,
           deaths_baseline = deaths,
           
           u5_dalys_baseline = u5_dalys,
           u5_cases_baseline = u5_cases,
           u5_severe_baseline = u5_severe) |>
    
    select(all_of(baseline_vars)) |> distinct()
  
  # separate out non baseline scenarios and merge
  base_IDs <- unique(none$file)
  
  averted <- output |> filter(!(file %in% base_IDs)) |>
    left_join(none |> select(-file, -scenario, -dose3), by = joinvars) |> 
    
    # calculate outcomes averted
    mutate(dalys_averted = daly_baseline - daly,
           u5_dalys_averted = u5_dalys_baseline - u5_dalys,
           cases_averted = cases_baseline - cases,
           u5_cases_averted = u5_cases_baseline - u5_cases,
           deaths_averted = deaths_baseline - deaths, 
           severe_averted = severe_baseline - sev_cases, 
           u5_severe_averted = u5_severe_baseline - u5_severe) |>
    
    # calculate outcomes averted per 1000 fully vaccinated people
    mutate(dalys_avertedper1000vax = ifelse(dose3 != 0, dalys_averted/dose3 *1000, 0),
           u5_dalys_avertedper1000vax = ifelse(dose3 != 0, u5_dalys_averted/dose3*1000, 0),
           cases_avertedper1000vax = ifelse(dose3 != 0, cases_averted/dose3*1000, 0),
           u5_cases_avertedper1000vax = ifelse(dose3 != 0, u5_cases_averted/dose3*1000, 0),
           deaths_avertedper1000vax = ifelse(dose3 != 0, deaths_averted/dose3*1000, 0),
           severe_avertedper1000vax = ifelse(dose3 != 0, severe_averted/dose3*1000, 0),
           u5_severe_avertedper1000vax = ifelse(dose3 != 0, u5_severe_averted/dose3*1000, 0)) |>
    # remove NaNs for when there were 0 doses given 
    mutate(dalys_avertedper1000vax = ifelse(is.nan(dalys_avertedper1000vax) | is.infinite(dalys_avertedper1000vax), NA, dalys_avertedper1000vax),
           u5_dalys_avertedper1000vax = ifelse(is.nan(u5_dalys_avertedper1000vax) | is.infinite(u5_dalys_avertedper1000vax), NA, u5_dalys_avertedper1000vax),
           cases_avertedper1000vax = ifelse(is.nan(cases_avertedper1000vax) | is.infinite(cases_avertedper1000vax), NA, cases_avertedper1000vax),
           u5_cases_avertedper1000vax = ifelse(is.nan(u5_cases_avertedper1000vax) | is.infinite(u5_cases_avertedper1000vax), NA, u5_cases_avertedper1000vax),
           deaths_avertedper1000vax = ifelse(is.nan(deaths_avertedper1000vax) | is.infinite(deaths_avertedper1000vax), NA, deaths_avertedper1000vax),
           severe_avertedper1000vax = ifelse(is.nan(severe_avertedper1000vax) | is.infinite(severe_avertedper1000vax), NA, severe_avertedper1000vax),
           u5_severe_avertedper1000vax = ifelse(is.nan(u5_severe_avertedper1000vax) | is.infinite(u5_severe_avertedper1000vax), NA, u5_severe_avertedper1000vax))
  
  return(averted)
}
