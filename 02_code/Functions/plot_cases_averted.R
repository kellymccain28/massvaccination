#' Plot cases averted  per age group
#' this will be for catch-up and routine strategies
#' it will use the data aggregated over the entire simulation 

plot_cases_averted <- function(HPCfolder, strategy, seas, catchupage = '5-15', total = FALSE){
  #' strategy can be Catch-up or Routine
  #' HPC folder is where the summarized data is saved and also indicates which vaccine
  #' seas is the seasonality - seasonal or perennial
  #' catchupage can be either 5-15 or 5-9
  
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized.rds'))
  
  
  df_plot <- df %>%
    mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                                 ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                        ifelse(PEVstrategy == 'mass', 'Mass',
                                               ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
    # Filter to strategy type
    filter(strategytype == strategy) %>%
    # Filter to be only 1 seasonality type 
    filter(seasonality == seas) %>%
    filter(PEV !='none') %>%
    select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper,
              dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper,
              dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper,
              dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
              dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
    mutate(alldoses = sum(across(starts_with("dose")), na.rm = TRUE)) %>%
    mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
           labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based, 12m booster with higher titre',
                              label_int == "AB - - 12m boost - -"  ~ 'Age-based, 12m booster',
                              label_int == "AB - - 18m - -" ~ 'Age-based, 18m booster',
                              label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                              label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y; 12m booster',
                              label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y; 12m, 10y booster',
                              label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y; 12m, 5y booster',
                              label_int == 'catch-up 5-9 - 12m boost - -' | label_int == 'catch-up 5-9 - 12m - -' ~ 'Catch-up to 5-9y; 12m booster',
                              label_int == 'catch-up 5-9 - 12m boost 10y -' | label_int == 'catch-up 5-9 - 12m 10y -' ~ 'Catch-up to 5-9y; 12m, 10y booster',
                              label_int == 'catch-up 5-9 - 12m boost 5y -' | label_int == 'catch-up 5-9 - 12m 5y -' ~ 'Catch-up to 5-9y; 12m, 5y booster',
                              label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; 1 annual booster',
                              label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; 4 annual boosters',
                              label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; annual boosters',
                              label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign; 1 annual booster',
                              label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; 4 annual boosters',
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; annual boosters')) %>%
    mutate(age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                         '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                         '5-15', '5-100','0-100'))) %>%
    mutate(across(c(cases_averted_lower, cases_averted_median, cases_averted_upper), 
                     ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(cases_averted_lower, cases_averted_median, cases_averted_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(cases_averted_lower, cases_averted_median, cases_averted_upper),
                  ~ .x/ 16, .names = "{.col}_peryr")
           ) 
  
  if(strategy == 'Catch-up'){
    df_plot <- df_plot %>% filter(PEVage == catchupage) 
  }
  
  # Set colors 
  colors <- brewer.pal(n = 6, name = 'Set1')
  
  if(total == FALSE){
    df_plot <- df_plot %>% filter(age_cuts == '1y split')
    titletext = paste0("Mean annual cases averted by age group per 1000 doses, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
    tot <- 'ages'
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(x = age_grp, y = cases_averted_median_perdose, fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = age_grp, ymin = cases_averted_lower_perdose, ymax = cases_averted_upper_perdose, color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      facet_grid(df_plot$pfpr, scales = 'free') + theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = 'Mean annual cases averted per dose per year',
           x = 'Age group', 
           fill = 'Vaccination strategy',
           title = titletext,
           caption = 'Cases averted are averaged over 21-year simulation.') +
      guides(color = 'none')
    
  } else if (total == TRUE){
    df_plot <- df_plot %>% filter(age_grp == '0-100')
    titletext <- paste0("Mean annual cases averted per 1000 doses, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
    tot <- 'total'
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(x = as.factor(pfpr), y = cases_averted_median_perdose, fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = as.factor(pfpr), ymin = cases_averted_lower_perdose, ymax = cases_averted_upper_perdose, color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      # facet_grid(df_plot$pfpr, scales = 'free') + 
      theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = 'Mean annual cases averted per 1000 doses',
           x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
           fill = 'Vaccination strategy',
           title = titletext,
           caption = 'Cases averted are averaged over 21-year simulation.') +
      guides(color = 'none')
    plt
  }
  
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_", strategy, "_", seas, ", ", catchupage,'_', tot,'.png'), width = 16, height = 9, units = 'in')
}

# Plot by age group
plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'seasonal', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'perennial', total = FALSE)

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'seasonal', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'perennial', total = FALSE)

plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'seasonal', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'perennial', catchupage = '5-9', total = FALSE)

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'seasonal', catchupage = '5-9', total = FALSE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'perennial', catchupage = '5-9', total = FALSE)

# Plot total 
plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'seasonal', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'perennial', total = TRUE)

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'seasonal', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'perennial', total = TRUE)

plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'seasonal', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Routine', 'perennial', catchupage = '5-9', total = TRUE)

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'seasonal', catchupage = '5-9', total = TRUE)
plot_cases_averted(HPCfolder = 'HPC_R21', 'Routine', 'perennial', catchupage = '5-9', total = TRUE)



# Compare cases averted by vaccine ---------------------------------------------------

compare_cases_averted <- function(strategy, seas, catchupage = '5-15'){
  #' strategy can be Catch-up or Routine
  #' HPC folder is where the summarized data is saved and also indicates which vaccine
  #' seas is the seasonality - seasonal or perennial
  #' catchupage can be either 5-15 or 5-9
  
  dfR21 <- readRDS(paste0(HPCpath, 'HPC_R21/summarized.rds'))
  dfRTSS <- readRDS(paste0(HPCpath, 'HPC_RTSS/summarized.rds'))
  df <- rbind(dfR21, dfRTSS)
  
  df_plot <- df %>%
    mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                                 ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                        ifelse(PEVstrategy == 'mass', 'Mass',
                                               ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
    # Filter to strategy type
    filter(strategytype == strategy) %>%
    # Filter to be only 1 seasonality type 
    filter(seasonality == seas) %>%
    filter(PEV !='none') %>%
    select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper,
              dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper,
              dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper,
              dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
              dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
    mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
           labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based, \n12m booster with higher titre',
                              label_int == "AB - - 12m boost - -"  ~ 'Age-based, \n12m booster',
                              label_int == "AB - - 18m - -" ~ 'Age-based, \n18m booster',
                              label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                              label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y; \n12m booster',
                              label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y; \n12m, 10y booster',
                              label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y; \n12m, 5y booster',
                              label_int == 'catch-up 5-9 - 12m boost - -' | label_int == 'catch-up 5-9 - 12m - -' ~ 'Catch-up to 5-9y; \n12m booster',
                              label_int == 'catch-up 5-9 - 12m boost 10y -' | label_int == 'catch-up 5-9 - 12m 10y -' ~ 'Catch-up to 5-9y; \n12m, 10y booster',
                              label_int == 'catch-up 5-9 - 12m boost 5y -' | label_int == 'catch-up 5-9 - 12m 5y -' ~ 'Catch-up to 5-9y; \n12m, 5y booster',
                              label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; \n1 annual booster',
                              label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years;\n 4 annual boosters',
                              label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years;\n annual boosters',
                              label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign; \n1 annual booster',
                              label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; \n4 annual boosters',
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign;\n annual boosters')) %>%
    filter(age_cuts == '1y split') %>%
    mutate(age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                                '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                                '5-15', '5-100','0-100')))
  if(strategy == 'Catch-up'){
    df_plot <- df_plot %>%
      filter(PEVage == catchupage) 
  }
  
  # Set colors 
  colors <- brewer.pal(n = 4, name = 'Set1')
  
  plt <- ggplot(df_plot) + 
    geom_col(aes(x = age_grp, y = cases_averted_median, fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
    geom_errorbar(aes(x = age_grp, ymin = cases_averted_lower, ymax = cases_averted_upper, color = labels),
                  position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
    facet_grid(df_plot$pfpr ~ df_plot$PEV, scales = 'free') + theme_bw() + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = 'Cases averted',
         x = 'Age group', 
         fill = 'Vaccination strategy',
         title = paste0("Total cases averted, ", 'RTSS v R21', ', ', strategy, ', ', seas, ", ", catchupage),
         caption = 'Cases averted are summed over 21-year simulation.') +
    guides(color = 'none')
  
  plt
  
  ggsave(paste0(HPCpath, '03_output/', "Other figures/", 'RTSSvR21', "_", strategy, "_", seas,"_", catchupage,'.png'), width = 16, height = 9, units = 'in')
}

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9')
