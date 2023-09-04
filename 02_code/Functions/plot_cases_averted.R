#' Plot cases averted  per age group
#' this will be for catch-up and routine strategies
#' it will use the data aggregated over the entire simulation 

plot_cases_averted <- function(HPCfolder, strategy, seas, catchupage = '5-15', total = FALSE, compare = '', variable = '_perdose'){
  #' strategy can be Catch-up or Routine or mass
  #' HPC folder is where the summarized data is saved and also indicates which vaccine
  #' seas is the seasonality - seasonal or perennial
  #' catchupage can be either 5-15 or 5-9
  #' compare can be '' which is compared to no int,'MDA','catchup12','AB'
  
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized.rds'))
  if(HPCfolder == 'HPC_RTSS'){
    df <- df %>% filter(PEV != 'R21')
    sim_length <- 9125
  } else if (HPCfolder == 'HPC_R21'){
    sim_length <- 7665
  }
  # df <- readRDS('M:/Kelly/massvacc/HPC_RTSS_wrongorder/Archive/Original summarized dfs/summarized.rds')
  
  
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
    # Keep only the median number of doses 
    select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper,
              dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper,
              dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper,
              dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
              dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
    mutate(alldoses = sum(across(starts_with("dose")), na.rm = TRUE)) %>%
    mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
           labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based, 12m booster',
                              label_int == "AB - - 12m boost - -"  ~ 'Age-based, 12m booster with higher titre',
                              label_int == "AB - - 18m - -" ~ 'Age-based, 18m booster',
                              label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                              label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y;\n12m booster',
                              label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y;\n12m, 10y booster',
                              label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y;\n12m,5y booster',
                              label_int == 'catch-up 5-9 - 12m boost - -' | label_int == 'catch-up 5-9 - 12m - -' ~ 'Catch-up to 5-9y;\n12m booster',
                              label_int == 'catch-up 5-9 - 12m boost 10y -' | label_int == 'catch-up 5-9 - 12m 10y -' ~ 'Catch-up to 5-9y;\n12m, 10y booster',
                              label_int == 'catch-up 5-9 - 12m boost 5y -' | label_int == 'catch-up 5-9 - 12m 5y -' ~ 'Catch-up to 5-9y;\n12m, 5y booster',
                              label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; 1 annual booster',
                              label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; 4 annual boosters',
                              label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; annual boosters',
                              label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign; 1 annual booster',
                              label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; 4 annual boosters',
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; annual boosters',
                              label_int == 'mass 5-100 3yrs - - 4 annual, no drop' ~ 'Mass campaign every 3 years; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 3yrs - - annual no drop' ~ 'Mass campaign every 3 years; \nannual boosters no drop',
                              label_int == 'mass 5-100 single - - 4 annual, no drop' ~ 'Single mass campaign; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters no drop')) %>%
    mutate(age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                         '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                         '5-15', '5-100','0-100'))) %>%
    mutate(across(c(cases_averted_lower, cases_averted, cases_averted_upper), 
                     ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper), 
                  ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper),
                  ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper),
                  ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper),
                  ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population 
           # AB
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper), 
                  ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper),
                  ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper),
                  ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper),
                  ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population 
           # Catch-up 5-9
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper), 
                  ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper),
                  ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper),
                  ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper),
                  ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population 
           # Catch-up 5-15
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper), 
                  ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper),
                  ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper),
                  ~ .x/ 16, .names = "{.col}_peryr"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper),
                  ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population
           # MDA
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper), 
                  ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper),
                  ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper),
                  ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper),
                  ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population 
           ) 
  
  if(strategy == 'Catch-up'){
    df_plot <- df_plot %>% filter(PEVage == catchupage) 
  } 
  
  if(compare == 'catchup12' & catchupage == '5-9'){
    ca_age = '59'
  } else if(compare == 'catchup12' & catchupage == '5-15'){
    ca_age = '515'
  } else { ca_age = ''}
  
  # Set colors 
  colors <- brewer.pal(n = 6, name = 'Set1')
  
  if(variable == '_perdoseperyr'){
  if(total == FALSE){
    df_plot <- df_plot %>% filter(age_cuts == '1y split')
    titletext = paste0("Mean annual cases averted by age group per 1000 doses, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
    tot <- 'ages'
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(x = age_grp, y = .data[[paste0(compare, "cases_averted", ca_age, "_perdoseperyr")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = age_grp, ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perdoseperyr")]], ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perdoseperyr")]], color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      facet_grid(df_plot$pfpr, scales = 'free') + theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = 'Mean cases averted per 1000 doses per year',
           x = 'Age group', 
           fill = 'Vaccination strategy',
           title = titletext,
           caption = 'Cases averted are averaged over simulation.') +
      guides(color = 'none')
    
  } else if (total == TRUE){
    df_plot <- df_plot %>% filter(age_grp == '0-100')
    titletext <- paste0("Mean annual cases averted per 1000 doses, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
    tot <- 'total'
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perdoseperyr")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perdoseperyr")]], ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perdoseperyr")]], color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      # facet_grid(df_plot$pfpr, scales = 'free') + 
      theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = 'Mean cases averted per 1000 doses per year',
           x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
           fill = 'Vaccination strategy',
           title = titletext,
           caption = 'Cases averted are averaged over simulation.') +
      guides(color = 'none')
    plt
  }
  } else if (variable == '_perFVC'){
    if(total == FALSE){
      dfpl <- df_plot%>% filter(age_cuts == '1y split')
      titletext = paste0("Cumulative cases averted per 1000 FVC, by age group, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
      tot <- 'ages'
      
      plt <- ggplot(dfpl) + 
        geom_col(aes(x = age_grp, y = .data[[paste0(compare, "cases_averted", ca_age, "_perFVC")]], fill = labels), 
                 position ='dodge', alpha = 0.7) + #, color = '#17556d'
        geom_errorbar(aes(x = age_grp, ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perFVC")]], 
                          ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perFVC")]], color = labels),
                      position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
        facet_grid(dfpl$pfpr, scales = 'free') + theme_bw() + 
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
        labs(y = 'Cumulative cases averted per 1000 FVC, by age group',
             x = 'Age group', 
             fill = 'Vaccination strategy',
             title = titletext) +
        guides(color = 'none')
      plt
    } else if (total == TRUE){
      dfpl <- df_plot%>% filter(age_grp == '0-100')
      titletext <- paste0("Cumulative cases averted per 1000 FVC, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
      tot <- 'total'
      
      plt <- ggplot(dfpl) + 
        geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perFVC")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
        geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perFVC")]], 
                          ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perFVC")]], color = labels),
                      position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
        # facet_grid(df_plot$pfpr, scales = 'free') + 
        theme_bw() + 
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
        labs(y = 'Cumulative cases averted per 1000 FVC',
             x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
             fill = 'Vaccination strategy',
             title = titletext) +
        guides(color = 'none')
      plt
    }
  
  } else if (variable == '_perpop'){
    if(total == FALSE){
      dfpl <- df_plot%>% filter(age_cuts == '1y split')
      titletext = paste0("Cumulative cases averted per 1000 population, by age group, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
      tot <- 'ages'
      
      plt <- ggplot(dfpl) + 
        geom_col(aes(x = age_grp, y = .data[[paste0(compare, "cases_averted", ca_age, "_perpop")]], fill = labels), 
                 position ='dodge', alpha = 0.7) + #, color = '#17556d'
        geom_errorbar(aes(x = age_grp, ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perpop")]], 
                          ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perpop")]], color = labels),
                      position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
        facet_grid(dfpl$pfpr, scales = 'free') + theme_bw() + 
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
        labs(y = 'Cumulative cases averted per 1000 population, by age group',
             x = 'Age group', 
             fill = 'Vaccination strategy',
             title = titletext) +
        guides(color = 'none')
      plt
    } else if (total == TRUE){
      dfpl <- df_plot%>% filter(age_grp == '0-100')
      titletext <- paste0("Cumulative cases averted per 1000 population, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
      tot <- 'total'
      
      plt <- ggplot(dfpl) + 
        geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perpop")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
        geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perpop")]], 
                          ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perpop")]], color = labels),
                      position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
        # facet_grid(df_plot$pfpr, scales = 'free') + 
        theme_bw() + 
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
        labs(y = 'Cumulative cases averted per 1000 population',
             x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
             fill = 'Vaccination strategy',
             title = titletext) +
        guides(color = 'none')
      plt
    }
    } 
  plt <- plt +
    theme(axis.title = element_text(size = 16),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_", variable, strategy, "_", seas, ", ", catchupage,'_', tot,'_', compare,'.png'), width = 16, height = 9, units = 'in')
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

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', total = TRUE, variable='_perpop')
plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', total = TRUE, variable='_perFVC')

# plot_cases_averted(HPCfolder = 'HPC_testRTSSR21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = FALSE)
# plot_cases_averted(HPCfolder = 'HPC_testRTSSR21', 'Catch-up','perennial', catchupage = '5-9', total = FALSE)
# plot_cases_averted(HPCfolder = 'HPC_testRTSSR21', 'Routine', 'seasonal', catchupage = '5-9', total = FALSE)
# plot_cases_averted(HPCfolder = 'HPC_testRTSSR21', 'Routine', 'perennial', catchupage = '5-9', total = FALSE)

# Cases averted compared to AB / catchup12
plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'seasonal', total = FALSE, compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'perennial', total = FALSE, compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'seasonal', total = FALSE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'perennial', total = FALSE , compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', catchupage = '5-9', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'seasonal', catchupage = '5-9', total = FALSE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'perennial', catchupage = '5-9', total = FALSE , compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', catchupage = '5-9', total = FALSE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'seasonal', catchupage = '5-9', total = FALSE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'perennial', catchupage = '5-9', total = FALSE , compare = 'catchup12')

# Total cases averted compared to AB / catchup12
plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'seasonal', total = TRUE, compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'perennial', total = TRUE, compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'seasonal', total = TRUE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'perennial', total = TRUE , compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_RTSS', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up','perennial', catchupage = '5-9', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'seasonal', catchupage = '5-9', total = TRUE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_RTSS', 'Catch-up', 'perennial', catchupage = '5-9', total = TRUE , compare = 'catchup12')

plot_cases_averted(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up','perennial', catchupage = '5-9', total = TRUE, compare = 'AB')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'seasonal', catchupage = '5-9', total = TRUE , compare = 'catchup12')
plot_cases_averted(HPCfolder = 'HPC_R21', 'Catch-up', 'perennial', catchupage = '5-9', total = TRUE , compare = 'catchup12')

# Compare cases averted by vaccine ---------------------------------------------------

compare_cases_averted <- function(strategy, seas, catchupage = '5-15', total = FALSE, compare = '', var = ""){
  #' strategy can be Catch-up or Routine
  #' HPC folder is where the summarized data is saved and also indicates which vaccine
  #' seas is the seasonality - seasonal or perennial
  #' catchupage can be either 5-15 or 5-9
  #' compare can be '', MDA, AB, catchup12
  #' var can be raw averted (""), per dose per year ('perdoseperyr'), per year('peryr'), or per dose ('perdose') 
  
  dfR21 <- readRDS(paste0(HPCpath, 'HPC_R21/summarized.rds'))
  dfRTSS <- readRDS(paste0(HPCpath, 'HPC_RTSS/summarized.rds'))
  df <- rbind(dfR21, dfRTSS)
  
  df_plot <- df %>%
    mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                                 ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                        ifelse(PEVstrategy == 'mass', 'Mass',
                                               ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
    # Filter to strategy type
    filter(strategytype == strategy) %>%#
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
    mutate(age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                                '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                                '5-15', '5-100','0-100'))) %>%
    mutate(across(c(cases_averted_lower, cases_averted, cases_averted_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(cases_averted_lower, cases_averted, cases_averted_upper),
                  ~ .x/ 16, .names = "{.col}_peryr"),
           # AB
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(ABcases_averted_lower, ABcases_averted, ABcases_averted_upper),
                  ~ .x/ 16, .names = "{.col}_peryr"),
           # Catch-up 5-9
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(catchup12cases_averted59_lower, catchup12cases_averted59, catchup12cases_averted59_upper),
                  ~ .x/ 16, .names = "{.col}_peryr"),
           # Catch-up 5-15
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(catchup12cases_averted515_lower, catchup12cases_averted515, catchup12cases_averted515_upper),
                  ~ .x/ 16, .names = "{.col}_peryr"),
           # MDA
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper), 
                  ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper),
                  ~ .x/ alldoses / 16 * 1000, .names = "{.col}_perdoseperyr"),
           across(c(MDAcases_averted_lower, MDAcases_averted, MDAcases_averted_upper),
                  ~ .x/ 16, .names = "{.col}_peryr")) 
  
  if(strategy == 'Catch-up'){
    df_plot <- df_plot %>%
      filter(PEVage == catchupage | PEVage == '-') 
  }
  
  if(compare == 'catchup12' & catchupage == '5-9'){
    ca_age = '59'
  } else if(compare == 'catchup12' & catchupage == '5-15'){
    ca_age = '515'
  } else { ca_age = ''}
  
  # Set colors 
  colors <- brewer.pal(n = 6, name = 'Set1')
  
  if(total == FALSE){
    if(var == ''){
      labely = "Total cases averted per age group"
    } else if(var == '_perdoseperyr'){
      labely = 'Mean cases averted per 1000 doses per year'
    } else if(var == '_perdose'){
      labely = 'Mean cases averted per 1000 doses'
    } else if(var == '_peryr'){
      labely = 'Mean cases averted per year'
    }
    
    tot <- 'ages'
    df_plot <- df_plot %>% filter(age_cuts == '1y split') 
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(x = age_grp, y = .data[[paste0(compare, "cases_averted", ca_age, var)]], 
                   fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = age_grp, ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower", var)]], 
                        ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper", var)]], color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      facet_grid(df_plot$pfpr ~ df_plot$PEV, scales = 'free') + theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = labely,
           x = 'Age group', 
           fill = 'Vaccination strategy',
           title = paste0('RTSS v R21', ', ', strategy, ', ', seas, ", ", catchupage),
           caption = 'Cases averted are averaged over 16 years of intervention.') +
      guides(color = 'none')
  } else if(total == TRUE){
    if(var == ''){
      labely = "Total cases averted"
    } else if(var == 'perdoseperyr'){
      labely = 'Mean cases averted per 1000 doses per year'
    } else if(var == 'perdose'){
      labely = 'Mean cases averted per 1000 doses'
    } else if(var == 'peryr'){
      labely = 'Mean cases averted per year'
    }
    
    tot <- 'total'
    df_plot <- df_plot %>% filter(age_grp == '0-100')
    
    plt <- ggplot(df_plot) + 
      geom_col(aes(as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, var)]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower", var)]],
                        ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper", var)]], color = labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      facet_wrap(~df_plot$PEV, scales = 'free') +
      theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = labely,
           x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), 
           fill = 'Vaccination strategy',
           title = paste0('RTSS v R21', ', ', strategy, ', ', seas, ", ", catchupage),
           caption = 'Cases averted are averaged over 16 years of intervention.') +
      guides(color = 'none')
  plt
  }
  
  ggsave(paste0(HPCpath, '03_output/', "Other figures/", 'RTSSvR21', "_", strategy, "_", seas,"_", catchupage, "_", tot, "_", compare,'_', var, '.png'), width = 16, height = 9, units = 'in')
}

## perdoseperyear
# By age
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'catchup12')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'catchup12')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'catchup12')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'catchup12')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'AB')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'AB')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'AB')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'AB')

# Now by total 
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', total = TRUE)

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'catchup12', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'catchup12', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'catchup12', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'catchup12', total = TRUE)

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'AB', total = TRUE)

# For routine 
# By age
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-15')
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-9')
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-15')
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-9')

compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-15', compare = 'AB')
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-9', compare = 'AB')
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-15', compare = 'AB')
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-9', compare = 'AB')

# Now by total 
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-15', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-9', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-15', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-9', total = TRUE)

compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-15', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'seasonal', catchupage = '5-9', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-15', compare = 'AB', total = TRUE)
compare_cases_averted(strategy = 'Routine', seas = 'perennial', catchupage = '5-9', compare = 'AB', total = TRUE)

## per year ----
# By age
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', var = 'peryr')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'catchup12', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'catchup12', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'catchup12', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'catchup12', var = 'peryr')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'AB', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'AB', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'AB', var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'AB', var = 'peryr')

# Now by total 
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', total = TRUE, var = 'peryr')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'catchup12', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'catchup12', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'catchup12', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'catchup12', total = TRUE, var = 'peryr')

compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-15', compare = 'AB', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'seasonal', catchupage = '5-9', compare = 'AB', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-15', compare = 'AB', total = TRUE, var = 'peryr')
compare_cases_averted(strategy = 'Catch-up', seas = 'perennial', catchupage = '5-9', compare = 'AB', total = TRUE, var = 'peryr')