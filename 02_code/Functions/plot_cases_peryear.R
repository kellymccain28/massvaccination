# plot cases per population over time 

plot_cases <- function(HPCfolder, strategy, seas, catchupage = '5-15', total = FALSE, compare = ''){
  #' strategy can be Catch-up or Routine or mass
  #' HPC folder is where the summarized data is saved and also indicates which vaccine
  #' seas is the seasonality - seasonal or perennial
  #' catchupage can be either 5-15 or 5-9
  #' compare can be '' which is compared to no int,'MDA','catchup12','AB'
  
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))
  df2 <- df %>% select(c(PEV, PEVstrategy, seasonality, starts_with('dose'), PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, PEVcov,MDA,
                  age_lower, age_upper, t, cases, cases_lower, cases_upper)) %>% distinct()
    
  df_plot <- df2 %>%
    # filter(PEV !='R21') %>%
    mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                                 ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                        ifelse(PEVstrategy == 'mass', 'Mass',
                                               ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
    # Filter to strategy type
    filter(strategytype == strategy) %>%
    # filter(pfpr == 0.03 | pfpr == 0.45) %>%
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
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; annual boosters',
                              label_int == 'mass 5-100 3yrs - - 4 annual, no drop' ~ 'Mass campaign every 3 years; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 3yrs - - annual no drop' ~ 'Mass campaign every 3 years; \nannual boosters no drop',
                              label_int == 'mass 5-100 single - - 4 annual, no drop' ~ 'Single mass campaign; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters no drop'),
           EPIextra_labels = case_when(EPIextra == '-'~ '12m only', 
                                       EPIextra == '5y' ~ '12m + 5y',
                                       EPIextra == '10y'~ '12m + 10y'),
           EPIextra_labels = factor(EPIextra_labels, levels = c('12m only', '12m + 5y', '12m + 10y'))) |>
    mutate(age_grp = paste0(age_lower, '-', age_upper),
           age_cuts = ifelse(age_grp %in% c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'), '1y split',
                             ifelse(age_grp %in% '0-100', 'all', 
                                    ifelse(age_grp %in% '5-15', 'school-aged', '5y split'))),
           age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                                '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100', 
                                                '0-100')),
           t = t - 5) |>
    # cases per 1000 (population was 100,000)
    mutate(casesper1000 = cases / 100,
           casesper1000_lower = cases_lower /100,
           casesper1000_upper = cases_upper / 100) %>% distinct()
  
  doses <- df_plot %>%
    pivot_longer(starts_with('dose'), 
                 names_to = "dose", 
                 values_to = "count") |>
    mutate(dose = str_replace(dose, '_median','')) |>
    mutate(dosetype = ifelse(dose =='dose1', 'Dose 1', 
                             ifelse(dose == 'dose2',"Dose 2", 
                                    ifelse(dose == 'dose3', 'Dose 3',
                                           ifelse(dose == 'dose4', 'First booster', 'Extra boosters')))),
           dosetype = ordered(dosetype, levels = c('Dose 1','Dose 2','Dose 3','First booster','Extra boosters'))) |>
    filter(age_cuts == '1y split') |>
    group_by(t, age_grp, dosetype, label_int, seasonality, pfpr, 
             PEV, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, labels) |>
    summarize(dosecount = mean(count)) %>% distinct()
              # prevalence_2_10 = mean(prevalence_2_10_median),
              # prevalence_0_100 = mean(prevalence_0_100_median)
     
  
  df_plot2 <- left_join(df_plot, doses, by = c('t', 'age_grp','label_int','seasonality', 'pfpr', 'PEV', 'PEVstrategy', 'PEVage', 
                                              'PEVrounds', 'EPIbooster', 'EPIextra', 'massbooster_rep', 'MDA'))
  
  if(strategy == 'Catch-up' & total == FALSE){
    df_plot3 <- df_plot2 %>% filter(PEVage == catchupage) 
  } 
  
  # Set colors 
  colors <- brewer.pal(n = 8, name = 'Set1')
  
  # New facet label names for pfpr variable
  pfpr.labs <- c("5%", '45%')
  names(pfpr.labs) <- c("0.05", "0.45")
  
  if(total == FALSE){
    titletext = paste0("Cases per 1000 people per year by age group, \n", df_plot$PEV[1], ', ', strategy, ' to children aged ', catchupage, ' years, ', seas)
    tot <- 'ages'
    plotdf <- df_plot3 %>%
      filter(t >= 0 & 
               age_cuts == '1y split' &
               age_upper < 21 & PEVage == catchupage & pfpr == 0.25)
    plt <- ggplot(plotdf) + #& pfpr %in% c(0.05, 0.45) 
      geom_col(aes(x = t, y = casesper1000, fill = EPIextra_labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = t, ymin = casesper1000_lower, 
                        ymax = casesper1000_upper, 
                        color = EPIextra_labels),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      # geom_col(aes(x = t, y = dosecount, fill = dosetype), alpha = 0.3) +
      facet_wrap(~age_grp, scales = 'free', nrow = 6) + 
      theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
      labs(y = 'Cases 1000 people over time, by age group',
           x = 'Year', 
           fill = 'Age-based booster timing',
           title = titletext,) +
      guides(color = 'none')
    plt
    
    ggplot(doses %>% filter(t>=0 & age_cuts=='1y split' & PEVage == catchupage & pfpr ==0.25)) +
      geom_col(aes(x = t, y = dosecount, fill = dosetype), alpha = 0.3) +
      facet_wrap(~age_grp, scales = 'free', nrow = 6) + 
      theme_bw() + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) 
      
    # Compare 5-9 and 5-15
    # plt <- ggplot(df_plot %>% filter(t > 4 & EPIextra =='-') %>% mutate(t = t-5)) + 
    #   geom_col(aes(x = t, y = .data[[paste0(compare, "cases_averted", ca_age, "_perdose")]], fill = PEVage), position ='dodge', alpha = 0.7) + #, color = '#17556d'
    #   geom_errorbar(aes(x = t, ymin = .data[[paste0(compare, "cases_averted", ca_age, "_lower_perdose")]], ymax = .data[[paste0(compare, "cases_averted", ca_age, "_upper_perdose")]], 
    #                     color = PEVage),
    #                 position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
    #   facet_grid(pfpr ~ age_grp, scales = 'free') + theme_bw() + 
    #   scale_fill_manual(values = colors) +
    #   scale_color_manual(values = colors) +
    #   scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
    #   labs(y = 'Cases averted per 1000 doses over time, by age group',
    #        x = 'Time', 
    #        fill = 'Vaccination strategy',
    #        title = "Cases averted over time by age group per 1000 doses, \nR21, Catch-up, seasonal") +
    #   guides(color = 'none')
    
  } else if (total == TRUE){
    # df_plot <- df_plot %>% filter(age_grp == '0-100')
    titletext <- paste0("Cases averted over time per 1000 doses, \n", df_plot$PEV[1], ' with 5y booster, ', strategy, ', ', seas)
    tot <- 'total'
    
    plt <- ggplot(df_plot %>% filter(age_grp=='0-100' & t>0 & EPIextra =='5y' &pfpr %in% c(0.05, 0.45))) + #
      geom_col(aes(x = t, y = casesper1000, fill = PEVage), position ='dodge', alpha = 0.7) + #, color = '#17556d'
      geom_errorbar(aes(x = t, ymin = casesper1000_lower, 
                        ymax = casesper1000_upper, color = PEVage),
                    position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
      # facet_grid(df_plot$pfpr, scales = 'free') + 
      theme_bw() + 
      facet_wrap(~pfpr, scales = 'free',
                 labeller = labeller(pfpr = pfpr.labs)) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_x_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), 
                                             breaks = NULL, labels = NULL)) +
      labs(y = 'Cases averted per 1000 doses by year',
           x = 'Year',
           fill = 'Catch-up age group',
           title = titletext,
           # caption = 'Cases averted are averaged over simulation.'
      ) +
      guides(color = 'none')
    plt
  }
  
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_", strategy, "_", seas, ", ", catchupage,'_', tot,'_', compare,'0-15_byyr.png'), width = 16, height = 8, units = 'in')
  
  }
