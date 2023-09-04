# Plot mass vaccination strategy alongside prevalence

plot_prev_vax <- function(seas, prev, outcome, strategy, MDAyn){
  # dfRTSS <- readRDS(paste0(HPCpath, 'HPC_RTSS/summbyyr.rds')) %>% 
  #   # removing no intervention scenarios - (will be just stochastically different from those in R21)
  #   filter(PEV !='none')
  dfR21 <- readRDS(paste0(HPCpath, 'HPC_R21/summbyyr.rds'))
  df <- dfR21%>% filter(MDA == MDAyn)
  # df <- rbind(dfR21, dfRTSS) %>% filter(MDA == MDAyn)
  
  df_plot <- df %>%
    filter(seasonality == seas) %>%
    filter(pfpr == prev) %>%
    select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper, 
              dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper, 
              dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper, 
              dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
              dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
    mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
           t = t - 5, # interventions given at 5 years
           labels = case_when(label_int == "AB - - 12m - -" | label_int == "AB - - 12m boost - -" ~ 'Age-based,\n12m booster',
                              label_int == "AB - - 18m - -" ~ 'Age-based, \n18m booster',
                              label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                              label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y; \n12m booster',
                              label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y;\n 12m, 10y booster',
                              label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y; \n12m, 5y booster',
                              label_int == 'catch-up 5-9 - 12m boost - -' | label_int == 'catch-up 5-9 - 12m - -' ~ 'Catch-up to 5-9y; \n12m booster',
                              label_int == 'catch-up 5-9 - 12m boost 10y -' | label_int == 'catch-up 5-9 - 12m 10y -' ~ 'Catch-up to 5-9y; \n12m, 10y booster',
                              label_int == 'catch-up 5-9 - 12m boost 5y -' | label_int == 'catch-up 5-9 - 12m 5y -' ~ 'Catch-up to 5-9y; \n12m, 5y booster',
                              label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; \n1 annual booster',
                              label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; \n4 annual boosters',
                              label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; \nannual boosters',
                              label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign;\n 1 annual booster',
                              label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; \n4 annual boosters',
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; \nannual boosters',
                              label_int == 'mass 5-100 3yrs - - 4 annual, no drop' ~ 'Mass campaign every 3 years; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 3yrs - - annual no drop' ~ 'Mass campaign every 3 years; \nannual boosters no drop',
                              label_int == 'mass 5-100 single - - 4 annual, no drop' ~ 'Single mass campaign; \n4 annual boosters no drop',
                              label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters no drop',
                              label_int == 'none - - - - -' & MDA == 1 ~ 'MDA only',
                              label_int == 'none - - - - -' & MDA == 0 ~ 'No interventions')) %>%
    filter(t > -2 ) 
  
  doses <- df_plot %>%
    pivot_longer(starts_with('dose'), 
                 names_to = "dose", 
                 values_to = "count") |>
    mutate(dose = str_replace(dose, '_median','')) |>
    group_by(t, dose, label_int, seasonality, pfpr, 
             PEV, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, labels) |>
    summarize(dosecount = sum(count),
              # prevalence_2_10 = mean(prevalence_2_10_median),
              # prevalence_0_100 = mean(prevalence_0_100_median)
              ) |> distinct() |>
    mutate(dosetype = ifelse(dose =='dose1', 'Dose 1', 
                             ifelse(dose == 'dose2',"Dose 2", 
                                    ifelse(dose == 'dose3', 'Dose 3',
                                           ifelse(dose == 'dose4', 'First booster', 'Extra boosters')))),
           dosetype = ordered(dosetype, levels = c('Dose 1','Dose 2','Dose 3','First booster','Extra boosters')))
  
  df_plot <- left_join(df_plot, doses, by = c('t', 'label_int','seasonality', 'pfpr', 'PEV', 'PEVstrategy', 'PEVage', 
                                              'PEVrounds', 'EPIbooster', 'EPIextra', 'massbooster_rep', 'MDA', 'labels'))
  
  # Set colors 
  fillcols <- c("#FBB4B9", "#CBC9E2", "#FCAE91", palette(rainbow(19)))
  linecols <- c("#F768A1", "#9E9AC8", "#FB6A4A")
  
  if(strategy == 'mass'){
    df_plot <- df_plot %>%
      filter(PEVstrategy == strategy | PEV == 'none') 
    doses <- doses %>%
      filter(PEVstrategy == strategy | PEV == 'none') 
  } else if (strategy == 'catch-up'){
    # Filter to scenarios to plot -- none, AB with 12m booster (boosted immunity for RTSS, normal for R21)
    df_plot <- df_plot %>%
      filter(PEVstrategy == strategy | PEV == 'none' | (PEVstrategy == 'AB' & ((PEV == 'R21' & EPIbooster == '12m') | (PEV == 'RTSS' & EPIbooster == '12m boost'))))
    doses <- doses %>%
      filter(PEVstrategy == strategy | PEV == 'none' | (PEVstrategy == 'AB' & ((PEV == 'R21' & EPIbooster == '12m') | (PEV == 'RTSS' & EPIbooster == '12m boost'))))
  }
  
  plt <- ggplot(df_plot) + 
    geom_ribbon(aes(x = t, ymin = prevalence_0_100_lower, ymax = prevalence_0_100_upper),
                    alpha = 0.6) +# fill = '#b9ccd3',color = as.factor(PEV)), 
    geom_line(aes(x = t, y = prevalence_0_100)) + #, color = '#17556d' color = as.factor(PEV))
    geom_col(aes(x = t, y = dosecount/50000, fill = dosetype), alpha = 0.8) +
    facet_wrap(~labels, scales = 'free') + theme_bw() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    scale_x_continuous(breaks = seq(min(df_plot$t), max(df_plot$t), by = 1)) +
    scale_y_continuous(sec.axis = sec_axis(~ .*50000 , name = 'Vaccine doses')) +
    labs(y = expression(italic(Pf)~PR[0-100]),
         x = 'Year', 
         title = str2expression(paste('RTSS', 'vs', 'R21', strategy, 'vaccination', sep = '~')),
         caption = paste0(seas, if(MDAyn == 1) ' with MDA' else " without MDA"),
         fill = 'Vaccine',
         color = 'Vaccine') #+
  plt
  
  # p <- ggplot(doses) +
  #   geom_col(aes(x = t, y = dosecount)) +
  #   facet_grid(~labels, scales = 'free') + theme_bw() +
  #   scale_fill_manual(values = fillcols) + 
  #   scale_color_manual(values = linecols) +
  #   scale_x_continuous(breaks = seq(min(df_plot$t), max(df_plot$t), by = 1)) +
  #   labs(y = out_label,
  #        x = 'Year', 
  #        title = str2expression(paste('RTSS', 'vs', 'R21', out_label, strategy, 'vaccination', sep = '~')),
  #        caption = paste0(seas, if(MDAyn == 1) ' with MDA' else " without MDA"),
  #        fill = 'Vaccine',
  #        color = 'Vaccine')
  # p
  # 
  # ggarrange(plt, p, nrow = 2)
  
  ggsave(paste0(HPCpath, '03_output/Other figures/', 'prev+vax_RTSSvR21', "_", outcome, "_", strategy, "_", seas, "_", MDAyn, '.png'), width = 16, height = 9, units = 'in')
}


plot_prev_vax(seas = 'seasonal', prev = 0.25, outcome = 'prevalence_0_100', strategy = 'mass', MDAyn = 0)
plot_prev_vax(seas = 'seasonal', prev = 0.25, outcome = 'prevalence_0_100', strategy = 'mass', MDAyn = 1)

plot_prev_vax(seas = 'seasonal', prev = 0.25, outcome = 'prevalence_0_100', strategy = 'catch-up', MDAyn = 0)
plot_prev_vax(seas = 'perennial', prev = 0.25, outcome = 'prevalence_0_100', strategy = 'catch-up', MDAyn = 0)
