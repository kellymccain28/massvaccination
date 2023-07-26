# Plot an annual outcome

RTSSdf_year <- readRDS(paste0(HPCpath, 'HPC_RTSS/summbyyr.rds'))
R21df_year <- readRDS(paste0(HPCpath, 'HPC_R21/summbyyr.rds'))
testdf_year <- readRDS(paste0(HPCpath, 'HPC_testRTSSR21/summbyyr.rds'))

plot_prevalence <- function(HPCfolder, outcome, strategy, seas, vaccine){
  #' outcome can be prevalence_2_10 or prevalence_0_100
  if(outcome == 'prevalence_2_10'){
    out_label <- expression(italic(Pf)~PR[2-10])
  } else if (outcome == 'prevalence_0_100'){
    out_label <- expression(italic(Pf)~PR[0-100])
  } 
  
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyr.rds'))
  
  df_plot <- df %>%
    filter(seasonality == seas) %>%
    select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper, 
              dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper, 
              dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper, 
              dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
              dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
    mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
           t = t - 5,
           labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based,\n12m booster',
                              label_int == "AB - - 18m - -" ~ 'Age-based, \n18m booster',
                              label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                              label_int == 'catch-up 5-15 - 12m boost - -' ~ 'Catch-up to 5-15y; \n12m booster',
                              label_int == 'catch-up 5-15 - 12m boost 10y -' ~ 'Catch-up to 5-15y;\n 12m, 10y booster',
                              label_int == 'catch-up 5-15 - 12m boost 5y -' ~ 'Catch-up to 5-15y; \n12m, 5y booster',
                              label_int == 'catch-up 5-9 - 12m boost - -' ~ 'Catch-up to 5-9y; \n12m booster',
                              label_int == 'catch-up 5-9 - 12m boost 10y -' ~ 'Catch-up to 5-9y; \n12m, 10y booster',
                              label_int == 'catch-up 5-9 - 12m boost 5y -' ~ 'Catch-up to 5-9y; \n12m, 5y booster',
                              label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; \n1 annual booster',
                              label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; \n4 annual boosters',
                              label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; \nannual boosters',
                              label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign;\n 1 annual booster',
                              label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; \n4 annual boosters',
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; \nannual boosters',
                              label_int == 'none - - - - -' & MDA == 1 ~ 'MDA only',
                              label_int == 'none - - - - -' & MDA == 0 ~ 'No interventions')) %>%
    filter(t > -2)  %>%
    filter(PEVstrategy == strategy | PEV == 'none') 
  
  # Set colors 
  fillcols <- c("#BDD7E7", "#FDBE85")
  linecols <- c("#6BAED6", "#FD8D3C")
  
  plt <- ggplot(df_plot %>% filter(PEV == vaccine | PEV == 'none')) + 
    geom_ribbon(aes(x = t, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]],
                    fill = as.factor(MDA)), alpha = 0.6) +# fill = '#b9ccd3',
    geom_line(aes(x = t, y = .data[[paste0(outcome, "_median")]],
                  color = as.factor(MDA))) + #, color = '#17556d'
    facet_grid(pfpr ~ labels, scales = 'free') + theme_bw() + 
    scale_fill_manual(values = fillcols) + 
    scale_color_manual(values = linecols) +
    scale_x_continuous(breaks = c(0, seq(2,max(df_plot$t), by = 2))) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", out_label, sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = out_label,
         x = 'Year', 
         title = str2expression(paste(df_plot$PEV[50], out_label, "mass", 'vaccination', sep = '~')),
         caption = seas,
         fill = 'MDA',
         color = 'MDA') #+
  plt
  
  
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_", outcome, "_", strategy, "_", seas, "_", vaccine, '.png'), width = 16, height = 9, units = 'in')
}

plot_prevalence('HPC_RTSS', 'prevalence_2_10', 'mass', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_2_10', 'mass', 'perennial', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_0_100', 'mass', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_0_100', 'mass', 'perennial', vaccine = 'RTSS')

plot_prevalence('HPC_R21', 'prevalence_2_10', 'mass', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_2_10', 'mass', 'perennial', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_0_100', 'mass', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_0_100', 'mass', 'perennial', vaccine = 'R21')

plot_prevalence('HPC_RTSS', 'prevalence_2_10', 'catch-up', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_2_10', 'catch-up', 'perennial', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_0_100', 'catch-up', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_RTSS', 'prevalence_0_100', 'catch-up', 'perennial', vaccine = 'RTSS')

plot_prevalence('HPC_R21', 'prevalence_2_10', 'catch-up', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_2_10', 'catch-up', 'perennial', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_0_100', 'catch-up', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_R21', 'prevalence_0_100', 'catch-up', 'perennial', vaccine = 'R21')

# test runs
plot_prevalence('HPC_testRTSSR21', 'prevalence_2_10', 'mass', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_testRTSSR21', 'prevalence_2_10', 'mass', 'perennial', vaccine = 'RTSS')
plot_prevalence('HPC_testRTSSR21', 'prevalence_0_100', 'mass', 'seasonal', vaccine = 'RTSS')
plot_prevalence('HPC_testRTSSR21', 'prevalence_0_100', 'mass', 'perennial', vaccine = 'RTSS')

plot_prevalence('HPC_testRTSSR21', 'prevalence_2_10', 'mass', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_testRTSSR21', 'prevalence_2_10', 'mass', 'perennial', vaccine = 'R21')
plot_prevalence('HPC_testRTSSR21', 'prevalence_0_100', 'mass', 'seasonal', vaccine = 'R21')
plot_prevalence('HPC_testRTSSR21', 'prevalence_0_100', 'mass', 'perennial', vaccine = 'R21')


# Compare R21 and RTSS ----------------------------------------------------
compare_prevalence <- function(outcome, strategy, seas, MDAyn){
  #' outcome can be prevalence_2_10 or prevalence_0_100
  if(outcome == 'prevalence_2_10'){
    out_label <- expression(italic(Pf)~PR[2-10])
  } else if (outcome == 'prevalence_0_100'){
    out_label <- expression(italic(Pf)~PR[0-100])
  } 
  
  dfRTSS <- readRDS(paste0(HPCpath, 'HPC_RTSS/summbyyr.rds')) %>% 
    # removing no intervention scenarios - (will be just stochastically different from those in R21)
    filter(PEV !='none')
  dfR21 <- readRDS(paste0(HPCpath, 'HPC_R21/summbyyr.rds'))
  
  df <- rbind(dfR21, dfRTSS) %>% filter(MDA == MDAyn)
  
  df_plot <- df %>%
    filter(seasonality == seas) %>%
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
                              label_int == 'none - - - - -' & MDA == 1 ~ 'MDA only',
                              label_int == 'none - - - - -' & MDA == 0 ~ 'No interventions')) %>%
    filter(t > -2 )
  
  # Set colors 
  fillcols <- c("#FBB4B9", "#CBC9E2", "#FCAE91")
  linecols <- c("#F768A1", "#9E9AC8", "#FB6A4A")
  
  if(strategy == 'mass'){
    df_plot <- df_plot %>%
      filter(PEVstrategy == strategy | PEV == 'none') 
  } else if (strategy == 'catch-up'){
    # Filter to scenarios to plot -- none, AB with 12m booster (boosted immunity for RTSS, normal for R21)
    df_plot <- df_plot %>%
      filter(PEVstrategy == strategy | PEV == 'none' | (PEVstrategy == 'AB' & ((PEV == 'R21' & EPIbooster == '12m') | (PEV == 'RTSS' & EPIbooster == '12m boost'))))
  }
  
  plt <- ggplot(df_plot) + 
    geom_ribbon(aes(x = t, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]],
                    fill = as.factor(PEV)), alpha = 0.6) +# fill = '#b9ccd3',
    geom_line(aes(x = t, y = .data[[paste0(outcome, "_median")]],
                  color = as.factor(PEV))) + #, color = '#17556d'
    facet_grid(pfpr ~ labels, scales = 'free') + theme_bw() +
    scale_fill_manual(values = fillcols) + 
    scale_color_manual(values = linecols) +
    scale_x_continuous(breaks = seq(min(df_plot$t),max(df_plot$t), by = 1)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", out_label, sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = out_label,
         x = 'Year', 
         title = str2expression(paste('RTSS', 'vs', 'R21', out_label, strategy, 'vaccination', sep = '~')),
         caption = paste0(seas, if(MDAyn == 1) ' with MDA' else " without MDA"),
         fill = 'Vaccine',
         color = 'Vaccine') #+
  plt
  
  ggsave(paste0(HPCpath, '03_output/Other figures/', 'RTSSvR21', "_", outcome, "_", strategy, "_", seas, "_", MDAyn, '.png'), width = 16, height = 9, units = 'in')
}

# Mass
compare_prevalence(outcome = 'prevalence_2_10', strategy = 'mass', seas = 'seasonal', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_2_10', strategy = 'mass', seas = 'seasonal', MDAyn = 1)

compare_prevalence(outcome = 'prevalence_2_10', strategy = 'mass', seas = 'perennial', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_2_10', strategy = 'mass', seas = 'perennial', MDAyn = 1)

compare_prevalence(outcome = 'prevalence_0_100', strategy = 'mass', seas = 'seasonal', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_0_100', strategy = 'mass', seas = 'seasonal', MDAyn = 1)

compare_prevalence(outcome = 'prevalence_0_100', strategy = 'mass', seas = 'perennial', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_0_100', strategy = 'mass', seas = 'perennial', MDAyn = 1)

# Catchup
compare_prevalence(outcome = 'prevalence_2_10', strategy = 'catch-up', seas = 'seasonal', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_2_10', strategy = 'catch-up', seas = 'perennial', MDAyn = 0)

compare_prevalence(outcome = 'prevalence_0_100', strategy = 'catch-up', seas = 'seasonal', MDAyn = 0)
compare_prevalence(outcome = 'prevalence_0_100', strategy = 'catch-up', seas = 'perennial', MDAyn = 0)

