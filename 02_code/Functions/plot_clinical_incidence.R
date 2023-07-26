# Function to plot clinical incidence 

plot_clin_inci <- function(df, HPCfolder, strategy, seas, vaccine){
  
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))
  
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
                              label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; annual boosters')) 
  
  # Set colors 
  fillcols <- c("#BDD7E7", "#FDBE85")
  linecols <- c("#6BAED6", "#FD8D3C")
  out_label <- expression(italic(Pf)~PR[2-10])
  
  plt <- ggplot(df_plot) + 
    geom_ribbon(aes(x = t, ymin = clinical_lower, ymax = clinical_upper), fill = "#FDBE85",
                   position = position_dodge(width = 0.9), alpha = 0.7) +
    geom_line(aes(x = t, y = clinical_median), color = "#FD8D3C") + #, color = '#17556d'
    facet_grid(pfpr ~ labels, scales = 'free') + theme_bw() + 
    # scale_fill_manual(values = fillcols) +
    # scale_color_manual(values = linecols) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", out_label, sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = 'Clinical incidence',
         x = 'Year', 
         fill = 'Vaccination strategy',
         title = paste0("Clinical incidence, ", df_plot$PEV[1], ', ', strategy)) #+
  plt
  
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_incidence_", strategy, "_", seas, '.png'), width = 16, height = 9, units = 'in')
}

# R21
HPCfolder <- 'HPC_R21'

plot_clin_inci(HPCfolder = 'HPC_R21', strategy = 'Catch-up', seas = 'seasonal', vaccine = 'R21')
