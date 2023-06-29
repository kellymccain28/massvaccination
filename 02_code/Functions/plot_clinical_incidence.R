# Function to plot clinical incidence 

plot_clin_inci <- function(df, HPCfolder, agegroup, strategy, seas, vaccine){
  
  # df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))
  
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
    mutate(age_grp = paste(age_lower, age_upper, sep = '-'),,
           age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                                '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                                '5-15', '5-100','0-100'))) %>%
    filter(age_grp == agegroup) %>%
    filter(PEV == vaccine)
  
  # Set colors 
  colors <- brewer.pal(n = 6, name = 'Set1')
  
  plt <- ggplot(df_plot) + 
    geom_col(aes(x = t, y = clinical_median, fill = labels), position ='dodge', alpha = 0.8) + #, color = '#17556d'
    geom_errorbar(aes(x = t, ymin = clinical_lower, ymax = clinical_upper, group = labels), color = 'black',
                   position = position_dodge(width = 0.9), width = 0.2) +
    facet_grid(df_plot$pfpr ~ labels, scales = 'free') + theme_bw() + 
    scale_fill_manual(values = colors) +
    # scale_color_manual(values = colors) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", out_label, sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = 'Clinical incidence',
         x = 'Year', 
         fill = 'Vaccination strategy',
         title = paste0("Clinical incidence, ", agegroup, " ", df_plot$PEV[1], ', ', strategy)) #+
  plt
  
  ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/", df_plot$PEV[1], "_incidence_", agegroup,"_", strategy, "_", seas, '.png'), width = 16, height = 9, units = 'in')
}

# R21
HPCfolder <- 'HPC_R21'
df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))
d <- df %>%
  mutate(age_grp = paste(age_lower, age_upper, sep = '-'),,
         age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                              '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                              '5-15', '5-100','0-100'))) 
  
yr1 <- d %>%
  filter(age_grp %in% c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21'))
yr5 <- d %>%
  filter(age_grp %in% c('0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100'))
extra <- d %>%
  filter(age_grp %in% c('5-15', '5-100','0-100'))

saveRDS(yr1, paste0(HPCpath, HPCfolder, '/summbyyrbyage_1y.rds'))
saveRDS(yr5, paste0(HPCpath, HPCfolder, '/summbyyrbyage_5y.rds'))
saveRDS(extra, paste0(HPCpath, HPCfolder, '/summbyyrbyage_extra.rds'))

plot_clin_inci(HPCfolder = 'HPC_R21', agegroup = '0-100', strategy = 'Catch-up', seas = 'seasonal', vaccine = 'R21')
