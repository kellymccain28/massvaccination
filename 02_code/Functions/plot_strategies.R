# Make plot to show vaccination strategies as examples 
source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))

# Read in data 
df_month <- readRDS(paste0(HPCpath, 'HPC_testRTSSR21/summbymonth.rds'))

summbymonth_long <- df_month %>%
  select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper, 
            dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper, 
            dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper, 
            dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
            dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) |>
  pivot_longer(starts_with('dose'), 
               names_to = "dose", 
               values_to = "count") |>
  mutate(dose = str_replace(dose, '_median',''),
         label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA)) |>
  filter(label_int != 'AB - - 12m boost - - 0') |> # Don't need to show both boosted and 12 month booster 
  group_by(t, dose, label_int, seasonality, pfpr, 
           PEV, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA) |>
  summarize(dosecount = sum(count),
            prevalence_2_10 = mean(prevalence_2_10_median),
            prevalence_0_100 = mean(prevalence_0_100_median)) |> distinct() |>
  mutate(dosetype = ifelse(dose =='dose1', 'Dose 1', 
                           ifelse(dose == 'dose2',"Dose 2", 
                                  ifelse(dose == 'dose3', 'Dose 3',
                                         ifelse(dose == 'dose4', 'First booster', 'Extra boosters')))),
         dosetype = ordered(dosetype, levels = c('Dose 1','Dose 2','Dose 3','First booster','Extra boosters')),
         strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                               ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                      ifelse(PEVstrategy == 'mass', 'Mass',
                                             ifelse(PEVstrategy == 'none','No vaccine', NA)))),
         labels = case_when(label_int == "AB - - 12m - - 0"  ~ 'Age-based, 12m booster',
                            label_int == "AB - - 18m - - 0" ~ 'Age-based, 18m booster',
                            label_int == 'hybrid - - seasonal - - 0' ~ 'Hybrid',
                            label_int == 'catch-up 5-15 - 12m boost - - 0' ~ 'Catch-up to 5-15y; 12m booster',
                            label_int == 'catch-up 5-15 - 12m boost 10y - 0' ~ 'Catch-up to 5-15y; 12m, 10y booster',
                            label_int == 'catch-up 5-15 - 12m boost 5y - 0' ~ 'Catch-up to 5-15y; 12m, 5y booster',
                            label_int == 'catch-up 5-9 - 12m boost - - 0' ~ 'Catch-up to 5-9y; 12m booster',
                            label_int == 'catch-up 5-9 - 12m boost 10y - 0' ~ 'Catch-up to 5-9y; 12m, 10y booster',
                            label_int == 'catch-up 5-9 - 12m boost 5y - 0' ~ 'Catch-up to 5-9y; 12m, 5y booster',
                            label_int == 'mass 5-100 3yrs - - - 0' ~ 'Mass campaign every 3 years; 1 annual booster',
                            label_int == 'mass 5-100 3yrs - - 4 annual 0' ~ 'Mass campaign every 3 years; 4 annual boosters',
                            label_int == 'mass 5-100 3yrs - - annual 0' ~ 'Mass campaign every 3 years; annual boosters',
                            label_int == 'mass 5-100 single - - - 0' ~ 'Single mass campaign; 1 annual booster',
                            label_int == 'mass 5-100 single - - 4 annual 0' ~ 'Single mass campaign; 4 annual boosters',
                            label_int == 'mass 5-100 single - - annual 0' ~ 'Single mass campaign; annual boosters',
                            )) %>%
  filter(t > 48 & t < 190 & seasonality == 'seasonal' & PEV == 'RTSS' & MDA != 1) |>
  mutate(t = t - 12*4) # remove the part where there is no vaccination


strategies <- c('Routine','Catch-up','Mass')

plot_doses <- function(df, strategy){
  if(strategy == 'Routine'){
    mult <- 3000
    df <- df %>% filter(t < 72)
  } else if (strategy == 'Catch-up'){
    mult <- 20000
  } else {
    mult <- 80000
  }
  
  primaryseries <- c("#6BAED6","#3182BD", "#08519C")
  boosters <- c("#31A354", "#006D2C")
  
  p <- ggplot(df %>% filter(strategytype == strategy)) + 
    # geom_ribbon(aes(x = t/12, ymax = prevalence_2_10*mult, ymin = 0, color = as.character(pfpr)), fill = "#BDD7E7", alpha = 0.5) +#, fill = "#BDD7E7"
    geom_col(aes(x = t/12, y = dosecount, fill = fct_relevel(dosetype, 'Dose 1','Dose 2','Dose 3','First booster','Extra boosters')), 
                 alpha = 0.7)+ # , position = 'dodge'
    facet_wrap(~labels, scales = 'free') + 
    theme_bw() +
    labs(y = 'Dose count: ', 
         x = 'Year', 
         title = paste0(strategy, " Vaccination"),
         fill = '',
         color = '') +
    scale_fill_manual(values = c(primaryseries, boosters),
                      labels = c('Dose 1','Dose 2','Dose 3','First booster','Extra boosters')) +
    # scale_color_manual(values = "#BDD7E7",
    #                    labels = expression(italic(Pf)~PR[2-10])) + 
    scale_x_continuous(breaks = seq(0,11, by = 1))
  
  print(p)
}

for (i in strategies){
  p <- plot_doses(summbymonth_long, i)
  assign(paste0("plot", which(i == strategies)), p)
}

# plot_grid(plot1+labs(tag="A"),
#           plot2+labs(tag="B"),
#           plot3+labs(tag="C"),
#           nrow=3,align="hv",rel_heights = c(0.9,1.4, 1.4))
ggsave(filename = paste0(HPCpath, "03_output/Other figures/Fig1_Routine.png"), plot1, bg = "white",width = 15, height=8)
ggsave(filename = paste0(HPCpath, "03_output/Other figures/Fig1_Catch-up.png"), plot2, bg = "white",width = 15, height=8)
ggsave(filename = paste0(HPCpath, "03_output/Other figures/Fig1_Mass.png"), plot3, bg = "white",width = 15, height=8)

# "#EFF3FF" "#BDD7E7" "#6BAED6" "#3182BD" "#08519C"
# "#FEE5D9" "#FCAE91" "#FB6A4A" "#DE2D26" "#A50F15"
# "#FEEDDE" "#FDBE85" "#FD8D3C" "#E6550D" "#A63603"
# "#F1EEF6" "#BDC9E1" "#74A9CF" "#2B8CBE" "#045A8D"
# "#FEEBE2" "#FBB4B9" "#F768A1" "#C51B8A" "#7A0177"
# "#EDF8E9" "#BAE4B3" "#74C476" "#31A354" "#006D2C"
# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00"
# "#F2F0F7" "#CBC9E2" "#9E9AC8" "#756BB1" "#54278F"

# ggplot(summbymonth_long) + 
#   geom_line(aes(x = t, y = prevalence_2_10), color = '#17556d') + 
#   facet_wrap(~label_int + pfpr + seasonality) + theme_bw() + 
#   # scale_x_continuous(breaks = seq(min(summbyyr$t),max(summbyyr$t), by = 1)) +
#   labs(y = 'PfPR 2-10', 
#        x = 'Month', 
#        title = paste0('Initial PfPR = ')) 

ggplot(summbymonth_long %>% filter(grepl('RTSS', int_ID))) + 
  geom_col(aes(x = t, y = dosecount, fill = dose) )+ 
  facet_wrap(~int_ID, scales = 'free') + theme_bw() + 
  # scale_x_continuous(breaks = seq(min(summbyyr$t),max(summbyyr$t), by = 1)) +
  labs(y = 'Dose count: RTSS', x = 'Month', title = paste0('Initial PfPR = ', summbyyr$pfpr, ", ", summbyyr$seasonality)) 
plt