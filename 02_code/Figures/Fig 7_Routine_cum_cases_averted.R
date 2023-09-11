# Fig 7 - routine cumulative cases averted 

# Plot cumulative cases averted per population (based on function in plot_cases_averted.R and is first plot in ESA)

HPCfolder = 'HPC_R21'
sim_length <- 7665


df <- readRDS(paste0(HPCpath, 'HPC_RTSS/summarized_draws.rds')) %>%
  outcomes_averted()

df_plot <- df %>%
  mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                               ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                      ifelse(PEVstrategy == 'mass', 'Mass',
                                             ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
  # Filter to strategy type
  filter(strategytype == 'Routine') %>%
  # Filter to be only 1 seasonality type 
  filter(seasonality == 'seasonal') %>%
  filter(PEV !='none') %>%
  mutate(alldoses = sum(across(starts_with("dose")), na.rm = TRUE)) %>%
  mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
         labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based,\n12m booster',
                            label_int == "AB - - 12m boost - -"  ~ 'Age-based,\n12m booster with higher titre',
                            label_int == "AB - - 18m - -" ~ 'Age-based,\n18m booster',
                            label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                            label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y;\n12m booster',
                            label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y;\n12m, 10y booster',
                            label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y;\n12m, 5y booster',
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
  mutate(across(cases_averted, 
                ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
         across(cases_averted, 
                ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
         across(cases_averted,
                ~ .x/ alldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
         across(cases_averted,
                ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
         across(cases_averted,
                ~ .x/ 100, .names = "{.col}_perpop")) %>%
  group_by(labels,age_grp, age_upper, age_lower, t, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%
  # Get median, and 95% CrI for each of the variables
  summarize(across(c(clinical:prop_n, starts_with('dose'), alldoses,
                     sevcases, cases,
                     prevalence_2_10, prevalence_0_100,
                     contains('averted'), contains('baseline')),
                   list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        median = ~quantile(.x, 0.5, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}") ) %>%
  rename_with(.fn = \(x)sub("_median","", x)) 


colors <- brewer.pal(n = 6, name = 'Set1')

dfpl <- df_plot%>% filter(age_grp == '0-100')
titletext <- paste0("Cumulative cases averted per 1000 fully vaccinated children (FVC)")#, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
tot <- 'total'
compare <- ''

plt <- ggplot(dfpl) + 
  geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perFVC")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
  geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_perFVC_lower")]], 
                    ymax = .data[[paste0(compare, "cases_averted", ca_age, "_perFVC_upper")]], color = labels),
                position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
  # facet_grid(df_plot$pfpr, scales = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = 'Cumulative cases averted per 1,000 FVC',
       x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
       fill = 'Vaccination strategy',
       title = titletext) +
  guides(color = 'none') + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.2, 'cm')) 
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/Fig8_", "RTSS_perpop_routine_seasonal", '_', 'total_', compare,'.png'), width = 16, height = 9, units = 'in')
