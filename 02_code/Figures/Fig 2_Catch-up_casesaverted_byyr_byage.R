#Fig 2 - cases averted over tiem by age group per 1000 doses (based on plot_cases_Averted_byyr.R)

df <- readRDS(paste0(HPCpath, 'HPC_R21', '/summbyyrbyagelimited_draws.rds')) %>% 
  outcomes_averted(byyear = TRUE)


df_plot <- df %>%
  mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                               ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                      ifelse(PEVstrategy == 'mass', 'Mass',
                                             ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
  # Filter to strategy type
  filter(strategytype == 'Catch-up') %>%
  # Filter to be only 1 seasonality type 
  filter(seasonality == 'seasonal') %>%
  filter(PEV !='none') %>%
  # select(-c(dose1_lower, dose1_upper, dose2_lower, dose2_upper, dose3_lower, dose3_upper, dose4_lower, dose4_upper,
  #           dose5_lower, dose5_upper, dose6_lower, dose6_upper, dose7_lower, dose7_upper, dose8_lower, dose8_upper,
  #           dose9_lower, dose9_upper, dose10_lower, dose10_upper, dose11_lower, dose11_upper, dose12_lower, dose12_upper,
  #           dose13_lower, dose13_upper, dose14_lower, dose14_upper, dose15_lower, dose15_upper, dose16_lower, dose16_upper,
  #           dose17_lower, dose17_upper, dose18_lower, dose18_upper, dose19_lower, dose19_upper)) %>%
  mutate(alldoses = sum(across(starts_with("dose")), na.rm = TRUE)) %>%
  mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
         labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based, 12m booster',
                            label_int == "AB - - 12m boost - -"  ~ 'Age-based, 12m booster with higher titre',
                            label_int == "AB - - 18m - -" ~ 'Age-based, 18m booster',
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
                            label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters no drop'),
         EPIextra_labels = case_when(EPIextra == '-'~ '12m only', 
                                     EPIextra == '5y' ~ '12m + 5y',
                                     EPIextra == '10y'~ '12m + 10y'),
         EPIextra_labels = factor(EPIextra_labels, levels = c('12m only', '12m + 5y', '12m + 10y'))) %>%
  mutate(across(cases_averted, 
                ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
         across(cases_averted, 
                ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
         # AB
         across(ABcases_averted, 
                ~ .x / alldoses * 1000, .names = "{.col}_perdose"),
         across(ABcases_averted, 
                ~ .x / dose3 * 1000, .names = "{.col}_perFVC")) |>
  mutate(age_grp = paste0(age_lower, '-', age_upper),
         age_grp = factor(age_grp, levels = c('0-5','5-10','10-15','0-100')),
         t = t - 5) %>%
  filter(PEVage == '5-15')%>%
  filter(t >= 0 & pfpr %in% c(0.05, 0.45) & age_grp !='0-100') %>%
  group_by(EPIextra_labels, labels,age_grp, age_upper, age_lower, t, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%
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

compare <- 'AB'
ca_age = ''

# New facet label names for pfpr variable
pfpr.labs <- c("5%", '45%')
names(pfpr.labs) <- c("0.05", "0.45")

titletext = paste0("Cases averted over time by age group per 1000 doses")#, " , \n", df_plot$PEV[1], ', ', strategy, ' to children aged ', catchupage, ' years, ', seas, ' setting')
tot <- 'ages'
variable <- '_perdose'
# variable <- '_perFVC'

anno <- crossing(label = c('5y booster','10y booster'),
                   pfpr = c(0.05, 0.45),
                   age_grp = c('0-5','5-10','10-15'))
anno$age_grp <- factor(anno$age_grp, levels = c('0-5','5-10','10-15'))
anno$y <- ifelse(anno$pfpr == 0.05, 130, 300)
anno$x <- ifelse(anno$label == '5y booster', 5.5, 10.5)

plt <- ggplot(df_plot %>% mutate(age_grp = factor(age_grp, levels = c('0-5','5-10','10-15','0-100')))) + 
  geom_col(aes(x = t, y = .data[[paste0(compare, "cases_averted", ca_age, variable)]], fill = EPIextra_labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
  geom_errorbar(aes(x = t, ymin = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_lower")]], 
                    ymax = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_upper")]], 
                    color = EPIextra_labels),
                position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
  geom_vline(aes(xintercept = 5), linetype = 2, alpha = 0.4) + 
  geom_vline(aes(xintercept = 10), linetype = 2, alpha = 0.4) +
  geom_text(data = anno, aes(x = x, y = y, label = label), angle = 90) +
  facet_grid(pfpr ~ age_grp, scales = 'free',
             labeller = labeller(pfpr = pfpr.labs)) + 
  theme_bw() + 
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
  labs(y = paste0('Cases averted per 1000 doses over time, by age group'),
       x = 'Year', 
       fill = 'Age-based\nbooster timing',
       title = titletext,) +
  guides(color = 'none')+ 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.2, 'cm'))
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/Fig2_", "R21_", "catchup_", seas, ", ", '515_total_', compare,'0.03and0.45byyr.png'), width = 16, height = 8, units = 'in')


