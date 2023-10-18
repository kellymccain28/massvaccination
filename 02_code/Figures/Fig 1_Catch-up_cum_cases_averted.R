# Plot cumulative cases averted per population (based on function in plot_cases_averted.R and is first plot in ESA)
source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))
source(paste0(HPCpath, '02_code/Functions/outcomes_averted.R'))
source(paste0(HPCpath, '02_code/Functions/add_labels_groups.R'))

HPCfolder = 'HPC_R21'
# HPCfolder ='Old demography/HPC_R21'
# df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws.rds'))%>%
#   outcomes_averted() 
df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws_CU_routine.rds'))%>%
  outcomes_averted() 
df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_CU_routine.rds'))
sim_length <- df$sim_length[1]

df_plot <- df %>%
  add_labels_groups() %>%
  # Filter to strategy type
  filter(strategytype == 'Catch-up' | strategytype == 'Routine') %>%
  # Filter to be only 1 seasonality type 
  filter(seasonality == 'seasonal') %>%
  filter(PEV !='none') %>%
  mutate(across(cases_averted, 
                ~ .x / totaldoses * 1000, .names = "{.col}_perdose"),
         across(cases_averted, 
                ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
         across(cases_averted,
                ~ .x/ totaldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
         across(cases_averted,
                ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
         across(cases_averted,
                ~ .x/ 100, .names = "{.col}_perpop"), # per 1000 population 
         # AB
         across(ABcases_averted, 
                ~ .x / totaldoses * 1000, .names = "{.col}_perdose"),
         across(ABcases_averted, 
                ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
         across(ABcases_averted,
                ~ .x/ totaldoses / (sim_length-5*365)/365 * 1000, .names = "{.col}_perdoseperyr"),
         across(ABcases_averted,
                ~ .x/ (sim_length-5*365)/365, .names = "{.col}_peryr"),
         across(ABcases_averted,
                ~ .x/ 100, .names = "{.col}_perpop")) %>% # per 1000 population 
  group_by(labels,age_grp, age_cuts, age_upper, age_lower, t, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%
  # Get median, and 95% CrI for each of the variables
  summarize(across(c(clinical:prop_n, starts_with('dose'), totaldoses,
                     sevcases, cases,
                     prevalence_2_10, prevalence_0_100,
                     contains('averted'), contains('baseline')),
  list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
       median = ~quantile(.x, 0.5, na.rm = TRUE),
       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
  .names = "{.col}_{.fn}") ) %>%
  rename_with(.fn = \(x)sub("_median","", x))  %>%
  # select out the dose_upper and lower 
  select(-matches("^dose[1-9][0-9]+_upper$|^dose[1-9]+_upper$|^dose[1-9]+_lower$|^dose[1-9][0-9]+_lower$")) 
 
saveRDS(df_plot, paste0(HPCpath, HPCfolder, '/Fig1df_summarized_CU_routine.rds'))

# colors <- brewer.pal(n = 6, name = 'Set1')
col515 <- c('#78281F','#B03A2E', '#E74C3C', '#F1948A')
col59 <- c('#4A235A', '#6C3483', '#A569BD', '#D2B4DE')
col5m15y <- c('#154360', '#1F618D', '#2980B9', '#7FB3D5')
col5m3y <- c('#0E6251', '#148F77', '#1ABC9C', '#76D7C4')
col5m5y <- c('#7E5109', '#B9770E', '#F39C12', '#F8C471')
colsAB <- c('#186A3B', '#239B56', '#58D68D', '#82E0AA')
colsSVhybrid <- c('#283747', '#85929E')
colors <- c(col515, col59, col5m15y, col5m3y, col5m5y, colsAB, colsSVhybrid)


# filter out only those that have no extra boosters
dfpl <- df_plot%>% 
  filter(age_grp == '0-100' & EPIextra == '-' & PEVstrategy != 'SV' & PEVstrategy != 'hybrid') %>%
  mutate(labels = factor(labels, levels = c("5m-3y;\n12m booster", "5m-5y;\n12m booster", "5-9y;\n12m booster", "5-15y;\n12m booster", "5m-15y;\n12m booster", "Age-based;\n12m booster")))

titletext <- paste0("Cumulative cases averted per 1000 population")#, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
tot <- 'total'
compare <- 'AB'
compare <- ''
ca_age <- ''


plt <- ggplot(dfpl) + 
  geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perpop")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
  geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_lower")]], 
                    ymax = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_upper")]], color = labels),
                position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
  # facet_grid(df_plot$pfpr, scales = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = CUcols) +
  scale_color_manual(values = CUcols) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
  labs(y = 'Cumulative cases averted per 1000 population',
       x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
       fill = 'Vaccination strategy',
       title = titletext,
       subtitle = paste0('Compared to ', if(compare == 'AB'){'age-based vaccination'} else {'no vaccination'})) +
  guides(color = 'none') + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        # legend.key.size = unit(1.2, 'cm')
        ) 
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/Fig1_", "R21_perpop_catchup_seasonal", '_', 'total', compare,'.png'), width = 16, height = 9, units = 'in')


################################################################################################
# Same plot but with only scenarios w/ 12m booster ----
dfpl_noextra <- df_plot %>% filter(age_grp == '0-100') %>% filter(EPIextra == '-')
titletext <- paste0("Cumulative cases averted per 1000 population")#, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
tot <- 'total'
compare <- 'AB'
# compare <- ''
ca_age <- ''

CUcols <- c('#B03A2E','#6C3483','#1F618D','#148F77','#B9770E','#239B56','#283747', 'tan','#85929E')

plt <- ggplot(dfpl_noextra) + 
  geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perpop")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
  geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_lower")]], 
                    ymax = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_upper")]], color = labels),
                position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
  theme_bw() + 
  scale_fill_manual(values = CUcols) +
  scale_color_manual(values = CUcols) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
  labs(y = 'Cumulative cases averted per 1000 population',
       x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
       fill = 'Vaccination strategy',
       title = titletext,
       subtitle = paste0('Compared to ', if(compare == 'AB'){'age-based vaccination'} else {'no vaccination'})) +
  guides(color = 'none') + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.8, 'cm')
  ) 
plt

################################################################################################
# Same plot but with only age-based scenarios ----
dfpl_AB <- df_plot %>% filter(age_grp == '0-100') %>% filter(PEVstrategy== 'AB' )
titletext <- paste0("Cumulative cases averted per 1000 population")#, \n", df_plot$PEV[1], ', ', strategy, ', ', seas)
tot <- 'total'
# compare <- 'AB'
compare <- ''
ca_age <- ''

CUcols <- c('#B03A2E','#6C3483','#1F618D','#148F77','#B9770E','#239B56','#283747', 'tan','#85929E')

plt <- ggplot(dfpl_AB) + 
  geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(compare, "cases_averted", ca_age, "_perpop")]], fill = labels), position ='dodge', alpha = 0.7) + #, color = '#17556d'
  geom_errorbar(aes(x = as.factor(pfpr), ymin = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_lower")]], 
                    ymax = .data[[paste0(compare, "cases_averted", ca_age, "_perpop_upper")]], color = labels),
                position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.7) +
  theme_bw() + 
  scale_fill_manual(values = CUcols) +
  scale_color_manual(values = CUcols) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
  labs(y = 'Cumulative cases averted per 1000 population',
       x = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')),
       fill = 'Vaccination strategy',
       title = titletext,
       subtitle = paste0('Compared to ', if(compare == 'AB'){'age-based vaccination'} else {'no vaccination'})) +
  guides(color = 'none') + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.8, 'cm')
  ) 
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/R21_perpop_AB_seasonal", '_', 'total_', compare,'.png'), width = 16, height = 9, units = 'in')
