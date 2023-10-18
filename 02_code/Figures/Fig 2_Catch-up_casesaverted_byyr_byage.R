#Fig 2 - cases averted over tiem by age group per 1000 doses (based on plot_cases_Averted_byyr.R)
source(paste0(HPCpath, '02_code/Functions/outcomes_averted.R'))
source(paste0(HPCpath, '02_code/Functions/add_labels_groups.R'))

HPCfolder <- 'Old demography/HPC_R21'
HPCfolder <- 'HPC_R21'
# df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyagelimited_draws_CU_routine.rds')) %>% 
#   outcomes_averted(byyear = TRUE)
df <- readRDS(paste0(HPCpath, HPCfolder, '/summbyyrbyage.rds'))

df_plot <- df %>%
  select(c(PEVstrategy, PEV, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, seasonality, pfpr, MDA,
           starts_with('dose'), age_lower, age_upper, n, t, int_ID,
           clinical:prop_n, sevcases, cases, prevalence_2_10, prevalence_0_100, contains('averted'))) %>%
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
                ~ .x / n * 1000, .names = "{.col}_perpop"),
         # AB
         # across(ABcases_averted, 
         #        ~ .x / totaldoses * 1000, .names = "{.col}_perdose"),
         # across(ABcases_averted, 
         #        ~ .x / dose3 * 1000, .names = "{.col}_perFVC"),
         # across(ABcases_averted, 
         #        ~ .x / n * 1000, .names = "{.col}_perpop")
         ) |>
  mutate(age_grp = factor(age_grp, levels = c('0-5','5-10','10-15','0-100')),
         t = t - 5) %>%
  filter(t >= 0 & pfpr %in% c(0.05, 0.45) & age_grp !='0-100')# %>%
  # group_by(EPIextra_labels, scen_labels, labels,age_grp, age_upper, age_lower, t, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%
  # # Get median, and 95% CrI for each of the variables
  # summarize(across(c(clinical:prop_n, starts_with('dose'), totaldoses,
  #                    sevcases, cases,
  #                    prevalence_2_10, prevalence_0_100,
  #                    contains('averted')),#, contains('baseline')
  #                  list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
  #                       median = ~quantile(.x, 0.5, na.rm = TRUE),
  #                       upper = ~quantile(.x, 0.975, na.rm = TRUE)),
  #                  .names = "{.col}_{.fn}") ) %>%
  # rename_with(.fn = \(x)sub("_median","", x)) 

# saveRDS(df_plot, paste0(HPCpath, HPCfolder, '/CU_routine_byyrbyage.rds'))
compare <- 'AB'
compare <- ''
ca_age = ''

# New facet label names for pfpr variable
pfpr.labs <- c("5%", '45%')
names(pfpr.labs) <- c("0.05", "0.45")


titletext = paste0("Cases averted over time by age group per 1000 population")#, " , \n", df_plot$PEV[1], ', ', strategy, ' to children aged ', catchupage, ' years, ', seas, ' setting')
tot <- 'ages'
variable <- '_perdose'
variable <- '_perFVC'
variable <- '_perpop'

CUcols <- c('#B03A2E','#6C3483','#1F618D','#148F77','#B9770E','#239B56','#283747', 'tan','#85929E')

anno <- crossing(label = c('5y booster','10y booster'),
                   pfpr = c(0.05),#, 0.45),
                   age_grp = c('0-5','5-10','10-15'))
anno$age_grp <- factor(anno$age_grp, levels = c('0-5','5-10','10-15'))
anno$y <- ifelse(anno$pfpr == 0.05, 135, 1600)
anno$x <- ifelse(anno$label == '5y booster', 5.5, 10.5)

plt <- ggplot(df_plot %>%
                filter(PEVage == '5-15')#EPIextra == '5y+10y')
              )+
  # this ribbon needs to be updated to be the correctly calculated upper and lower limits (once I re-run the processing function and calc it prior to summarizing by int id)
  geom_ribbon(aes(x = t, ymin = cases_averted_lower/n * 1000, ymax = cases_averted_upper / n * 1000, fill = labels), alpha = 0.2) +
  # geom_ribbon(aes(x = t, ymin = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_lower")]],
  #                 ymax = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_upper")]], color = labels), alpha = 0.2) +
  # geom_errorbar(aes(x = t, ymin = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_lower")]], 
  #                   ymax = .data[[paste0(compare, "cases_averted", ca_age,  variable,"_upper")]], 
  #                   fill = labels), alpha = 0.2,
  #               position = position_dodge(width = 0.9), width = 0.35, linewidth = 0.6
  #             ) +
  geom_line(aes(x = t, y = .data[[paste0(compare, "cases_averted", ca_age, variable)]], color = labels),
            linewidth = 1.2
            # position ='dodge', alpha = 0.7
  ) + #, color = '#17556d'
  geom_vline(aes(xintercept = 5), linetype = 2, alpha = 0.4) + 
  geom_vline(aes(xintercept = 10), linetype = 2, alpha = 0.4) +
  geom_hline(aes(yintercept = 0), linetype = 3) +
  geom_text(data = anno, aes(x = x, y = y, label = label), angle = 90) +
  facet_grid(pfpr ~ age_grp, scales = 'free',
             labeller = labeller(pfpr = pfpr.labs)) + 
  theme_bw() + 
  scale_fill_manual(values = CUcols) +
  scale_color_manual(values = CUcols) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", expression(italic(Pf)~PR[2-10]), sep = '~')), breaks = NULL, labels = NULL)) +
  labs(y = paste0('Cases averted per 1000 population over time, by age group'),
       x = 'Year', 
       fill = 'Vaccination\nstrategy',
       color = 'Vaccination\nstrategy',
       title = titletext,
       subtitle = "Compared to no vaccination") +
  # guides(color = 'none')+ 
  theme(axis.title = element_text(size = 16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/Fig2_", "R21_", "catchup_seasonal_per1000popLINES", ", ", compare,'0.03and0.45byyr.png'), width = 16, height = 8, units = 'in')


