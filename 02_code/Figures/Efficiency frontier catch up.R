# Efficiency frontier
source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))
source(paste0(HPCpath, '02_code/Functions/outcomes_averted.R'))
source(paste0(HPCpath, '02_code/Functions/add_labels_groups.R'))

# scatter plot of cases averted and number of doses given 
HPCfolder <- 'Old demography/HPC_R21'
HPCfolder <- 'HPC_R21'

df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws_CU_routine.rds'))%>%
  outcomes_averted()
saveRDS(df, paste0(HPCpath, HPCfolder, '/summarized_CU_routine.rds'))

df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_CU_routine.rds'))

df_plot <- df %>%
  add_labels_groups() %>%
  filter(age_cuts == 'all') %>%
  # Filter to strategy type
  filter(strategytype == 'Catch-up' | strategytype == 'Routine') %>%
  group_by(int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA) %>%
  # mutate_at(vars(starts_with('dose'), n), median, na.rm = TRUE)%>%
  # mutate_at(vars(sevcases, cases, cases_averted, severe_averted), sum, na.rm = TRUE) %>%
  summarize(across(c(n, contains('dose'), 
                     cases, sevcases, cases_averted, severe_averted,
                     ABcases_averted, ABsevere_averted),
                   list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        median = ~quantile(.x, 0.5, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}") ) %>%
  # rename those with _median to be just the variable name 
  rename_with(.fn = \(x)sub("_median","", x)) %>%
  select(c(n, int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, labels,
           contains('dose'), starts_with('sevcases'), starts_with('cases'), starts_with('cases_averted'), starts_with('severe_averted'),
           -matches("^dose[1-9][0-9]+_upper$|^dose[1-9]+_upper$|^dose[1-9]+_lower$|^dose[1-9][0-9]+_lower$"))) %>%
  # distinct() %>%
  # Filter to be only 1 seasonality type 
  # filter(seasonality == 'seasonal') %>%
  group_by(pfpr, seasonality) %>%
  arrange(totaldoses, .by_group = TRUE) %>% 
  mutate(mincases = cummin(cases),
         minsev = cummin(sevcases),
         maxCA = cummax(cases_averted),
         maxSA = cummax(severe_averted)) %>% ungroup() %>%
  mutate(mincases = ifelse(cases == mincases, 1, 0),
         minsev = ifelse(sevcases == minsev, 1, 0),
         maxCA = ifelse(cases_averted == maxCA,1, 0),
         maxSA = ifelse(severe_averted == maxSA, 1, 0))
  
saveRDS(df_plot, paste0(HPCpath, HPCfolder, '/Efficiencyfrontier_df_summarized_CU_routine.rds'))

df_plot <- df_plot %>% filter(!(PEVstrategy %in% c('hybrid', 'SV'))) 

col515 <- c('#0E6251', '#148F77', '#1ABC9C', '#76D7C4')
col59 <- c('#4A235A', '#6C3483', '#A569BD', '#D2B4DE')
col5m15y <- c('#154360', '#1F618D', '#2980B9', '#7FB3D5')
col5m3y <- c('#7E5109', '#B9770E', '#F39C12', '#F8C471')
col5m5y <- c('#78281F','#B03A2E', '#E74C3C', '#F1948A')
colsAB <- c('#186A3B', '#239B56', '#58D68D', '#82E0AA')
colsSVhybrid <- c('#283747', '#85929E')
colors <- c(col515, col59, col5m15y, col5m3y, col5m5y, colsAB, colsSVhybrid)

# Cases and severe cases averted 
ggplot(df_plot %>% filter(maxCA == 1)) +
  geom_line(aes(x = totaldoses, y = cases_averted/n*1000), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = cases_averted/n*1000, color = labels, shape = PEVstrategy), size = 4) +
  geom_point(data = df_plot %>% filter(maxCA == 0), aes(x = totaldoses, y = cases_averted/n*1000, color = labels, shape = PEVstrategy), alpha = 0, size = 2) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 15, 18)) +
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw() +
  labs(x = 'Total doses',
       y = 'Clinical cases averted per 1000 people',
       color = 'Vaccination strategy',
       shape = 'Strategy type')
ggsave(paste0(HPCpath, '/03_output/HPC_R21/casesAVERTEDbytotaldoses_clean.png'), width = 12, height = 6)

ggplot(df_plot %>% filter(maxSA == 1)) +
  geom_line(aes(x = totaldoses, y = severe_averted/n*1000), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = severe_averted/n*1000, shape = PEVstrategy,color = labels), size = 4) +
  geom_point(data =df_plot %>% filter(maxSA == 0), aes(x = totaldoses, y = severe_averted/n*1000, color = labels, shape = PEVstrategy), alpha = 0,size = 2) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 15, 18)) +
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw() + 
  labs(x = 'Total doses',
       y = 'Severe cases averted per 1000 people',
       color = 'Vaccination strategy',
       shape = 'Strategy type')
ggsave(paste0(HPCpath, '/03_output/HPC_R21/sevcasesAVERTEDbytotaldoses_clean.png'), width = 12, height = 6)

# Cases and severe cases
ggplot(df_plot %>% filter(mincases == 1)) +
  geom_line(aes(x = totaldoses, y = cases/n), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = cases/n, color = labels, shape = PEVstrategy), size = 4) +
  geom_point(data =df_plot %>% filter(mincases == 0), aes(x = totaldoses, y = cases/n, color = labels, shape = PEVstrategy), alpha = 0, size = 2) +
  scale_x_continuous(limits = c(200000, 370000)) + 
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 15, 18)) +
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw() + 
  labs(x = 'Total doses',
       y = 'Clinical cases per person',
       color = 'Vaccination strategy',
       shape = 'Strategy type')
ggsave(paste0(HPCpath, '/03_output/HPC_R21/casesbytotaldoses_clean.png'), width = 12, height = 6)

ggplot(df_plot %>% filter(minsev == 1)) +
  geom_line(aes(x = totaldoses, y = sevcases/n), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = sevcases/n, color = labels, shape = PEVstrategy), size = 4) +
  geom_point(data =df_plot %>% filter(minsev == 0), aes(x = totaldoses, y = sevcases/n, color = labels, shape = PEVstrategy), alpha = 0, size = 2) +
  scale_x_continuous(limits = c(200000, 370000)) + 
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(17, 16, 15, 18)) +
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw() + 
  labs(x = 'Total doses',
       y = 'Severe cases per person',
       color = 'Vaccination strategy',
       shape = 'Strategy type')  
ggsave(paste0(HPCpath, '/03_output/HPC_R21/sevcasesbytotaldoses_clean.png'), width = 12, height = 6)

