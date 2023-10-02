# Efficiency frontier

# scatter plot of cases averted and number of doses given 
HPCfolder <- 'Old demography/HPC_R21'
df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws.rds'))%>%
  outcomes_averted()

df_plot <- df %>%
  group_by(int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA) %>%
  mutate_at(vars(starts_with('dose')), mean, na.rm = TRUE)%>%
  mutate_at(vars(sevcases, cases, cases_averted, severe_averted), sum, na.rm = TRUE) %>%
  select(c(int_ID, pfpr, seasonality, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA,
           starts_with('dose'), sevcases, cases, cases_averted, severe_averted)) %>%
  distinct() %>%
  mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid', 'Routine',
                               ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                      ifelse(PEVstrategy == 'mass', 'Mass',
                                             ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
  # Filter to strategy type
  filter(strategytype == 'Catch-up' | PEVstrategy == 'AB') %>%
  # Filter to be only 1 seasonality type 
  # filter(seasonality == 'seasonal') %>%
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
                            label_int == 'none - - - - -' ~ 'None')) %>%
  rowwise() %>%
  mutate(totaldoses = across(starts_with('dose')) %>% rowSums) 

# Cases and severe cases averted 
ggplot(df_plot) +
  geom_line(aes(x = totaldoses, y = cases_averted), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = cases_averted, color = labels), size = 3) +
  scale_x_continuous(limits = c(170000, 300000)) + 
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw()

ggplot(df_plot) +
  geom_line(aes(x = totaldoses, y = severe_averted), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = severe_averted, color = labels), size = 3) +
  scale_x_continuous(limits = c(170000, 300000)) + 
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw()

# Cases and severe cases
ggplot(df_plot) +
  geom_line(aes(x = totaldoses, y = cases), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = cases, color = labels), size = 3) +
  scale_x_continuous(limits = c(170000, 300000)) + 
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw()  
ggsave(paste0(HPCpath, '/03_output/HPC_R21/casesbytotaldoses.png'), width = 12, height = 6)

ggplot(df_plot) +
  geom_line(aes(x = totaldoses, y = sevcases), linewidth = 0.7) +
  geom_point(aes(x = totaldoses, y = sevcases, color = labels), size = 3) +
  scale_x_continuous(limits = c(170000, 300000)) + 
  facet_wrap(~pfpr + seasonality, scales = 'free') + 
  theme_bw()  
ggsave(paste0(HPCpath, '/03_output/HPC_R21/sevcasesbytotaldoses.png'), width = 12, height = 6)

