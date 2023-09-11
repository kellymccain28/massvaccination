
# percent of cases and deaths averted 
source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))
source(paste0(HPCpath, '02_code/Functions/outcomes_averted.R'))

p_averted <- function(HPCfolder){
  df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws.rds')) %>%
    outcomes_averted() %>%
    select(c(age_lower:MDA, PEVage, age_grp,
             cases, #cases_lower, cases_upper, 
             deaths, #deaths_lower, deaths_upper, 
             cases_baseline, #cases_baseline_lower, cases_baseline_upper, 
             deaths_baseline, #deaths_baseline_lower, deaths_baseline_upper, 
             cases_averted, #cases_averted_lower, cases_averted_upper, 
             deaths_averted))# deaths_averted_lower, deaths_averted_upper))
  
  d <- df %>%
    filter(age_grp == '0-100' & PEVstrategy == 'catch-up') %>%
    mutate(p_CA = cases_averted / cases_baseline * 100,
           p_DA = deaths_averted / deaths_baseline * 100) %>% 
    group_by(int_ID, pfpr, seasonality, PEVage, EPIextra, EPIbooster) %>%
    summarize(p_CA_med = round(median(p_CA, na.rm = T),1),
              p_DA_med = round(median(p_DA, na.rm = T),1)) %>%
    distinct() %>%
    select(pfpr, PEVage, EPIextra, EPIbooster, int_ID, starts_with('p_CA'), starts_with('p_DA'), seasonality) %>%
    mutate(ID = interaction(PEVage, EPIextra, sep = "\n"),
           ID = factor(ID, levels = c('5-15\n-', '5-15\n5y', '5-15\n10y',
                                      '5-9\n-', '5-9\n5y', '5-9\n10y')))
  
  return(d)
}

r21 <- p_averted(HPCfolder = 'HPC_R21')%>% filter(seasonality == 'seasonal')
# rtss <- p_averted(HPCfolder = 'HPC_RTSS')


ggplot(r21 ) +
  geom_tile(aes(x = ID, y = as.factor(pfpr), fill = p_CA_med), alpha = 0.7) +
  scale_fill_viridis_c(direction = -1,
                       breaks = c(min(r21$p_CA_med), 20, 30, max(r21$p_CA_med)),
                       labels = c(paste0(min(r21$p_CA_med), " (min)"), 20, 30, paste0(max(r21$p_CA_med), ' (max)'))) +
  theme_bw() + 
  # facet_wrap(~seasonality) +
  labs(title = 'Percent of cases averted over 16 years, by catch-up vaccination strategy',
       x = 'Catch-up age group and fifth dose',
       y = expression(italic(Pf)~PR[0-100]),
       fill = 'Percent of cases averted') +
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.2, 'cm')) 

ggsave(paste0(HPCpath, '03_output/HPC_R21/Fig3_Perc_cases_averted.png'), width = 14, height = 7)
