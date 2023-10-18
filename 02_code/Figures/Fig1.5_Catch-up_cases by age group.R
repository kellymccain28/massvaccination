# Figure 1.5 - cases per person by age by scenario (catch-up only )

source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))
source(paste0(HPCpath, '02_code/Functions/outcomes_averted.R'))
HPCfolder = 'Old demography/HPC_R21'
HPCfolder = 'HPC_R21'

df <- readRDS(paste0(HPCpath, HPCfolder, '/summarized_draws_CU_routine.rds'))%>%
  outcomes_averted()

sim_length <- df$sim_length[1]

df_plot <- df %>%
  mutate(strategytype = ifelse(PEVstrategy == "AB" | PEVstrategy == 'hybrid' | PEVstrategy == 'SV', 'Routine',
                               ifelse(PEVstrategy == 'catch-up', 'Catch-up',
                                      ifelse(PEVstrategy == 'mass', 'Mass',
                                             ifelse(PEVstrategy == 'none','No vaccine', NA))))) %>%
  # Filter to strategy type
  filter(strategytype == 'Catch-up' | strategytype == 'Routine' | PEVstrategy == 'none') %>%
  # Filter to be only 1 seasonality type 
  filter(seasonality == 'seasonal') %>%
  mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
         labels = case_when(label_int == "AB - - 12m - -"  ~ 'Age-based, 12m booster',
                            label_int == "AB - - 12m boost - -"  ~ 'Age-based, 12m booster with higher titre',
                            label_int == "AB - - 18m - -" ~ 'Age-based, 18m booster',
                            label_int == "AB - - 12m 10y -" ~ "Age-based;\n12m, 10y boosters",
                            label_int == 'AB - - 12m 5y -' ~ "Age-based;\n12m, 5y boosters",
                            label_int == 'AB - - 12m 5y+10y -' ~ "Age-based;\n12m, 5y, 10y boosters",
                            label_int == 'SV - - seasonal - -' ~ 'Seasonal vaccination',
                            label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                            label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ '5-15y;\n12m booster',
                            label_int == 'catch-up 5-15 - 12m boost 10y -' | label_int == 'catch-up 5-15 - 12m 10y -' ~ '5-15y;\n12m, 10y booster',
                            label_int == 'catch-up 5-15 - 12m boost 5y -' | label_int == 'catch-up 5-15 - 12m 5y -' ~ '5-15y;\n12m, 5y booster',
                            label_int == 'catch-up 5-15 - 12m boost 5y+10y -' | label_int == 'catch-up 5-15 - 12m 5y+10y -' ~ '5-15y;\n12m, 5y, 10y boosters',
                            label_int == 'catch-up 5-9 - 12m boost - -' | label_int == 'catch-up 5-9 - 12m - -' ~ '5-9y;\n12m booster',
                            label_int == 'catch-up 5-9 - 12m boost 10y -' | label_int == 'catch-up 5-9 - 12m 10y -' ~ '5-9y;\n12m, 10y booster',
                            label_int == 'catch-up 5-9 - 12m boost 5y -' | label_int == 'catch-up 5-9 - 12m 5y -' ~ '5-9y;\n12m, 5y booster',
                            label_int == 'catch-up 5-9 - 12m boost 5y+10y -' | label_int == 'catch-up 5-9 - 12m 5y+10y -' ~ '5-9y;\n12m, 5y, 10y boosters',
                            label_int == 'catch-up 5m-15y - 12m - -' ~ '5m-15y;\n12m booster',
                            label_int == 'catch-up 5m-15y - 12m 10y -' ~ '5m-15y;\n12m, 10y boosters',
                            label_int == 'catch-up 5m-15y - 12m 5y -' ~ '5m-15y;\n12m, 5y boosters',
                            label_int == 'catch-up 5m-15y - 12m 5y+10y -' ~ '5m-15y;\n12m, 5y, 10y boosters',
                            label_int == 'catch-up 5m-3y - 12m - -' ~ '5m-3y;\n12m booster',
                            label_int == 'catch-up 5m-3y - 12m 10y -' ~ '5m-3y;\n12m, 10y boosters',
                            label_int == 'catch-up 5m-3y - 12m 5y -' ~ '5m-3y;\n12m, 5y boosters',
                            label_int == 'catch-up 5m-3y - 12m 5y+10y -' ~ '5m-3y;\n12m, 5y, 10y boosters',
                            label_int == 'catch-up 5m-5y - 12m - -' ~ '5m-5y;\n12m booster',
                            label_int == 'catch-up 5m-5y - 12m 10y -' ~ '5m-5y;\n12m, 10y boosters',
                            label_int == 'catch-up 5m-5y - 12m 5y -' ~ '5m-5y;\n12m, 5y boosters',
                            label_int == 'catch-up 5m-5y - 12m 5y+10y -' ~ '5m-5y;\n12m, 5y, 10y boosters',
                            label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; 1 annual booster',
                            label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; 4 annual boosters',
                            label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; annual boosters',
                            label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign; 1 annual booster',
                            label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; 4 annual boosters',
                            label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; annual boosters',
                            label_int == 'mass 5-100 3yrs - - 4 annual, no drop' ~ 'Mass campaign every 3 years; \n4 annual boosters no drop',
                            label_int == 'mass 5-100 3yrs - - annual no drop' ~ 'Mass campaign every 3 years; \nannual boosters no drop',
                            label_int == 'mass 5-100 single - - 4 annual, no drop' ~ 'Single mass campaign; \n4 annual boosters no drop',
                            label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters no drop',
                            label_int == 'none - - - - -' ~ 'None')) %>%
  mutate(age_grp = factor(age_grp, levels = c('0-0','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20','20-21',
                                              '0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100',
                                              '5-15', '5-100','0-100'))) %>%
  mutate(across(cases,
                ~ .x/ n, .names = "{.col}_perpop")) %>% # per person
  group_by(labels,age_grp, age_cuts, age_upper, age_lower, t, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%
  # Get median, and 95% CrI for each of the variables
  summarize(across(c(cases, cases_perpop),
                   list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        median = ~quantile(.x, 0.5, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}") ) %>%
  rename_with(.fn = \(x)sub("_median","", x)) 


col515 <- c('#78281F','#B03A2E', '#E74C3C', '#F1948A')
col59 <- c('#4A235A', '#6C3483', '#A569BD', '#D2B4DE')
col5m15y <- c('#154360', '#1F618D', '#2980B9', '#7FB3D5')
col5m3y <- c('#0E6251', '#148F77', '#1ABC9C', '#76D7C4')
col5m5y <- c('#7E5109', '#B9770E', '#F39C12', '#F8C471')
colsAB <- c('#186A3B', '#239B56', '#58D68D', '#82E0AA')
colsSVhybrid <- c('#283747', '#85929E')
colors <- c(col515, col59, col5m15y, col5m3y, col5m5y, colsAB, colsSVhybrid)


dfpl <- df_plot %>%
  filter(pfpr %in% c(0.05, 0.25, 0.45)) %>%
  # filter(PEVage == '5-15' | PEVage == '-' | PEVstrategy == 'none') %>%
  filter(age_grp %in% c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18','18-19','19-20'))
titletext <- "Cases per person, by age group"
tot <- 'total'
compare <- ''
ca_age <- ''

pfpr.labs <- c("5%", "25%", '45%')
names(pfpr.labs) <- c("0.05","0.25", "0.45")

anno <- crossing(label = c('5y booster','10y booster'),
                 pfpr = c(0.05, 0.25, 0.45),
                 age_grp = c('0-5','5-10','10-15'))
anno$age_grp <- factor(anno$age_grp, levels = c('0-5','5-10','10-15'))
anno$y <- ifelse(anno$pfpr == 0.05, 500, 2000)
anno$x <- ifelse(anno$label == '5y booster', 5.5, 10.5)

CUcols <- c('#B03A2E','#6C3483','#1F618D','#148F77','#B9770E','#239B56','#283747', 'tan','#85929E')

plt <- ggplot(dfpl %>% filter(EPIextra=='-' & (PEVage == '5m-15y' | PEVage == '5m-5y' | PEVstrategy == 'none' | PEVstrategy == 'AB'))) + 
  geom_ribbon(aes(x = age_grp, ymin = cases_perpop_lower, ymax = cases_perpop_upper, group = labels,
  fill = labels), alpha = 0.1) +
  geom_line(aes(x = age_grp, y = cases_perpop, group = labels, color = labels, linetype = labels), linewidth = 1.5, alpha = 0.8) + #, color = '#17556d'
  facet_wrap(~ pfpr, scales = 'free',
             labeller = labeller(pfpr = pfpr.labs)) +
  theme_bw() + 
  scale_color_manual(values = CUcols) +
  scale_fill_manual(values = CUcols) +
  scale_linetype_manual(values = c(2,2, 1, 1)) +
  labs(y = 'Cases per person',
       x = 'Age group (years)',
       color = 'Vaccination strategy',
       fill = 'Vaccination strategy',
       linetype = 'Vaccination strategy') +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.text.x=element_text(angle=90),
        legend.position = 'bottom') 
plt

ggsave(paste0(HPCpath, '03_output/', HPCfolder, '/Fig1.5_R21_perperson_noextraboosters_seasonal_.png'), width = 9, height = 4, units = 'in')
