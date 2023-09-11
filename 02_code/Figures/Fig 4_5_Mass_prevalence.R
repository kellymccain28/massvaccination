# Figure 4 and 5 - prevalence over time - from plot_prevalence.R

R21df_year <- readRDS(paste0(HPCpath, 'HPC_R21/summbyyr_draws.rds')) #%>%
  # outcomes_averted(byyear = TRUE)

outcome = 'prevalence_0_100'
out_label <- expression(italic(Pf)~PR[0-100])
strategy <- 'mass'

df_plot <- R21df_year %>%
  filter(seasonality == 'seasonal') %>%
  mutate(label_int = paste(PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep),
         t = t - 5,
         labels = case_when(label_int == "AB - - 12m - -" | label_int == 'AB - - 12m boost - -' ~ 'Age-based,\n12m booster',
                            label_int == "AB - - 18m - -" ~ 'Age-based, \n18m booster',
                            label_int == 'hybrid - - seasonal - -' ~ 'Hybrid',
                            label_int == 'catch-up 5-15 - 12m boost - -' | label_int == 'catch-up 5-15 - 12m - -' ~ 'Catch-up to 5-15y; \n12m booster',
                            label_int == 'catch-up 5-15 - 12m boost 10y -' |label_int == 'catch-up 5-15 - 12m 10y -' ~ 'Catch-up to 5-15y;\n 12m, 10y booster',
                            label_int == 'catch-up 5-15 - 12m boost 5y -' |label_int == 'catch-up 5-15 - 12m 5y -' ~ 'Catch-up to 5-15y; \n12m, 5y booster',
                            label_int == 'catch-up 5-9 - 12m boost - -' |label_int == 'catch-up 5-9 - 12m - -' ~ 'Catch-up to 5-9y; \n12m booster',
                            label_int == 'catch-up 5-9 - 12m boost 10y -'|label_int == 'catch-up 5-9 - 12m 10y -' ~ 'Catch-up to 5-9y; \n12m, 10y booster',
                            label_int == 'catch-up 5-9 - 12m boost 5y -'|label_int == 'catch-up 5-9 - 12m 5y -' ~ 'Catch-up to 5-9y; \n12m, 5y booster',
                            label_int == 'mass 5-100 3yrs - - -' ~ 'Mass campaign every 3 years; \n1 booster',
                            label_int == 'mass 5-100 3yrs - - 4 annual' ~ 'Mass campaign every 3 years; \n4 annual boosters cov drop',
                            label_int == 'mass 5-100 3yrs - - annual' ~ 'Mass campaign every 3 years; \nannual boosters cov drop',
                            label_int == 'mass 5-100 single - - -' ~ 'Single mass campaign;\n 1 booster',
                            label_int == 'mass 5-100 single - - 4 annual' ~ 'Single mass campaign; \n4 annual boosters cov drop',
                            label_int == 'mass 5-100 single - - annual' ~ 'Single mass campaign; \nannual boosters cov drop',
                            label_int == 'mass 5-100 3yrs - - 4 annual, no drop' ~ 'Mass campaign every 3 years; \n4 annual boosters',
                            label_int == 'mass 5-100 3yrs - - annual no drop' ~ 'Mass campaign every 3 years; \nannual boosters',
                            label_int == 'mass 5-100 single - - 4 annual, no drop' ~ 'Single mass campaign; \n4 annual boosters',
                            label_int == 'mass 5-100 single - - annual no drop' ~ 'Single mass campaign; \nannual boosters',
                            label_int == 'none - - - - -' & MDA == 1 ~ 'MDA only',
                            label_int == 'none - - - - -' & MDA == 0 ~ 'No interventions'),
         labels = factor(labels, levels = c('Age-based,\n12m booster', 'Age-based, \n18m booster', 'Hybrid',
                                            'Catch-up to 5-15y; \n12m booster', 'Catch-up to 5-15y;\n 12m, 10y booster', 'Catch-up to 5-15y; \n12m, 5y booster',
                                            'Catch-up to 5-9y; \n12m booster', 'Catch-up to 5-9y; \n12m, 10y booster', 'Catch-up to 5-9y; \n12m, 5y booster',
                                            'Single mass campaign;\n 1 booster', 'Single mass campaign; \n4 annual boosters', 'Single mass campaign; \nannual boosters',
                                            'Mass campaign every 3 years; \n1 booster', 'Mass campaign every 3 years; \n4 annual boosters', 'Mass campaign every 3 years; \nannual boosters',
                                            'MDA only', 'No interventions',
                                            'Single mass campaign; \n4 annual boosters cov drop', 'Single mass campaign; \nannual boosters cov drop',
                                            'Mass campaign every 3 years; \n4 annual boosters cov drop', 'Mass campaign every 3 years; \nannual boosters cov drop'))) %>%
  filter(t > -2)  %>%
  filter(PEVstrategy == strategy | PEV == 'none') %>%
  arrange(PEV) %>%
  filter(grepl('no drop', massbooster_rep) | (PEVstrategy=='mass' & massbooster_rep=='-') | (PEV == 'none' & MDA == 1)) %>%
  group_by(labels, t, PEV, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA, label_int, pfpr, seasonality) %>%#age_grp, age_cuts, age_upper, age_lower,
  # Get median, and 95% CrI for each of the variables
  summarize(across(c(clinical:n, starts_with('dose'), 
                     sevcases, cases,
                     prevalence_2_10, prevalence_0_100,
                     contains('averted'), contains('baseline')),
                   list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        median = ~quantile(.x, 0.5, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}") ) %>%
  rename_with(.fn = \(x)sub("_median","", x))

# Set colors 
fillcols <- c("#BDD7E7", "#FDBE85")
linecols <- c("#6BAED6", "#FD8D3C")

# New facet label names for pfpr variable
pfpr.labs <- c("1%", "5%", '45%')
names(pfpr.labs) <- c("0.01", "0.05", "0.45")

rounds.labs <- c('Single round','Rounds every 3 years')
names(rounds.labs) <- c('single','3yrs')

massboost.labs <- c('Single booster','4 annual boosters','Annual booster')
names(massboost.labs) <- c('-','4 annual, no drop','annual no drop')

plotit <- function(massrounds){
  plt <- ggplot(df_plot %>% filter((PEV == 'R21' | PEV == 'none') & 
                                     (PEVrounds == massrounds | PEV == 'none'),
                                   pfpr %in% c(0.01, 0.05, 0.45))) + 
    geom_ribbon(aes(x = t, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]],
                    fill = as.factor(MDA)), alpha = 0.5) +
    geom_line(aes(x = t, y = .data[[outcome]],
                  color = as.factor(MDA)), linewidth = 1) + 
    theme_bw() + 
    scale_fill_manual(values = fillcols, labels = c('Vaccine only', 'MDA + vaccine')) + 
    scale_color_manual(values = linecols, labels = c('Vaccine only', 'MDA + vaccine')) +
    scale_x_continuous(breaks = c(0, seq(2, max(df_plot$t), by = 2))) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = str2expression(paste("Baseline ", 'cases', 'averted', sep = '~')), breaks = NULL, labels = NULL)) +
    labs(y = out_label,
         x = 'Year', 
         fill = '',
         color = '' ) + 
    facet_grid(pfpr ~ labels, scales = 'free',
               labeller = labeller(pfpr = pfpr.labs)) +
    theme(axis.title = element_text(size = 18),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          legend.key.size = unit(1.2, 'cm')) 
  return(plt)
}

single <- plotit('single')
ggsave(paste0(HPCpath, '03_output/', 'HPC_R21', "/Fig4_", df_plot$PEV[1], "_single", outcome, "_", strategy, "_", 'seasonal', "_R21.png"), width = 14, height = 7, units = 'in')

repeated <- plotit('3yrs')
ggsave(paste0(HPCpath, '03_output/', HPCfolder, "/Fig5_", df_plot$PEV[1], "_3yrs", outcome, "_", strategy, "_", seas, "_R21.png"), width = 14, height = 7, units = 'in')



# Get min median prevalence (Table S2)
dfmin <- df_plot %>% 
  ungroup() %>%
  mutate(massbooster_rep = str_replace(massbooster_rep, ' no drop',''),
         massbooster_rep = str_replace(massbooster_rep, ',', ''),
         massbooster_rep = str_replace(massbooster_rep, '-', '1 booster'),
         ID = interaction(PEVrounds, massbooster_rep, sep = "\n"),
         ID = factor(ID, levels = c("single\n1 booster", "single\n4 annual", "single\nannual",
                                    "3yrs\n1 booster", "3yrs\n4 annual", "3yrs\nannual")),
         MDA_lab = ifelse(MDA == 1 & PEV == 'R21', 'MDA + vaccine', 
                          ifelse(MDA ==1 & PEV =='none', 'MDA only','Vaccine only')),
         massbooster_rep = ifelse(MDA_lab == 'MDA only', '-',massbooster_rep)) %>%
  group_by(pfpr, seasonality, PEVrounds, massbooster_rep, MDA_lab) %>%
  filter((PEV == 'R21') & (PEVrounds == 'single' | PEVrounds == '3yrs') | (MDA == 1 & PEV == 'none')) %>%
  summarise(min0_100 = round(min(prevalence_0_100)*100,2) ) %>% ungroup() %>%
  select(c(pfpr, MDA_lab, everything()), -seasonality, `Minimum PfPR0-100 (%)` = min0_100) %>%
  arrange(pfpr, desc(MDA_lab), desc(PEVrounds)) 
  
print(dfmin, noSpaces = T) |> write.table("clipboard", sep = "\t", row.names = FALSE)
