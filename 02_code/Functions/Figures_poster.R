# Figures for poster
source(paste0(path, '02_code/Figures.R'))
output <- readRDS(paste0(HPCpath,'HPC_5355/HPC_summarized/aggregatedoversim_1_105_5y.rds')) %>%
  filter(RTSSrounds %in% c('single', 'none')) %>%
  mutate(label = paste0(RTSS, " to ", RTSSage, ', ', RTSSrounds))

plot_AVERT <- function(outcome, prev){
  if(outcome == 'dalys_averted'){
    out_label <- 'DALYs averted'
  } else if(outcome == 'cases_averted'){
    out_label <- 'Cases averted'
  } else if(outcome == 'deaths_averted'){
    out_label <- 'Deaths averted'
  } else if(outcome == 'severe_averted'){
    out_label <- 'Severe cases averted'
  } else if(outcome == 'dalys_avertedper1000vax'){
    out_label <- 'DALYs averted per 1000 \nfully vaccinated people'
  } else if(outcome == 'cases_avertedper1000vax'){
    out_label <- 'Cases averted per 1000 fully vaccinated people'
  } else if(outcome == 'deaths_avertedper1000vax'){
    out_label <- 'Deaths averted per 1000 fully vaccinated people'
  } else if(outcome == 'severe_avertedper1000vax'){
    out_label <- 'Severe cases averted per \n1000 fully vaccinated people'
  } 
  output$vaccination_label <- ifelse(output$RTSS == 'SVmass+EPI', 'EPI + seasonal mass vaccination',
                                     ifelse(output$RTSS == 'SVmass+hybrid', 'Hybrid + seasonal mass vaccination',
                                            ifelse(output$RTSS == 'mass+EPI', 'EPI + non-seasonal mass vaccination', 
                                                   ifelse(output$RTSS == 'hybrid', 'Hybrid', output$RTSS))))
  output$int_ID_lab <- paste0(output$vaccination_label, " to ", output$RTSSage)
  output$int_ID_lab <- ifelse(output$int_ID_lab == 'EPI to young children', 'EPI',
                              ifelse(output$int_ID_lab == 'Hybrid to young children', 'Hybrid', output$int_ID_lab))
  output$int_ID_lab <- ordered(output$int_ID_lab, levels = c('EPI','Hybrid','EPI + seasonal mass vaccination to everyone','EPI + seasonal mass vaccination to school-aged',
                                                             'Hybrid + seasonal mass vaccination to everyone', 'Hybrid + seasonal mass vaccination to school-aged'))
  
  lab_col <- c("#450061", "#7C00AD", "#AD0911", "#E3595F","#48610A", "#BAFA19")#, "#7CAE00", "#1D471D")
  
  plt <- ggplot(output %>% filter(pfpr == prev & seasonality == 'seasonal' & RTSS != 'mass+EPI') %>% filter(age_grp %in% c('0-5', '5-10', '10-15', '15-20', '20-25', '25-30'))) +
    geom_col(aes(x = age_grp, y = .data[[paste0(outcome, "_median")]], group = int_ID_lab, fill = int_ID_lab), 
             color = 'black', position = 'dodge') + 
    geom_linerange(aes(x = age_grp, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]], group = int_ID_lab), 
                   position = position_dodge(width = 0.9), linewidth = 0.3) +
    theme_bw() + 
    # facet_wrap(~seasonality) +
    scale_fill_manual(values = lab_col, name = 'Vaccination strategies') +
    labs(y = out_label, x = 'Age group', title = paste0('Initial PfPR = ', prev)) +
    theme(axis.text.x=element_text(angle=90))
  
  return(plt)
}

# Plot of cases averted per 1000 fully vaccinated people initial pfpr = 0.03
p1 <- plot_AVERT(outcome = 'cases_avertedper1000vax', prev = 0.03) +
  theme(legend.position = 'none')

# Plot of cases averted per 1000 fully vaccinated people initial pfpr = 0.65
p2 <- plot_AVERT(outcome = 'cases_avertedper1000vax', prev = 0.65) +
  theme(legend.position = 'none')

casesavertedplot <- p1 | p2
casesavertedplot
ggsave(paste0(path, '03_output/Figures/casesavertedplot_forposter.png'), width = 15, height = 8, units = 'in')

# Plot of deaths averted per 1000 fully vaccinated people initial pfpr = 0.03
p3 <- plot_AVERT(outcome = 'deaths_avertedper1000vax', prev = 0.03)+
  theme(legend.position = 'none')

# Plot of deaths averted per 1000 fully vaccinated people initial pfpr = 0.65
p4 <- plot_AVERT(outcome = 'deaths_avertedper1000vax', prev = 0.65)+
  theme(legend.position = c(0.55, 0.7))

deathssavertedplot <- p3 | p4
deathssavertedplot
ggsave(paste0(path, '03_output/Figures/deathsavertedplot_forposter.png'), width = 15, height = 8, units = 'in')

legend <- ggpubr::get_legend(p4)
legend
ggsave(paste0(path, '03_output/Figures/legend_forposter.png'), plot = legend, width = 15, height = 8, units = 'in')

# table  
agg <- output |>
  ungroup() |>
  filter((pfpr == 0.03 | pfpr==0.65) & seasonality == 'seasonal' & RTSS != 'mass+EPI') %>% 
  filter(age_grp %in% c('0-5', '5-10', '10-15', '15-20', '20-25', '25-30')) %>%
  select(c(int_ID, cases_avertedper1000vax_lower:cases_avertedper1000vax_upper,
           deaths_avertedper1000vax_lower:deaths_avertedper1000vax_upper,
           cases_averted_lower:cases_averted_upper,
           deaths_averted_lower:deaths_averted_upper)) |>
  group_by(int_ID) |>
  summarise(across(c(cases_avertedper1000vax_lower:deaths_avertedper1000vax_upper, 
                     cases_averted_lower:deaths_averted_upper), sum)) |>
  distinct() |>
  tidyr::separate(int_ID, into = c('PfPR','Seasonality','Vaccination strategy','Coverage','Age group vaccinated','Rounds','Fifth'), sep="_") |>
  group_by(PfPR) |>
  arrange(PfPR) |>
  mutate(across(c(cases_avertedper1000vax_lower:deaths_avertedper1000vax_upper), as.numeric)) |>
  mutate(maxcases = ifelse(cases_averted_median == max(cases_averted_median), 1, 0),
         maxdeath = ifelse(deaths_averted_median == max(deaths_averted_median), 1, 0),
         maxcases1000 = ifelse(cases_avertedper1000vax_median == max(as.numeric(cases_avertedper1000vax_median)), 1, 0),
         maxdeath1000 = ifelse(deaths_avertedper1000vax_median == max(as.numeric(deaths_avertedper1000vax_median)), 1, 0)) |>
  mutate(across(c(cases_avertedper1000vax_lower:deaths_avertedper1000vax_upper,
                  cases_averted_lower:deaths_averted_upper), ~format(round(.x, digits = 1), big.mark = ","))) |>
  mutate(`Cases averted per 1000 vaccinated` = paste0(cases_avertedper1000vax_median, " (", cases_avertedper1000vax_lower, ", ", cases_avertedper1000vax_upper, ")"),
         `Deaths averted per 1000 vaccinated` = paste0(deaths_avertedper1000vax_median, " (", deaths_avertedper1000vax_lower, ", ", deaths_avertedper1000vax_upper, ")"),
         `Cases averted` = paste0(cases_averted_median, " (", cases_averted_lower, ", ", cases_averted_upper, ")"),
         `Deaths averted` = paste0(deaths_averted_median, " (", deaths_averted_lower, ", ", deaths_averted_upper, ")"),) 
  
 

maxcases <- which(agg$maxcases == 1) 
maxdeath <- which(agg$maxdeath == 1)
maxcases1000 <- which(agg$maxcases1000 == 1) 
maxdeath1000 <- which(agg$maxdeath1000 == 1) 

border_style = officer::fp_border(color="black", width=1)

agg_out <- agg |>
  ungroup() |>
  select(PfPR:Fifth, `Cases averted per 1000 vaccinated`:`Deaths averted`, -Coverage, -Rounds, -Seasonality, -Fifth) |>
  flextable()|>
  # width(j=c(1,2,3,4,5), width = 0.7) |>
  # width(j=c(6,7,8,9), width = 2) |>
  # width(j=c(10,11,12,13), width = 2) |>
  fontsize(i = 1, size = 12, part = "header") |>  # adjust font size of header
  bold(i = 1, bold = TRUE, part = "header")  |>   # adjust bold face of header
  border_remove() |>
  theme_booktabs() |>
  # vline(j = c(5, 9), border = border_style) |>
  hline(i = c(6), border = border_style) |>
  # highlight max values for max cases averted per 1000 vaccinated
  bg(j = 4, i = maxcases1000, part = "body", bg = "#91c293") |>
  # highlight values for max deaths averted per 1000 vaccinated
  bg(j = 5, i = maxdeath1000, part = "body", bg = "#91c293") |>
  bg(j = 6, i = maxcases, part = "body", bg = "#91c293") |>
  bg(j = 7, i = maxdeath, part = "body", bg = "#91c293") 

agg_out

save_as_image(agg_out, path = paste0(path, "03_output/Figures/table_outcomes_averted_forposter.png"))
