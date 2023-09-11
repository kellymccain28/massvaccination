# Fig 6 - % decrease in prevalence - from Results_misc_prevalence.R

# Calculate % decrease in prevalence 

R21df_year <- readRDS(paste0(HPCpath, 'HPC_R21/summbyyr_draws.rds'))

df <- R21df_year %>%
  filter(grepl('no drop', massbooster_rep) |  massbooster_rep=='-') %>%
  filter(PEVstrategy == 'mass') %>%
  arrange(t) %>%
  group_by(identifier, pfpr, seasonality, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA) %>%
  # Find min prevalence per group and % decrease - 
  summarize(min_prev_0_100 = min(prevalence_0_100),
            min_prev_2_10 = min(prevalence_2_10),
         final1_prev_0_100 = last(prevalence_0_100),
         final2_prev_0_100 = nth(prevalence_0_100, -2),
         baseline1_prev_0_100 = first(prevalence_0_100),
         baseline2_prev_0_100 = nth(prevalence_0_100, 2)) %>%
  mutate(baseline_prev_0_100 = mean(baseline1_prev_0_100, baseline2_prev_0_100, trim = 0, na.rm=TRUE),
         final_prev_0_100 = mean(final1_prev_0_100, final2_prev_0_100, trim = 0, na.rm = TRUE),
         perc_d_toend_prev_0_100 = (baseline_prev_0_100 - final_prev_0_100) / baseline_prev_0_100* 100,
         perc_decr_prev_0_100 = (baseline_prev_0_100 - min_prev_0_100) / baseline_prev_0_100* 100,
         massbooster_rep = str_replace(massbooster_rep, ' no drop',''),
         massbooster_rep = str_replace(massbooster_rep, ',', ''),
         massbooster_rep = str_replace(massbooster_rep, '-', '1 booster'),
         ID = interaction(PEVrounds, massbooster_rep, sep = "\n"),
         ID = factor(ID, levels = c("single\n1 booster", "single\n4 annual", "single\nannual",
                                    "3yrs\n1 booster", "3yrs\n4 annual", "3yrs\nannual")),
         MDA_lab = ifelse(MDA == 1, 'MDA + vaccine', 'Vaccine only')) %>%
  group_by(ID, MDA_lab, pfpr, seasonality, PEVstrategy, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA) %>%
  summarize(baseline_prev_0_100_med = median(baseline_prev_0_100, na.rm = TRUE),
            final_prev_0_100_med = median(final_prev_0_100, na.rm = TRUE),
            perc_decr_prev_0_100_med = median(perc_decr_prev_0_100, na.rm = TRUE),
            perc_d_toend_prev_0_100_med = median(perc_d_toend_prev_0_100, na.rm = TRUE),
            min_prev_0_100_med = median(min_prev_0_100, na.rm = TRUE),
            min_prev_2_10_med = median(min_prev_2_10, na.rm = TRUE))  %>%
  filter(seasonality == 'seasonal')


ggplot(df) +
  geom_tile(aes(x = ID, y = as.factor(pfpr), fill = perc_d_toend_prev_0_100_med), alpha = 0.7) +
  facet_grid(~ MDA_lab) +
  scale_fill_viridis_c(direction = -1,
                       breaks = c(round(min(df$perc_d_toend_prev_0_100_med),0), seq(25, 100, 25)),
                       limits = c(round(min(df$perc_d_toend_prev_0_100_med),0), 100),
                       labels = c(paste0(round(min(df$perc_d_toend_prev_0_100_med),0)," (min)"), 25, 50, 75, '100 (max)')) +
  theme_bw() + 
  labs(title = 'Percent reduction of whole population parasite prevalence over 16 years,\nby mass vaccination strategy',
       x = 'Mass vaccination rounds and boosters',
       y = str2expression(paste('Baseline ', expression(italic(Pf)~PR[2-10]), sep = '~')),
       fill = str2expression(paste('Percent', 'reduction', expression(italic(Pf)~PR[0-100]), sep = '~'))) + 
  theme(axis.title = element_text(size = 18),
        plot.title = element_text(size = 22),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1, 'cm')) 

ggsave(paste0(HPCpath, '03_output/HPC_R21/Fig6_Perc_reduction_pfpr_0_100.png'), width = 14, height = 7)


# # Look at min prevalences 
# dfmin <- df %>% filter(pfpr %in% c(0.01,0.05,0.45)) %>%
#   ungroup() %>%
#   select(-c(ID, seasonality, EPIbooster, EPIextra, PEVstrategy, PEVage, MDA, perc_decr_prev_0_100_med, perc_d_toend_prev_0_100_med)) %>%
#   arrange(desc(MDA_lab)) %>%
#   mutate(min_prev_0_100_med = round(min_prev_0_100_med*100, 2),
#          min_prev_2_10_med = round(min_prev_2_10_med*100,2),
#          baseline_prev_0_100_med = round(baseline_prev_0_100_med * 100,2),
#          final_prev_0_100_med = round(final_prev_0_100_med * 100,2)) %>%
#   select(c(pfpr, everything()), -baseline_prev_0_100_med, -final_prev_0_100_med, `Minimum PfPR0-100 (%)` = min_prev_0_100_med)
# 
# print(dfmin, noSpaces = T) |> write.table("clipboard", sep = "\t", row.names = FALSE)
