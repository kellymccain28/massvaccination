# Calculate the time to when prevalence returns to 90% of initial prevalence

#' input: processed datasets by year - paste0(HPCpath,'./HPC_summbyyr/outcomes_byyear_1_',length(index),'.rds')
#' process: takes dataset and then calculates initial pfpr and 90% threshold
#' output:
#' 

df <- readRDS(paste0(HPCpath,'HPC_summbyyr/outcomes_byyear_1_204.rds'))

d <- df %>%
  group_by(int_ID) %>%
  mutate(initial = first(prev_byyear_median),
         threshold = 0.90 * initial,
         above_threshold = ifelse(prev_byyear_median >= threshold, 1, 0),
         # get first instance of rebound within each group 
         reach_threshold = ifelse(above_threshold == 1 & lag(above_threshold == 0),1,0),
         # change NAs (first obs in each group) to 0
         reach_threshold = ifelse(is.na(reach_threshold), 0, reach_threshold),
         x_int = ifelse(reach_threshold == 1, year, NA),
         int_ID = paste0(seasonality, '_', RTSS,"_", RTSSage,"_", booster_rep,"_", MDAtiming)) 

# table showing the time to rebound 
d %>% ungroup() %>%
  filter(!is.na(x_int) & year < 11) %>%
  select(pfpr, RTSS, RTSSage, RTSSrounds, MDAtiming, time_to_threshold = x_int) %>%
  arrange(RTSS) %>%
  flextable() 
  

p0.01 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter( pfpr == 0.01)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw()+ 
  labs(title = paste0('Pfpr: 0.01'),
       ylab = 'PfPR 2-10')


p0.03 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter( pfpr == 0.03)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw()+ 
  labs(title = paste0('Pfpr: 0.03'),
       ylab = 'PfPR 2-10')

p0.05 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter(pfpr == 0.05)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw()+ 
  labs(title = paste0('Pfpr: 0.05'),
       ylab = 'PfPR 2-10')

p0.25 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter(pfpr == 0.25)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw()+ 
  labs(title = paste0('Pfpr: 0.25'),
       ylab = 'PfPR 2-10')

p0.45 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter(pfpr == 0.45)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw()+ 
  labs(title = paste0('Pfpr: 0.45'),
       ylab = 'PfPR 2-10')

p0.65 <- ggplot(d %>% filter(RTSSrounds == 'single') %>% filter(pfpr == 0.65)) +
  geom_line(aes(y = prev_byyear_median, x = year), color = '#4381C1', linewidth = 1.5) +
  geom_line(aes(y = threshold, x = year), color = 'grey40', linetype = 4) +
  geom_vline(aes(xintercept = x_int), linetype = 3, linewidth = 1.5, color = '#BEA2C2') +
  scale_x_continuous(breaks = seq(0, 21, by = 3)) +
  facet_wrap(~ int_ID) + 
  theme_bw() + 
  labs(title = paste0('Pfpr: 0.65'),
       ylab = 'PfPR 2-10') 

p0.01
p0.03
p0.05
p0.25
p0.45
p0.65
# in higher transmission areas, it never goes below 80% at all 
# in lower transmission areas, it goes below 80% and then rarely goes back above; even rarer (obviously) to go back above 90%
# rebounds more quickly and decreases less in scenarios where there is a single round of mass vaccination and it's to school-aged children