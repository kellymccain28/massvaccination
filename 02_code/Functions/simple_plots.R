plot_clin <- function(){
  g <- ggplot(data = output, aes(x = timestep / year, y = clinical_incidence)) + 
    geom_line(col = "grey80") +
    stat_smooth(col = "darkblue", se = FALSE) +
    xlab("year") +
    ylab("Clinical incidence \n (per 1000 children aged 0-5)") +
    theme_bw()
  
  if (any(names(output) == "n_rtss_epi_dose_1" | any(names(output) == "n_rtss_mass_dose_1"))) {
   g + 
      geom_vline(xintercept = timesteps/year, lty = 2)
  } else {
    g
  }
}

plot_sev <- function(){
  g <- ggplot(data = output, aes(x = timestep / year, y = severe_incidence)) + 
    geom_line(col = "grey80") +
    stat_smooth(col = "darkblue", se = FALSE) +
    xlab("year") +
    ylab("Incidence of severe malaria \n (per 1000 children aged 0-5)") +
    theme_bw()
  
  if (any(names(output) == "n_rtss_epi_dose_1" | any(names(output) == "n_rtss_mass_dose_1"))) {
    g + 
      geom_vline(xintercept = timesteps/year, lty = 2)
  } else {
    g
  }
}

plot_doses <- function(){
  if(any(names(output) == 'n_rtss_epi_dose_1')) {
    doses <- output |> 
      select(timestep, n_rtss_epi_dose_1:n_rtss_epi_booster_1) |>
      pivot_longer(n_rtss_epi_dose_1:n_rtss_epi_booster_1, names_to = "dose", values_to = "count") |>
      mutate(month = ceiling(timestep / month)) |>
      group_by(month, dose) |>
      summarize(count = sum(count)) 
  } else if (any(names(output) == 'n_rtss_mass_dose_1')) {
    doses <- output |> 
      select(timestep, n_rtss_mass_dose_1:n_rtss_mass_booster_1) |>
      pivot_longer(n_rtss_mass_dose_1:n_rtss_mass_booster_1, names_to = "dose", values_to = "count") |>
      group_by(month, dose) |>
      summarize(count = sum(count))
  }
  
  # plot doses
  ggplot(data = doses) + 
    geom_col(aes(x = month, y = count, fill = dose)) +
    geom_vline(xintercept = timesteps/month, lty = 2) +  
    xlab("month") +
    ylab("doses") +
    theme_bw() +
    theme(legend.position="bottom") 
}

plot_prev <- function(){
  g <- ggplot(data = output, aes(x = timestep / year, y = prevalence)) + 
    geom_line(col = "grey80") +
    stat_smooth(col = "darkblue", se = FALSE) +
    xlab("year") +
    ylab("Prevalence \n (per 1000 children aged 0-5)") +
    theme_bw()
  
  if (any(names(output) == "n_rtss_epi_dose_1" | any(names(output) == "n_rtss_mass_dose_1"))) {
    g + 
      geom_vline(xintercept = timesteps/year, lty = 2)
  } else {
    g
  }
}

get_diff <- function(){
  out_rtss <- output %>%
    select(-c(starts_with('n_rtss'))) %>%
    mutate(timestep = 0)
  
  diff <- output_no_int - out_rtss 
  diff <- diff %>%
    select(c(timestep, starts_with('n_inc_clin'), starts_with('n_det'),starts_with('n_')), -n_inc_severe_0_1825) %>%
    mutate_at(vars(n_0_1825:n_730_3650,
                   n_detect_730_3650), mean, na.rm = TRUE) %>%
    mutate_at(vars(n_inc_clinical_0_1825:n_inc_clinical_7300_7665), sum, na.rm = TRUE) %>%
    select(-c(timestep, n_detect_730_3650)) %>%
    distinct() %>%
    pivot_longer(cols = c(n_0_1825:n_730_3650,
                          n_inc_clinical_0_1825:n_inc_clinical_7300_7665),
                 names_to = c('age'), values_to = c('value')) %>%
    mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
           inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
           age = gsub('n_inc_clinical_', '', age),                         # combining age vars
           age = gsub('n_', '', age),
           age = gsub('_', '-', age)) %>%
    group_by(age) %>%
    select(-c(value)) %>%
    mutate_at(vars(n:inc_clinical), sum, na.rm = TRUE) %>% distinct() %>%
    separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
    mutate(age_lower = as.numeric(age_lower)/365,
           age_upper = as.numeric(age_upper)/365,
           inc = inc_clinical / n,
           cases = inc_clinical)
}



######### testing
doses <- output |> 
  mutate(month = ceiling(timestep / 30)) |>
  select(month, starts_with('n_rtss')) |>
  pivot_longer(starts_with('n_rtss'), names_to = "dose", values_to = "count") |>
  group_by(month, dose) |>
  summarize(count = sum(count))

# plot doses
dosesplot <- ggplot(data = doses) + 
  geom_col(aes(x = month, y = count, fill = dose)) +
  # geom_vline(xintercept = (timesteps-warmup) / month, lty = 2) + # start of vaccination campaign
  # geom_vline(xintercept = (peak / month) + c(0,year, year*2, year*3, year*6, year*9)/month, lty = 1) + # peak transmission 
  # annotate("text", x = 25, y = 6000, label = "peak transmission") +
  # facet_wrap(~ model, nrow = 3, ncol = 2) + 
  xlab("month") +
  ylab("doses") +
  theme_bw()
dosesplot


# plot clinical incidence
clin_inc <- output %>%
  mutate(month = ceiling(timestep/30)) %>%
  # group_by(month) %>%
  # summarize(n_inc_clinical_0_1825_mo = mean(n_inc_clinical_0_1825),
  #           n_0_1825_mo_mean = mean(n_0_1825),
  #           mo_clin_inc_0_1825 = 1000 * n_inc_clinical_0_1825_mo / n_0_1825_mo_mean,
  #           # yr_clinical_incidence_0_1825 = sum(clinical_incidence_0_1825)
  # ) %>%
  ggplot(aes(x = month, y = n_inc_clinical_0_1825 / n_0_1825*1000)) +
  geom_line(col = "grey80") +
  # geom_vline(xintercept = program_start / year, col = "red") +
  stat_smooth(col = "darkblue", se = FALSE) +
  # facet_wrap(~ model, nrow = 3, ncol = 2) + 
  xlab("month") +
  ylab("Annual clinical incidence \n (per 1000 children aged 0-5)") +
  theme_bw()
clin_inc

prev <- output %>%
  # filter(timestep < year*6) %>%
  mutate(month = ceiling(timestep/30)) %>%
  group_by(month) %>%
  summarize(mo_n_detect_730_3650 = mean(n_detect_730_3650),
            mo_n_730_3650 = mean(n_730_3650),
            pfpr2_10 = mo_n_detect_730_3650 / mo_n_730_3650) 

ggplot(prev, aes(y = pfpr2_10, x = month)) +
  # geom_line() +
  stat_smooth(se = FALSE) +
  geom_vline(xintercept = 7, col = 'red') +
  xlab("month") +
  ylab("PfPR2-10") +
  scale_x_continuous(breaks = seq(-180, 146, by = 12)) +
  ylim(c(0.4,0.6)) +
  theme_bw()
