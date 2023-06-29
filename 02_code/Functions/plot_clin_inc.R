plot_clin_inc <- function(){
  starting_EIRs <- c(5,20,50)
  doses <- output %>% 
    group_by(year, month, n_rtss_mass_dose_1, n_rtss_mass_dose_2, n_rtss_mass_dose_3, n_rtss_mass_booster_1) %>%
    summarize() %>%
    pivot_longer(!c(year, month), names_to = 'dose', values_to = 'n') %>%
    filter(n!=0) %>% ungroup() %>% select(month) %>% pull
  
  return(ggplot(data = output, aes(x = month, y = inc*1000)) + 
    geom_line(col = 'grey80') +
    stat_smooth(col = 'darkblue',se = T) +
    geom_vline(xintercept = doses, col = "red") +
    facet_wrap(~age_grp, labeller = as_labeller(age_names), scales = "free") +
    xlab("Month after warmup period") +
    ylab("Clinical incidence \n (per 1000 children") +
    theme_bw()
  )
}