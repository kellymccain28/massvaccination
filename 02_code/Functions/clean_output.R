clean_output <- function(){
  output <- output %>%
    filter(timestep > warmup) %>%
    mutate(timestep = timestep - warmup,
           year = ceiling(timestep/year),
           month = ceiling(timestep/month)) %>%
    select(c(year, month, 
             starts_with('n_inc_cl'), starts_with('p_inc_cl'),
             # starts_with('n_inc_sev'), starts_with('p_inc_sev'),
             starts_with('n_rtss'), starts_with('n_'),
             -n_bitten, p_detect_730_3650, -timestep)) %>%
    group_by(month, year) %>%
    mutate_at(vars(n_0_153.3:n_7300_9125, n_730_3650,
                   n_detect_730_3650, p_detect_730_3650), mean, na.rm = TRUE) %>%
    mutate_at(vars(n_inc_clinical_0_153.3:p_inc_clinical_7300_9125,
                   # n_inc_severe_0_153.3:p_inc_severe_7300_9125,
                   n_treated, n_infections), sum, na.rm = TRUE) %>%
    pivot_longer(cols = c(n_0_153.3:n_7300_9125,
                          n_inc_clinical_0_153.3:n_inc_clinical_7300_9125,
                          # n_inc_severe_0_153.3:n_inc_severe_7300_9125
                          ),
                 names_to = c('age'), values_to = c('value')) %>%
    mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
           inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
           # inc_severe = ifelse(grepl('n_inc_severe', age), value, NA), # creating var for inc_severe
           age = gsub('n_inc_clinical_', '', age),                         # combining age vars
           # age = gsub('n_inc_severe_', '', age),
           age = gsub('n_', '', age),
           age = gsub('_', '-', age)) %>% 
    group_by(age, year, month) %>%
    select(-c(value)) %>%
    mutate_at(vars(n:inc_clinical), sum, na.rm = TRUE) %>% # consolidate
    distinct() %>% ungroup() %>%
    separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
    mutate(age_lower = as.numeric(age_lower)/365,
           age_upper = as.numeric(age_upper)/365,
           age_grp = paste(age_lower, age_upper, sep = '-'),
           age_grp = factor(age_grp, levels =c('0-0.42', '0.42-5','5-10','10-15','15-20','20-25')),
           inc = inc_clinical / n,
           # sev = inc_severe / n,
           cases = inc_clinical)
}