# Figures for paper 
# with an outputted dataset over x years that have clinical incidence stratified by year 0-20 years of age and 0-5 years

# Doses
# i <- sample(1:1836, 20, replace = FALSE) # to randomly sample 20 indices from list of runs to plot doses
# lapply(i, pltd)



plot_doses <- function(df, seas, prevalence){
#' input: dataset by month (output from process_by_month)
#' output: plot with doses and prevalence 
  df2 <- df %>%
    tidyr::separate(int_ID, into = c('pfpr','seasonality','RTSS','RTSScov','RTSSage','RTSSrounds','fifth', 'MDA','MDAtiming'), 
                    sep="_") %>%
    mutate(int_ID = paste(RTSS, RTSScov, RTSSage, sep = '_'),
           prev_med = n_detect_median / n_median,
           prev_lower = n_detect_lower / n_lower,
           prev_upper = n_detect_upper / n_upper) %>%
    filter(MDAtiming =='none', seasonality == seas, pfpr == prevalence, RTSSrounds != 'every 3 years')
  
  RTSS <- df2$RTSS[1]
  RTSSage <- df2$RTSSage[1]
  RTSSrounds <- df2$RTSSrounds[1]
  seasonality <- df2$seasonality[1]
  
  # plot doses
  dosesplot <- ggplot(data = df2) + 
    geom_col(aes(x = month, y = dosecount_median, fill = dose)) +
    geom_line(aes(x = month, y = prev_med * 200000), color= 'darkgreen', linewidth = 1.3) +
    facet_wrap(~int_ID, scales = 'free') +
    labs(x = "Month",
         y = "Number of doses",
         fill = 'Dose',
         title = paste0(seasonality, ', ', pfpr)) +
    scale_y_continuous(name = "Number of doses",
                       sec.axis = sec_axis(trans = ~./200000,
                                           name = 'Prevalence 2-10')) + #2*max(df2$dosecount_median)
    theme_bw() 
  
  dosesplot
  # ggsave(paste0(path, '03_output/Figures/Doses/dosesplot',x,'.png'), width = 8, height = 4)
}
# plot_doses(df, seasonality = 'seasonal', prevalence = 0.03)

######
# Plot annual outcomes including parasite prevalence2-10  (n_detect_730_3650/n_730_3650)
# input: './HPC_summbyyr/pfpr_byyear_1_x.rds'; prev_seas which is a crossing of each variable you want with prev in 1st col and seas in 2nd
# values possible for outcome: prev_byyear, inc_0_1825, inc_0_91.25, inc_1825_365000, sev_1825_365000, inc_1825_5475, sev_1825_5475
plot_annual_outcome <- function(prev, seas_type, outcome){
  if(outcome == 'prev_byyear'){
    out_label <- 'PfPR 2-10'
  } else if(outcome == 'inc_0_1825'){
    out_label <- 'Incidence per 1000 people aged 0-5y'
  } else if(outcome == 'sev_0_1825'){
    out_label <- 'Severe incidence per 1000 people aged 0-5y'
  }else if(outcome == 'inc_0_91.25'){
    out_label <- 'Incidence per 1000 people aged 0-4 m'
  } else if(outcome == 'inc_1825_36500'){
    out_label <- 'Incidence per 1000 people aged 5-100y'
  } else if(outcome == 'sev_1825_36500'){
    out_label <- 'Severe incidence per 1000 people aged 5-100y'
  } else if(outcome == 'inc_1825_5475'){
    out_label <- 'Incidence per 1000 people aged 5-15y'
  } else if(outcome == 'sev_1825_5475'){
    out_label <- 'Severe incidence per 1000 people aged 5-15y'
  }

  df_plot <- output %>% filter(pfpr == prev) %>% filter(year < 18) %>% 
    filter(RTSS != 'SV') %>%
    filter(seasonality == seas_type) %>%
    mutate(facet_lab = paste0(RTSS,"_", RTSSage,"_", booster_rep,"_", MDAtiming),
           year = year - 2) # so that year0 is = program start akka 1 year
  
  plt <- ggplot(df_plot) + 
    geom_ribbon(aes(x = year, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]]), fill = '#b9ccd3', alpha = 0.7) +
    geom_line(aes(x = year, y = .data[[paste0(outcome, "_median")]]), color = '#17556d') + 
    facet_wrap(~facet_lab) + theme_bw() + 
    scale_x_continuous(breaks = seq(min(df_plot$year),max(df_plot$year), by = 1)) +
    labs(y = out_label, x = 'Year', title = paste0('Initial PfPR = ', prev, ", ", seas_type)) #+
    # geom_vline(aes(xintercept = 1.5))
  plt
  
  
  ggsave(paste0(path, '03_output/Figures/', outcome, 'annual_', prev, "_", seas_type, '.png'), width = 16, height = 9, units = 'in')
}

#####
# Plot outcomes averted (Deaths and DALYs, severe cases, clinical cases)
plot_out_averted <- function(outcome, prev, seas, rtss){
  # input: paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'.rds')
  # possible values for outcome: dalys_averted, cases_averted, deaths_averted, severe_averted #u5_dalys_averted, u5_severe_averted, u5_cases_averted
  # output: plot of specified outcome averted by age group
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
      out_label <- 'Cases averted per 1000 \nfully vaccinated people'
    } else if(outcome == 'deaths_avertedper1000vax'){
      out_label <- 'Deaths averted per 1000 \nfully vaccinated people'
    } else if(outcome == 'severe_avertedper1000vax'){
      out_label <- 'Severe cases averted per \n1000 fully vaccinated people'
    } 
    output$int_ID_lab <- paste0(output$RTSS, " to ", output$RTSSage, ', ', output$booster_rep, '; ', output$MDAtiming)#str_sub(output$int_ID, 6, -3)
  
    plt <- ggplot(output %>% filter(pfpr == prev) %>% filter(seasonality == seas) %>% filter(grepl(rtss, RTSS))) +
      geom_col(aes(x = age_grp, y = .data[[paste0(outcome, "_median")]], group = int_ID_lab, fill = int_ID_lab), 
               color = 'black', position = 'dodge') + 
      geom_linerange(aes(x = age_grp, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]], group = int_ID_lab), 
                     position = position_dodge(width = 0.9), linewidth = 0.3) +
      # facet_wrap(~int_ID) + 
      theme_bw() + 
      labs(y = out_label, x = 'Age group', title = paste0('Initial PfPR = ', prev, ', ', seas)) +
      theme(axis.text.x=element_text(angle=90),
            legend.title=element_blank(),
            legend.position = c(0.8, 0.8))
    
    ggsave(paste0(path, '03_output/Figures/', i, prev, seas, rtss, '.png'), width = 15, height = 10, units = 'in')
    return(plt)
}


# Plot total outcomes averted over whole simulation or annually ; grouped by vaccine strategy and faceted by seasonality 
plot_total_averted <- function(df, outcome, total = TRUE, labtitle){
  #' input: dataset processed by HPC_summ -- paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'.rds')
  #' process: outcomes can be dalys_averted, cases_averted, deaths_averted, severe_averted, also per1000vax
  #' output: plot of total outcomes averted grouped by vaccine strategy and faceted by seasonality 
  
  # df <- readRDS(paste0(HPCpath,'HPC_5355/HPC_summarized/aggregatedoversim_1_105.rds'))
  lab_col <- c("#450061", "#7C00AD", # routine (EPI or hybrid)
               "#AD0911", "#E3595F", "#FA55B0", "#A83284", # mass + EPI
               "#87B3F5", "#6580A8", "#394EA8", "#283675", # SVmass + EPI
               "#87F6FA", "#66D1B8", "#29A8AD", "#005E61", # SV mass + hybrid
               "#48610A", "#BAFA19",
               "#F5E57B", "#FFD83E")
  # first need to group the dataset by int_ID
  df2 <- df %>%
    select(!starts_with("u5")) %>%
    dplyr::mutate(int_ID = paste(RTSS, RTSScov, RTSSage, RTSSrounds, MDAcov, sep = '_')) %>%
    filter(RTSS !="none")
  
  if (total == TRUE){
    df2 <- df2 %>%
      group_by(int_ID, seasonality, pfpr) %>%
      # get sum over all of the years
      summarise(across(c(dalys_averted_lower:severe_avertedper1000vax_upper), sum)) %>% 
      distinct() %>% ungroup()
    
    ylab <- 'Total outcomes averted'
    lab <- 'total'
  } else if (total == FALSE) { # mean annual outcomes averted
    df2 <- df2 %>%
      group_by(int_ID, seasonality, pfpr) %>%
      # get mean over all years of simulation 
      summarise(across(c(dalys_averted_lower:severe_avertedper1000vax_upper), mean)) %>% 
      distinct() %>% ungroup()
    
    ylab <- 'Mean annual outcomes averted'
    lab <- 'annual_mean'
  }
  
  plt <- ggplot(df2) +
    geom_col(aes(x = as.factor(pfpr), y = .data[[paste0(outcome, "_median")]], fill = int_ID), 
             color = 'black', position = 'dodge') + 
    geom_linerange(aes(x = as.factor(pfpr), ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]], group = int_ID), 
                   position = position_dodge(width = 0.9), linewidth = 0.3) +
    facet_wrap(~seasonality) + 
    scale_fill_manual(values = lab_col) +
    labs(y = ylab,
         x = "Initial PfPR",
         title = labtitle,
         fill = 'Intervention ID') +
    theme_bw()
  
  ggsave(paste0(path, '03_output/Figures/', lab, outcome, 'averted.png'), width = 15, height = 10, units = 'in')
  return(plt)
}

plot_cases_averted_by_age <- function(df){
  
}

# Incidence by 1 year age group over time ----
plot_inci_age <- function(df){
  #' input: HPC_summ output paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'.rds')
  #'  age_split: 1yr or 5yr
  #' output: line plot of x = year, y = incidence, grouped by age group  
  
  # df <- readRDS(paste0(HPCpath,'HPC_5355/HPC_summbyyr/outcomes_byyear_1_105.rds'))
  ggplot(df) +
    geom_line(aes(x = year, y = inc_clinical, group = age_grp, color = age_grp))
}


# Function to plot all runs of a specific scenario to see if there are any reaching elimination 
plot_all_runs <- function(){
  #' input: dataset with all draws combined together 
  #' process: plot the overall prevalence over time per scenario grouped by drawID
  #' output: plot with all 51 runs on same plot to visualize the number of runs that resulted in elimination  
  df <- readRDS(paste0(HPCpath, "HPC_summbyyrbyage/output_byyear_byage_draws.rds"))
  
  ggplot(df %>% filter(int_ID == "0.01_perennial_EPI_0.8_young children_none_0_0_0_none")) + 
    geom_line(aes(x = year, y = n_detect_0_36500/100000 * 1000, group = as.factor(drawID), color = as.factor(drawID))) +
    theme_bw() +
    theme(legend.position = 'none') 
  
  ggplot(d <- df %>% filter(int_ID == "0.01_perennial_mass+EPI_0.8_school-aged_single_0_0_0_none")) + 
    geom_line(aes(x = year, y = n_detect_0_36500/100000 * 1000, group = as.factor(drawID), color = as.factor(drawID))) +
    theme_bw() +
    theme(legend.position = 'none') 
  
  ggplot(d <- df %>% filter(int_ID == "0.01_perennial_mass+EPI_0.8_everyone_single_0_annual_0_none")) + 
    geom_line(aes(x = year, y = n_detect_0_36500/100000 * 1000, group = as.factor(drawID), color = as.factor(drawID))) +
    theme_bw() +
    theme(legend.position = 'none') 
  ggplot(d <- df %>% filter(int_ID == "0.03_seasonal_SVmass+EPI_0.8_everyone_single_0_annual_0_none")) + 
    geom_line(aes(x = year, y = n_detect_0_36500/100000 * 1000, group = as.factor(drawID), color = as.factor(drawID))) +
    theme_bw() +
    theme(legend.position = 'none') +
    scale_y_continuous(breaks = seq(0, 40, by = 0.5))
  ggplot(d <- df %>% filter(int_ID == "0.01_seasonal_SVmass+EPI_0.8_everyone_single_0_annual_0.8_during")) + 
    geom_line(aes(x = year, y = n_detect_0_36500/100000 * 1000, group = as.factor(drawID), color = as.factor(drawID))) +
    theme_bw() +
    theme(legend.position = 'none') +
    scale_y_continuous(breaks = seq(0, 16, by = 0.5))
}
