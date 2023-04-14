# Figures for paper 
# with an outputted dataset over x years that have clinical incidence stratified by year 0-20 years of age and 0-5 years

# Doses
# i <- sample(1:1836, 20, replace = FALSE) # to randomly sample 20 indices from list of runs to plot doses
# lapply(i, pltd)

plot_doses <- function(){
  out <- readRDS(paste0(HPCpath, 'HPC/raw_modelrun_', format(x, scientific = FALSE), '.rds'))
  
  doses <- out |> 
    # mutate(#timestep = timestep#-warmup,
    #   month = ceiling(timestep / 30)) |>
    filter(timestep >0) %>%
    select(month, starts_with('n_rtss')) |>
    pivot_longer(starts_with('n_rtss'), names_to = "dose", values_to = "count") |>
    group_by(month, dose) |>
    summarize(count = sum(count))

  # plot doses
  dosesplot <- ggplot(data = doses) + 
    geom_col(aes(x = month/12, y = count, fill = dose)) +
    xlab("month") +
    ylab("doses") +
    theme_bw() 
  
  dosesplot
  ggsave(paste0(path, '03_output/Figures/Doses/dosesplot',x,'.png'), width = 8, height = 4)
}


######
# Plot annual outcomes including parasite prevalence2-10  (n_detect_730_3650/n_730_3650)
# input: './HPC_summbyyr/pfpr_byyear_1_x.rds'
# values possible for outcome: prev_byyear, inc_0_1825, inc_0_91.25, inc_1825_365000, sev_1825_365000, inc_1825_5475, sev_1825_5475
plot_annual_outcome <- function(prev, outcome){
  if(outcome == 'prev_byyear'){
    out_label <- 'PfPR 2-10'
  } else if(outcome == 'inc_0_1825'){
    out_label <- 'Incidence per 1000 people aged 0-5y'
  } else if(outcome == 'inc_0_91.25'){
    out_label <- 'Incidence per 1000 people aged 0-4 m'
  } else if(outcome == 'inc_1825_365000'){
    out_label <- 'Incidence per 1000 people aged 5-100y'
  } else if(outcome == 'sev_1825_365000'){
    out_label <- 'Severe incidence per 1000 people aged 5-100y'
  } else if(outcome == 'inc_1825_5475'){
    out_label <- 'Incidence per 1000 people aged 5-15y'
  } else if(outcome == 'sev_1825_5475'){
    out_label <- 'Severe incidence per 1000 people aged 5-15y'
  }
  
  df_plot <- output %>% filter(pfpr == prev) %>% filter(year < 6)
  
  plt <- ggplot(df_plot) + 
    geom_ribbon(aes(x = year, ymin = .data[[paste0(outcome, "_lower")]], ymax = .data[[paste0(outcome, "_upper")]]), fill = '#b9ccd3', alpha = 0.7) +
    geom_line(aes(x = year, y = .data[[paste0(outcome, "_median")]]), color = '#17556d') + 
    facet_wrap(~int_ID) + theme_bw() + 
    labs(y = out_label, x = 'Year', title = paste0('Initial PfPR = ', prev)) 
  # if(outcome == 'prev_byyear'){
  #   plt <- plt + geom_vline(aes(xintercept = df_plot[df_plot$prev_byyear_median >= 0.9 * max(df_plot$pfpr)]))
  # }
  
  ggsave(paste0(path, '03_output/Figures/', outcome, 'annual_', prev, '.png'), width = 16, height = 9, units = 'in')
}

#####
# Plot outcomes averted (Deaths and DALYs, severe cases, clinical cases)
# input: paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'.rds')
# possible values for outcome: dalys_averted, cases_averted, deaths_averted, severe_averted #u5_dalys_averted, u5_severe_averted, u5_cases_averted
plot_out_averted <- function(outcome, prev, seas, rtss){
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
    output$int_ID_lab <- paste0(output$RTSS, " to ", output$RTSSage, ', ', output$RTSSrounds)#str_sub(output$int_ID, 6, -3)
  
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
