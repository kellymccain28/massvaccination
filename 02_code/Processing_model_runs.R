# Processing the model runs from the HPC 

source('02_code/packages_data.R')
source(paste0(HPCpath, '02_code/Functions/HPC_processing.R'))
source(paste0(path, '02_code/Functions/calc_deaths_dalys.R'))


# Set up cluster ----------------------------------------------------------
setwd(HPCpath)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")#'kem22','Q:','//qdrive.dide.ic.ac.uk/homes','Q:')

config <- didehpc::didehpc_config(credentials = list(
  username = Sys.getenv("DIDE_USERNAME"),
  password = Sys.getenv("DIDE_PASSWORD")),
  workdir = HPCpath,
  shares = share,
  use_rrq = FALSE,
  cores = 1,
  cluster = "fi--dideclusthn", #"fi--didemrchnb",
  template = '8Core',#"32Core","GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core", "32Core"
  parallel = FALSE) 

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, '02_code/Functions/HPC_processing.R')),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

obj <- didehpc::queue_didehpc(ctx, config = config)
###############################################################################################
scenarios <- readRDS(paste0(HPCpath,'03_output/parameters_torun.rds'))
index <- c(1:nrow(scenarios)) # runs
index <- index[which(index$x > 9830),]
# Process model runs by 1 year age group  ---------------------------------

# Create function to aggregate the raw model runs over the whole simulation by ***1 year*** age group (HPC_summ()), add 
# mortality rate variables, uncertainty for outcomes, and DALY components, and save it
process_summ_runs <- function(x){
  start <- Sys.time()
  clean_summ <- HPC_summ(x) |>
    mortality_rate() |>
    outcome_uncertainty() |>
    daly_components() 
  
  print(paste0('Saving run ', x))
  saveRDS(clean_summ, paste0(HPCpath,'HPC_summarized/run_summ_', x, '_', clean_summ$int_ID[1], '.rds'))
  end <- Sys.time()
  print(end-start)
}
# index <- which(combo$MDAcov>0)

map_dfr(index, process_summ_runs)

# Process model runs by 5 year age group  ---------------------------------
# Create function to aggregate the raw model runs over the whole simulation by ***5 year*** age group (HPC_summ_5y()), add 
# mortality rate variables, uncertainty for outcomes, and DALY components, and save it
process_summ_runs_5y <- function(x){
  clean_summ <- HPC_summ_5y(x) |>
    mortality_rate() |>
    outcome_uncertainty() |>
    daly_components() 
  
  print(paste0('Saving run ', x))
  saveRDS(clean_summ, paste0(HPCpath,'HPC_summarized_5y/run_summ_', x, '_', clean_summ$int_ID[1], '.rds'))
}
# index <- index[which(index >9830),]
map_dfr(index, process_summ_runs_5y)

# Import and combine summarized files ------------------------------------------------
# Import the files that have been aggregated above by process_summ_runs to whole sim with 1 line per age group
files <- list.files(path = paste0(HPCpath, "HPC_summarized/"), pattern = "run_summ_*", full.names = TRUE) 
files_5y <- list.files(path = paste0(HPCpath, "HPC_summarized_5y/"), pattern = "run_summ_*", full.names = TRUE) 
dat_list <- lapply(files, function (x) readRDS(x))
dat_list_5y <- lapply(files_5y, function (x) readRDS(x))

# Bind the files together and add group ID for intervention scenario
dalyoutput <- bind_rows(dat_list) |>
  group_by(pfpr, seasonality, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
  mutate(group = cur_group_id()) |> ungroup()#, fill = TRUE, idcol = "identifier")
dalyoutput_5y <- bind_rows(dat_list_5y) |>
  group_by(pfpr, seasonality, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming, int_ID) |>
  mutate(group = cur_group_id()) |> ungroup()#, fill = TRUE, idcol = "identifier")
  
# save output
saveRDS(dalyoutput, paste0(HPCpath, "03_output/dalyoutput_draws.rds"))
saveRDS(dalyoutput_5y, paste0(HPCpath, "03_output/dalyoutput_draws_5y.rds"))

# Calculate outcomes averted ----------------------------------------------

source(paste0(path, "02_code/Functions/outcomes_averted.R"))
dalyoutput <- readRDS(paste0(HPCpath, "03_output/dalyoutput_draws.rds"))
dalyoutput_5y <- readRDS(paste0(HPCpath, "03_output/dalyoutput_draws_5y.rds"))
output <- outcomes_averted(dalyoutput)
output_5y <- outcomes_averted(dalyoutput_5y)

# save final output both on local machine and on shared drive
saveRDS(output, paste0(path, "03_output/scenarios_draws.rds"))
saveRDS(output, paste0(HPCpath, "scenarios_draws.rds"))

saveRDS(output_5y, paste0(path, "03_output/scenarios_draws_5y.rds"))
saveRDS(output_5y, paste0(HPCpath, "scenarios_draws_5y.rds"))


# Get credible intervals  -------------------------------------------------
# for incidence, severe incidence, incidence per 1000 doses, dalys, deaths
df <- readRDS(paste0(HPCpath, "scenarios_draws.rds"))
df_5y <- readRDS(paste0(HPCpath, "scenarios_draws_5y.rds"))
index <- seq(1, max(df$group), 1)

source(paste0(HPCpath,'02_code/Functions/HPC_processing.R'))
output <- map_dfr(index, get_cr_int, df)
output_5y <- map_dfr(index, get_cr_int, df_5y)

saveRDS(output, paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'.rds'))
saveRDS(output_5y, paste0(HPCpath,'./HPC_summarized/aggregatedoversim_1_',length(index),'_5y.rds'))


# Prep to plot outcomes averted ---------------------------------------------------
source(paste0(path, '02_code/Figures.R'))
output <- readRDS(paste0(HPCpath,'/HPC_summarized/aggregatedoversim_1_',length(index),'.rds'))
# output <- readRDS(paste0(HPCpath,'/HPC_summarized/aggregatedoversim_1_',length(index),'_5y.rds'))
prev <- unique(output$pfpr)

# write function to save plots 
save_averted_plots <- function(prev, type, seas, rtss){
  legend <- get_legend(
    plot1 + 
      theme(legend.position = "bottom",
            legend.title = element_blank())
  )
  top_row <- plot_grid(plot1 + theme(legend.position = 'none'), 
                       plot2 + theme(legend.position = 'none'), 
                       plot3 + theme(legend.position = 'none'), 
                       plot4 + theme(legend.position = 'none'), ncol = 2)
  p <- plot_grid(top_row, legend,
                      ncol = 1, rel_heights = c(1, 0.1))
  ggsave(paste0(path, '03_output/Figures/outcomesaverted', type, seas, rtss, prev, '.png'), width = 15, height = 10, units = 'in')
}

# Create plots for raw outcomes averted ------
outcomes_list = c('cases_averted', 'deaths_averted', 'dalys_averted', 'severe_averted')
prev <- unique(output$pfpr)
seas <- unique(output$seasonality)
rtss <- c('EPI','hybrid')

# plot outcomes averted 
for(pr in prev){
  for(s in seas){
    for(r in rtss){
      for(i in outcomes_list){
        print(paste0('Outcome: ', i, ' Seas: ', s, ' Prev: ', pr, ' RTSS: ', r))
        p <- plot_out_averted(outcome = i, prev = pr, seas = s, rtss = r)
        assign(paste0("plot", which(i == outcomes_list)), p)
      }
      print(paste0('Saving plot ', s, ' ', pr, ' ', r))
      save_averted_plots(prev = pr, type = 'raw', seas = s, rtss = r)
    }
  }
}

# Create plots for raw outcomes averted PER 1000 fully vaccinated individuals------
outcomes_list = c('cases_avertedper1000vax', 'deaths_avertedper1000vax', 'dalys_avertedper1000vax', 'severe_avertedper1000vax')

# plot outcomes averted 
for(pr in prev){
  for(s in seas){
    for(r in rtss){
      for(i in outcomes_list){
        print(paste0('Outcome: ', i, ' Seas: ', s, ' Prev: ', pr, ' RTSS: ', r))
        p <- plot_out_averted(outcome = i, prev = pr, seas = s, rtss = r)
        assign(paste0("plot", which(i == outcomes_list)), p)
      }
      print(paste0('Saving plot ', s, ' ', pr, ' ', r))
      save_averted_plots(prev = pr, type = 'per1000', seas = s, rtss = r)
    }
  }
}

 # delete files that produced plots for perennial et hybrid which doesn't exist
to_be_deleted <- list.files(paste0(path, "03_output/Figures/"), pattern = "perennialhybrid")
file.remove(paste0(path, '03_output/Figures/',to_be_deleted))


# Table of total number of outcomes averted ------------------------------------
agg <- output |>
  ungroup() |>
  select(c(int_ID, dalys_averted_lower:u5_severe_avertedper1000vax_upper)) |>
  group_by(int_ID) |>
  summarise(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), sum)) |>
  distinct() |>
  tidyr::separate(int_ID, into = c('PfPR','Seasonality','Vaccination strategy','Coverage','Age group vaccinated','Rounds','Fifth', 'MDA','MDAtiming'), sep="_") |>
  group_by(PfPR, Seasonality) |>
  arrange(PfPR, Seasonality) |>
  mutate(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), as.numeric)) |>
  mutate(maxcases = ifelse(cases_averted_median == max(cases_averted_median), 1, 0),
         maxsev = ifelse(severe_averted_median == max(severe_averted_median), 1, 0),
         maxdaly = ifelse(dalys_averted_median == max(dalys_averted_median), 1, 0),
         maxdeath = ifelse(deaths_averted_median == max(deaths_averted_median), 1, 0),
         maxcases1000 = ifelse(cases_avertedper1000vax_median == max(cases_avertedper1000vax_median), 1, 0),
         maxsev1000 = ifelse(severe_avertedper1000vax_median == max(severe_avertedper1000vax_median), 1, 0),
         maxdaly1000 = ifelse(dalys_avertedper1000vax_median == max(dalys_avertedper1000vax_median), 1, 0),
         maxdeath1000 = ifelse(deaths_avertedper1000vax_median == max(deaths_avertedper1000vax_median), 1, 0)) |>
  mutate(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), ~format(round(.x, digits = 1), big.mark = ","))) |>
  mutate(`Cases averted` = paste0(cases_averted_median, " (", cases_averted_lower, ", ", cases_averted_upper, ")"),
         `Severe cases averted` = paste0(severe_averted_median, " (", severe_averted_lower, ", ", severe_averted_upper, ")"),
         `DALYs averted` = paste0(dalys_averted_median, " (", dalys_averted_lower, ", ", dalys_averted_upper, ")"),
         `Deaths averted` = paste0(deaths_averted_median, " (", deaths_averted_lower, ", ", deaths_averted_upper, ")"),
         `Cases averted per 1000 vaccinated` = paste0(cases_avertedper1000vax_median, " (", cases_avertedper1000vax_lower, ", ", cases_avertedper1000vax_upper, ")"),
         `Severe cases averted per 1000 vaccinated` = paste0(severe_avertedper1000vax_median, " (", severe_avertedper1000vax_lower, ", ", severe_avertedper1000vax_upper, ")"),
         `DALYs averted per 1000 vaccinated` = paste0(dalys_avertedper1000vax_median, " (", dalys_avertedper1000vax_lower, ", ", dalys_avertedper1000vax_upper, ")"),
         `Deaths averted per 1000 vaccinated` = paste0(deaths_avertedper1000vax_median, " (", deaths_avertedper1000vax_lower, ", ", deaths_avertedper1000vax_upper, ")")) 

maxcases <- which(agg$maxcases == 1) 
maxsev <-   which(agg$maxsev == 1) 
maxdaly <- which(agg$maxdaly == 1) 
maxdeath <- which(agg$maxdeath == 1)
maxcases1000 <- which(agg$maxcases1000 == 1) 
maxsev1000 <- which(agg$maxsev1000 == 1) 
maxdaly1000 <- which(agg$maxdaly1000 == 1)
maxdeath1000 <- which(agg$maxdeath1000 == 1) 

border_style = officer::fp_border(color="black", width=1)

agg_out <- agg |>
  select(PfPR:MDAtiming, `Cases averted`:`Deaths averted per 1000 vaccinated`, -Coverage, -Fifth) |>
  # filter(Seasonality == 'perennial') |>
  flextable()|>
  # autofit()|>
  width(j=c(1,2,3,4,5), width = 0.7) |>
  width(j=c(6,7,8,9), width = 2) |>
  width(j=c(10,11,12,13), width = 2) |>
  fontsize(i = 1, size = 12, part = "header") |>  # adjust font size of header
  bold(i = 1, bold = TRUE, part = "header")  |>   # adjust bold face of header
  border_remove() |>
  theme_booktabs() |>
  vline(j = c(7, 11), border = border_style) |>
  hline(i = c(10, 40, 50, 80, 90, 120, 130, 140, 170), border = border_style) |>
  # highlight max values for max cases averted
  bg(j = 8, i = maxcases, part = "body", bg = "#91c293") |>
  # highl8ight values for max severe cases averted
  bg(j = 9, i = maxsev, part = "body", bg = "#91c293") |>
  # highlight values for max dalys averted
  bg(j = 10, i = maxdaly, part = "body", bg = "#91c293") |>
  # highlight values for max deaths averted
  bg(j = 11, i = maxdeath, part = "body", bg = "#91c293") |>
  # highlight max values for max cases averted per 1000 vaccinated
  bg(j = 12, i = maxcases1000, part = "body", bg = "#91c293") |>
  # highlight values for max severe cases averted per 1000 vaccinated
  bg(j = 13, i = maxsev1000, part = "body", bg = "#91c293") |>
  # highlight values for max dalys averted per 1000 vaccinated
  bg(j = 14, i = maxdaly1000, part = "body", bg = "#91c293") |>
  # highlight values for max deaths averted per 1000 vaccinated
  bg(j = 15, i = maxdeath1000, part = "body", bg = "#91c293") 
  
agg_out

save_as_image(agg_out, path = paste0(path, "03_output/Figures/table_outcomes_averted_210.png"))


# Table of max values over all prevalences for each of the outcome --------
agg <- output |>
  ungroup() |>
  select(c(int_ID, dalys_averted_lower:u5_severe_avertedper1000vax_upper)) |>
  group_by(int_ID) |>
  summarise(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), sum)) |>
  distinct() |>
  tidyr::separate(int_ID, into = c('PfPR','Seasonality','Vaccination strategy','Coverage','Age group vaccinated','Rounds','Fifth', 'MDA','MDAtiming'), sep="_") |>
  # group_by(Seasonality) |>
  # arrange(Seasonality) |>
  mutate(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), as.numeric)) |>
  mutate(maxcases = ifelse(cases_averted_median == max(cases_averted_median), 1, 0),
         maxsev = ifelse(severe_averted_median == max(severe_averted_median), 1, 0),
         maxdaly = ifelse(dalys_averted_median == max(dalys_averted_median), 1, 0),
         maxdeath = ifelse(deaths_averted_median == max(deaths_averted_median), 1, 0),
         maxcases1000 = ifelse(cases_avertedper1000vax_median == max(cases_avertedper1000vax_median), 1, 0),
         maxsev1000 = ifelse(severe_avertedper1000vax_median == max(severe_avertedper1000vax_median), 1, 0),
         maxdaly1000 = ifelse(dalys_avertedper1000vax_median == max(dalys_avertedper1000vax_median), 1, 0),
         maxdeath1000 = ifelse(deaths_avertedper1000vax_median == max(deaths_avertedper1000vax_median), 1, 0)) |>
  mutate(across(c(dalys_averted_lower:u5_severe_avertedper1000vax_upper), ~format(round(.x, digits = 1), big.mark = ","))) |>
  mutate(`Cases averted` = paste0(cases_averted_median, " (", cases_averted_lower, ", ", cases_averted_upper, ")"),
         `Severe cases averted` = paste0(severe_averted_median, " (", severe_averted_lower, ", ", severe_averted_upper, ")"),
         `DALYs averted` = paste0(dalys_averted_median, " (", dalys_averted_lower, ", ", dalys_averted_upper, ")"),
         `Deaths averted` = paste0(deaths_averted_median, " (", deaths_averted_lower, ", ", deaths_averted_upper, ")"),
         `Cases averted per 1000 vaccinated` = paste0(cases_avertedper1000vax_median, " (", cases_avertedper1000vax_lower, ", ", cases_avertedper1000vax_upper, ")"),
         `Severe cases averted per 1000 vaccinated` = paste0(severe_avertedper1000vax_median, " (", severe_avertedper1000vax_lower, ", ", severe_avertedper1000vax_upper, ")"),
         `DALYs averted per 1000 vaccinated` = paste0(dalys_avertedper1000vax_median, " (", dalys_avertedper1000vax_lower, ", ", dalys_avertedper1000vax_upper, ")"),
         `Deaths averted per 1000 vaccinated` = paste0(deaths_avertedper1000vax_median, " (", deaths_avertedper1000vax_lower, ", ", deaths_avertedper1000vax_upper, ")")) 

border_style = officer::fp_border(color="black", width=1)

agg <- agg |>
  filter(maxcases == 1 | maxsev == 1 |maxdaly==1|maxdeath==1|maxcases1000==1|maxsev1000==1|maxdaly1000==1|maxdeath1000==100) 

maxcases <- which(agg$maxcases == 1) 
maxsev <-   which(agg$maxsev == 1) 
maxdaly <- which(agg$maxdaly == 1) 
maxdeath <- which(agg$maxdeath == 1)
maxcases1000 <- which(agg$maxcases1000 == 1) 
maxsev1000 <- which(agg$maxsev1000 == 1) 
maxdaly1000 <- which(agg$maxdaly1000 == 1)
maxdeath1000 <- which(agg$maxdeath1000 == 1) 

agg_out <- agg |>
  select(PfPR:Fifth, `Cases averted`:`Deaths averted per 1000 vaccinated`, -Coverage, -Fifth) |>
  flextable()|>
  width(j=c(1,2,3,4,5), width = 0.7) |>
  width(j=c(6,7,8,9), width = 2) |>
  width(j=c(10,11,12,13), width = 2) |>
  fontsize(i = 1, size = 12, part = "header") |>  # adjust font size of header
  bold(i = 1, bold = TRUE, part = "header")  |>   # adjust bold face of header
  border_remove() |>
  theme_booktabs() |>
  vline(j = c(7, 11), border = border_style) |>
  # hline(i = c(4), border = border_style) |>
  # highlight max values for max cases averted
  bg(j = 6, i = maxcases, part = "body", bg = "#91c293") |>
  # highlight values for max severe cases averted
  bg(j = 7, i = maxsev, part = "body", bg = "#91c293") |>
  # highlight values for max dalys averted
  bg(j = 8, i = maxdaly, part = "body", bg = "#91c293") |>
  # highlight values for max deaths averted
  bg(j = 9, i = maxdeath, part = "body", bg = "#91c293") |>
  # highlight max values for max cases averted per 1000 vaccinated
  bg(j = 10, i = maxcases1000, part = "body", bg = "#91c293") |>
  # highlight values for max severe cases averted per 1000 vaccinated
  bg(j = 11, i = maxsev1000, part = "body", bg = "#91c293") |>
  # highlight values for max dalys averted per 1000 vaccinated
  bg(j = 12, i = maxdaly1000, part = "body", bg = "#91c293") |>
  # highlight values for max deaths averted per 1000 vaccinated
  bg(j = 13, i = maxdeath1000, part = "body", bg = "#91c293") 

agg_out
save_as_image(agg_out, path = paste0(path, "03_output/Figures/table_outcomes_averted_maxoverall_210.png"))


# Plot total outcomes ----
source(paste0(path, '02_code/Figures.R'))
df <- readRDS(paste0(HPCpath,'/HPC_summarized/aggregatedoversim_1_',length(index),'.rds'))

casesav <- plot_total_averted(df, outcome = "cases_averted", total = TRUE, labtitle = 'Total cases averted')
deathsav <- plot_total_averted(df, outcome = "deaths_averted", total = TRUE, labtitle = 'Total deaths averted')
dalysav <- plot_total_averted(df, outcome = "dalys_averted", total = TRUE, labtitle = 'Total DALYs averted')
severeav <- plot_total_averted(df, outcome = "severe_averted", total = TRUE, labtitle = 'Total severe cases averted')

casesavper1000 <- plot_total_averted(df, outcome = "cases_avertedper1000vax", total = TRUE, labtitle = 'Total cases averted per 1000 fully vaccinated people')
deathsavper1000 <- plot_total_averted(df, outcome = "deaths_avertedper1000vax", total = TRUE, labtitle = 'Total deaths averted per 1000 fully vaccinated people')
dalysavper1000 <- plot_total_averted(df, outcome = "dalys_avertedper1000vax", total = TRUE, labtitle = 'Total DALYs averted per 1000 fully vaccinated people')
severeavper1000 <- plot_total_averted(df, outcome = "severe_avertedper1000vax", total = TRUE, labtitle = 'Total severe cases averted per 1000 fully vaccinated people')


# ANNUALLY ###################################################################################################################################################
# Processing the output by year ------------------------------------------

## Run the summarize function to get annual prevalence
scenarios <- readRDS(paste0(HPCpath,'03_output/parameters_torun.rds'))
index <- c(1:nrow(scenarios)) # runs

# source(paste0(HPCpath, '02_code/')

# run a test with the first scenario
t <- obj$enqueue_bulk(1:2, process_runs_byyr) #936-svmass+hybrid #545
t$status()
#t$results()

# submit jobs, 100 as a time
sjob <- function(x, y){
  t <- obj$enqueue_bulk(index[x:y], process_runs_byyr)
  return(1)
}

map2_dfr(seq(0, length(index)- 100, 100),
         seq(99, length(index), 100),
         sjob)
#Submit last jobs
# index <- 5300:5355# 1800:1836
obj$enqueue_bulk(index, process_runs_byyr)

# map_dfr(index, process_runs_byyr)


# Get 95% CrI for prevalence by year ----------------------------------------------
# Get the file names for the annual summarized outputs and import/group by scenario ID to get 95% credible intervals 
# for prevalence
files <- list.files(path = paste0(HPCpath, "HPC_summbyyr/"), pattern = "run_summ*", full.names = TRUE) 
files <- data.frame(filenames = files) |>
  mutate(file = filenames) |>
  tidyr::separate(filenames, into = c('HPC',"run","summ",'draw','pfpr','seas','RTSS','RTSScov','RTSSage','RTSSrounds','fifth'), sep="_") %>%
  group_by(pfpr, seas, RTSS, RTSScov, RTSSage, RTSSrounds, fifth) |>
  mutate(group = cur_group_id()) 

index <- seq(1, max(files$group),1)

# # submit jobs, 10 as a time
# sjob <- function(x, y){
#   t <- obj$enqueue_bulk(index[x:y], get_annu_cr_int)
#   return(1)
# }
# 
# obj$enqueue_bulk(index, get_annu_cr_int)

source(paste0(HPCpath,'02_code/Functions/HPC_processing.R'))
output <- map_dfr(index, get_annu_cr_int)

saveRDS(output, paste0(HPCpath,'./HPC_summbyyr/outcomes_byyear_1_',length(index),'.rds'))


# Plotting prevalence -----------------------------------------------------
source(paste0(path, '02_code/Figures.R'))
# over time by year
output <- readRDS(paste0(HPCpath,'HPC_summbyyr/outcomes_byyear_1_',length(index),'.rds'))
prev <- unique(output$pfpr)

lapply(prev, outcome = 'prev_byyear', plot_annual_outcome)

# Plotting clinical incidence 0-5 ---------------------------------------------
lapply(prev, outcome = 'inc_0_1825', plot_annual_outcome)

# Plotting incidence among people 5-100 years ---------------------------------
lapply(prev, outcome = 'inc_1825_365000', plot_annual_outcome)

# Plotting severe incidence among people 5-100 years  ---------------------------------
lapply(prev, outcome = 'sev_1825_365000', plot_annual_outcome)

# Plotting incidence among people 5-15 years ---------------------------------
lapply(prev, outcome = 'inc_1825_5475', plot_annual_outcome)

# Plotting severe incidence among people 5-15 years  ---------------------------------
lapply(prev, outcome = 'sev_1825_5475', plot_annual_outcome)




####################################################################################################################################################

# Process model runs to show doses ----------------------------------------
scenarios <- readRDS(paste0(HPCpath,'03_output/parameters_torun.rds'))
index <- c(1:nrow(scenarios)) # runs
index <- 1:5355

# Save summarized df for each of the runs
lapply(index, summ_fordose)

#get credible intervals and combine together
files <- list.files(path = paste0(HPCpath, "HPC_5355/HPC_bydose/"), pattern = "run_bydose_*", full.names = TRUE) 
# files_missing <- strsplit(files, "_") 
# files_miss <- lapply(files_missing, `[[`, 5) %>% unlist() %>% as.numeric() 
# sequence <- min(files_miss):max(files_miss)
# index <- sequence[!sequence %in% files_miss] # these are the missing files (but they are missing because no vaccination )

dat_list <- lapply(files, function (y) readRDS(y))

# Bind the files together and add group ID for intervention scenario
out <- bind_rows(dat_list) |>
  group_by(int_ID, month, dose) |> 
  mutate(group = cur_group_id()) |> ungroup()#, fill = TRUE, idcol = "identifier")

# save output
saveRDS(out, paste0(HPCpath, "03_output/output_5355_preMDA/output_bydoses_summ_draws.rds"))

out <- readRDS(paste0(HPCpath, "03_output/output_5355_preMDA/output_bydoses_summ_draws.rds")) |>
  group_by(int_ID, month, dose) |> 
  mutate(group = cur_group_id()) |> ungroup()
index <- 1:max(out$group)

summarized_cr_doses <- lapply(index, get_cr_int_doses)
summarized_cr_doses_agg <- bind_rows(summarized_cr_doses)

saveRDS(summarized_cr_doses_agg, paste0(HPCpath,'HPC_5355/HPC_bydose/outcomes_byyear_1_',length(index),'.rds'))

# Plot vaccination strategies and doses ----
df <- readRDS(paste0(HPCpath,'HPC_5355/HPC_bydose/outcomes_byyear_1_',length(index),'.rds'))
plot_doses(df, seas = 'seasonal', prevalence = 0.03)
plot_doses(df, seas = 'perennial', prevalence = 0.03)

####################################################################################################################################################

# Process model runs to get 1 line per age group per year ----------------------------------------
## Run the summarize function to get annual prevalence
scenarios <- readRDS(paste0(HPCpath,'03_output/parameters_torun.rds'))
# index <- c(1:nrow(scenarios)) # runs
index <- 408:5355
#######################
########################## **** need to replace rtss with pev in function***********
#######################
map_dfr(index, process_age_yr)

# run a test with the first scenario
t <- obj$enqueue_bulk(1:2, process_age_yr) 
t$status()
#t$results()

# submit jobs, 100 as a time
sjob <- function(x, y){
  t <- obj$enqueue_bulk(index[x:y], process_age_yr)
  return(1)
}

map2_dfr(seq(0, length(index)- 100, 100),
         seq(99, length(index), 100),
         sjob)
#Submit last jobs
index <- 5300:5355# 1800:1836
obj$enqueue_bulk(index, process_age_yr)

# Get credible intervals for the runs by age/by year
files <- list.files(path = paste0(HPCpath, "HPC_5355/HPC_summbyyrbyage/"), pattern = "run_summ_*", full.names = TRUE) 
dat_list <- lapply(files, function (y) readRDS(y))

# Bind the files together and add group ID for intervention scenario
out <- bind_rows(dat_list) |>
  group_by(int_ID, year, age_grp) |> 
  mutate(group = cur_group_id()) |> ungroup()#, fill = TRUE, idcol = "identifier")

# save output
saveRDS(out, paste0(HPCpath, "03_output/output_5355_preMDA/output_byyear_byage_draws.rds"))

# get credible intervals 
out <- readRDS(paste0(HPCpath, "03_output/output_5355_preMDA/output_byyear_byage_draws.rds"))
index <- 1:max(out$group)
map_dfr(index, get_age_yr_cr)

