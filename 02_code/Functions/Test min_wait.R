# Test case to see if I can vaccinate a different population every 3 years (months to test)

# read in dataframe of all scenario combinations
scenarios <- readRDS(paste0(path, '03_output/scenarios_torun.rds'))

# pull one scenario at a time
data <- scenarios[264, ]# scenario with vaccination every 3 years and no MDA

# assign values
population = 1000
seasonality = data$seasonality
seas_name = data$seas_name
pfpr = data$pfpr
warmup = 1*365
sim_length = 9*365
RTSS = data$RTSS
RTSScov = data$RTSScov
RTSSage = data$RTSSage
RTSSrounds = data$RTSSrounds
fifth = data$fifth
booster_rep = data$booster_rep
MDAcov = data$MDAcov
MDAtiming = data$MDAtiming
ID = data$ID
drawID = data$drawID

year <- 365
month <- year / 12

# starting parameters ----------
params <- get_parameters(list(
  human_population = population,
  model_seasonality = TRUE,
  # rainfall fourier parameters
  g0 = unlist(seasonality)[1],
  g = unlist(seasonality)[2:4],
  h = unlist(seasonality)[5:7],
  individual_mosquitoes = FALSE))

# outcome definitions ----------
# Set clinical incidence rendering 
params$clinical_incidence_rendering_min_ages = c(seq(0, 95, by = 5)*year, 0)
params$clinical_incidence_rendering_max_ages = c(seq(5, 100, by = 5)*year, 100*year) 

# Set severe incidence rendering 
params$severe_incidence_rendering_min_ages = c(seq(0, 95, by = 5)*year, 0) 
params$severe_incidence_rendering_max_ages = c(seq(5, 100, by = 5)*year, 100*year) 

# prevalence 2-10 year olds
params$prevalence_rendering_min_ages = c(2 * year, 0 * year)
params$prevalence_rendering_max_ages = c(10 * year, 100 * year)

# demography ----------
flat_demog <- read.table(paste0(path,'/01_data/Flat_demog.txt')) # from mlgts
ages <- round(flat_demog$V3 * year) # top of age bracket
deathrates <- flat_demog$V5 / 365   # age-specific death rates


params <- set_demography(
  params,
  agegroups = ages,
  timesteps = 0,
  deathrates = matrix(deathrates, nrow = 1)
)

# RTSS ----------
program_start <- 1 * month

boost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) # coverage from 10.1016/S2214-109X(22)00416-8

min_ages = 5*year
max_ages = 100*year

# mass ----------

params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
peak <- peak_season_offset(params)

massboosters <- round(c(1 * year + 2 * month)) 
boost_cov <- c(RTSScov * 0.8, rep(RTSScov*0.8*0.9, length(massboosters)-1)) # coverage from 10.1016/S2214-109X(22)00416-8

# First set the EPI strategy 
EPIboosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month)) # from phase III R21 trial
pevtimesteps <- warmup + program_start # starting when warmup ends
EPIboost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) 

params <- set_pev_epi(
  parameters = params,
  profile = rtss_profile,
  timesteps = pevtimesteps,
  coverages = RTSScov,
  age = round(5 * month),
  min_wait = 10*year,
  booster_timestep = EPIboosters,
  booster_coverage = EPIboost_cov,
  booster_profile = list(rtss_booster_profile),
  seasonal_boosters = FALSE)

# Get timing for mass vaccination rounds
first <- round(warmup + program_start + (peak - month * 3.5), 0) # start a year later

pevtimesteps <- c(first, first + seq(3 * year, sim_length, 3 * year))

# Next, set the mass vaccination strategy in addition to routine vaccination
# if(RTSSage %in% c('young children', 'all children', 'under 5s')){
#   min_wait = 3 * year  # I have no idea what the min wait would actually be - if vaccinated around 6 m, then next would be around 3.5-4 yrs old.
# } else {
#   min_wait = 1 * month
# }

proppop_notpregnant <- 1 - 0.078/2 # from DHS data - see Get_pregnancy_rate.R

# add mass vaccination for taking into account pregnant women
params <- set_mass_pev(
  parameters = params,
  profile = rtss_profile,
  timesteps = pevtimesteps, 
  coverages = rep(RTSScov * proppop_notpregnant, length(pevtimesteps)), 
  min_ages = min_ages,
  max_ages = max_ages,
  min_wait = 0*year,
  booster_timestep = massboosters, # timesteps following initial vaccination 
  booster_profile = rep(list(rtss_booster_profile), length(massboosters)),
  booster_coverage = boost_cov#rep(boost_cov, length(massboosters)) # prop of vaccinated pop who will receive booster vaccine
)


params <- set_parameter_draw(params, drawID)



# save as data.frame
data$params <- list(params)


outtest <- run_simulation(timesteps = warmup+sim_length, parameters = params)

doses <- outtest |> 
  mutate(timestep = timestep-warmup,
    month = ceiling(timestep / 30)) |>
  filter(timestep >0) %>%
  select(month, starts_with('n_pev')) |>
  pivot_longer(starts_with('n_pev'), names_to = "dose", values_to = "count") |>
  group_by(month, dose) |>
  summarize(count = sum(count))

# plot doses
dosesplot <- ggplot(data = doses) + 
  geom_col(aes(x = month, y = count, fill = dose)) +
  xlab("month") +
  ylab("doses") +
  theme_bw() 
# scale_x_continuous(breaks = seq(180,550, by =20))
dosesplot
