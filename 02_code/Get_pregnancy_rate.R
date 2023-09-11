# Troubleshooting removal of pregnant women

source("02_code/packages_data.R")
# devtools::install_github("ropensci/rdhs")
library(rdhs)

# What are the tags
tags <- dhs_tags()

indicators <- dhs_indicators()
id <- indicators[grepl("Percentage of the de facto household population age ", indicators$Definition), ]
ids <- id[11:27,'IndicatorId']
#FE_FRTY_W_NPG # number of women 15-49
preg <- dhs_data(countryIds = c('GM','BF','MW','MZ','SN'), indicatorIds = "FE_FRTY_W_PRG", surveyYearStart = 2015,breakdown = "background")
preg_country <- dhs_data(countryIds = c('GM','BF','MW','MZ','SN'),indicatorIds = "FE_FRTY_W_PRG", surveyYearStart = 2015)

pop <- dhs_data(countryIds = c('GM','BF','MW','MZ','SN'), indicatorIds = ids, surveyYearStart = 2015,breakdown = "background")
pop_country <- dhs_data(countryIds = c('GM','BF','MW','MZ','SN'),indicatorIds = ids, surveyYearStart = 2015) |>
  group_by(Indicator) |>
  summarize(mean_pop = mean(Value)) |> ungroup() |>
  filter(grepl("15", Indicator)| grepl("20", Indicator)| grepl("25", Indicator)| grepl("30", Indicator)| grepl("35", Indicator)| grepl("40", Indicator)|
         grepl("45", Indicator))
sum(pop_country$mean_pop)

mean(preg_country$Value) # 7.767
median(preg_country$Value) # 7.9
# let's say 7.8% of women between 15 and 49 are pregnant 
summ <- preg_country %>% group_by(CountryName) %>% summarize(mean_preg = mean(Value))
mean(summ$mean_preg)
## make a call with no arguments
# sc <- dhs_survey_characteristics()
# sc[grepl("Preg" , sc$SurveyCharacteristicName) | grepl("Wome" , sc$SurveyCharacteristicName), ]

## set up your credentials
# set_rdhs_config(email = "k.mccain22@imperial.ac.uk",
#                 project = "Understanding the public health impact of different methods of implementation of current and new vaccines against P.falciparum malaria in sub-Saharan Africa")

# lets find all the surveys that fit our search criteria
# survs <- dhs_surveys(surveyCharacteristicIds = 3,
#                      countryIds = c('GM','BF','MW','MZ','SN'),
#                      surveyType = "DHS",
#                      surveyYearStart = 2010)

# and lastly use this to find the datasets we will want to download and let's download the flat files (.dat) datasets (have a look in the dhs_datasets documentation for all argument options, and fileformat abbreviations etc.)
# datasets <- dhs_datasets(surveyIds = survs$SurveyId, 
#                          fileFormat = "flat")
# # the first time we call this function, rdhs will make the API request
# library(microbenchmark)
# microbenchmark::microbenchmark(dhs_surveys(surveyYear = 2010),times = 1)
# # download datasets
# downloads <- get_datasets(datasets$FileName)
# 
# bfbr <- readRDS(downloads$BFBR62FL)

# Parameterizing model 
year <- 365
month <- 30
sim_length <- 3 * year #12 for warmup and simlength
human_population <- 10000
starting_EIR <- 20
program_start <- round(0.5 * year)
warmup <- 2 * year


params <- get_parameters(list(
  human_population = human_population,
  model_seasonality = TRUE,
  #These seasonal params are from the seasonal values in Hillary's code (locally malariavax2G_2/02_code/B_PfPR_EIR_match.R)
  g0 = 0.285505,
  g = c(-0.325352,-0.0109352,0.0779865),
  h = c(-0.132815,0.104675,-0.013919),
  individual_mosquitoes = FALSE
)
)

# Set prevalence rendering
params$prevalence_rendering_min_ages = 2 * year
params$prevalence_rendering_max_ages = 10 * year

# Set clinical incidence rendering 
params$clinical_incidence_rendering_min_ages = c(0, 5*month)
params$clinical_incidence_rendering_max_ages = c(5 * year, 17*month)


peak <- peak_season_offset(params)

# One round of vaccination 
timesteps <- round(warmup + (1*year) + (peak - month * 3.5), 0)
boosters <- round(c(1*year))
params$rtss_doses <- round(c(0, 1 * month, 2 * month))

RTSScov <- 0.8
proppop_notpregnant <- 1 - 0.078/2
propwomen <- 0.5

params <- set_mass_rtss(
  params,
  timesteps = timesteps+5, 
  coverages = RTSScov, 
  min_wait = 1*month, 
  min_age = 5*month,
  max_age = 100*year, 
  boosters = boosters, #timesteps following initial vaccination 
  booster_coverage = rep(0.64, length(boosters)) # prop of vaccinated pop who will receive booster vaccine
)

# To increase booster titre up to the same as 3rd dose  (only for seasonal and hybrid)
params$rtss_cs_boost <- c(6.37008, 0.35) 

params <- set_equilibrium(params, starting_EIR)

output <- run_simulation(sim_length + warmup, params)

doses <- output |> 
  mutate(month = ceiling(timestep / 30)) |>
  select(month, starts_with('n_rtss')) |>
  pivot_longer(starts_with('n_rtss'), names_to = "dose", values_to = "count") |>
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
