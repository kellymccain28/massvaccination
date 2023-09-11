# Run scenarios in malariasimulation on HPC
### This script first sets up the HPC. 
### Next, creates a set of scenarios that we want to test and then uses the `generate_params` 
### function to create a list of parameters for each scenario and outputs this as a file. 
### The function `run_simulation` is run for each scenario in this parameter list.

source(paste0('C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/02_code/packages_data.R'))

## Model set up ------------------------------------------------------------------------------------------------
# year
year <- 365

# population
population <- 100000 # increased from 50,000 

# run time
warmup <- 15 * year       # 
sim_length <- 25*year  # value > 0  # should update from 21 to 25 years so that there is at least a 10 year follow up after 10y booster

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- c(0, 942,  40, 541, 497, 877, 697, 400, 450, 806, 
            600, 670, 363, 838, 478, 403, 375, 335, 598, 142,
            919, 444, 986, 659,  71, 457, 891, 188, 432, 975, 
            488, 867, 538, 912, 534, 215, 540, 866, 613, 973, 
            917, 937, 931, 296, 835, 328, 147, 701, 889, 708, 888)# c(0, sample(1:1000, 50))#c(0, 1:50)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.01, 0.03, 0.05, 0.25, 0.45, 0.65)
# pfpr <- c(0.01, 0.02, 0.03)

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- tibble(seasonality, seas_name)

seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- tibble(seasonality, seas_name)

seas_name <- 'perennial'
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
s3 <- tibble(seasonality, seas_name)

stable <- bind_rows( s2, s3)#s1,

# vectors
# list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

# INTERVENTIONS ----

# treatment coverage (baseline)
treatment <- c(0.45) 

# MDA
MDA <- c(0, 1)

# MDA timing
# MDAtiming <- c('none', 'during') 

# SMC: 0, 1
SMC <- c(0) 

# Vaccine 
PEV <- c('none','RTSS','R21')

# PEV strategy
PEVstrategy <- c('none', 'AB', 'hybrid', 'catch-up', 'mass')#SVmass+AB', 'SVmass+hybrid', 'mass+AB') #, 'SV'

# PEV coverage
PEVcov <- c(0, 0.8)

# RTS,S boost in immunogenicity for 4th dose (for AB with booster returning to same leve and for hybrid)
# RTSSboost <- c(0, 1)

# PEV age group
PEVage <- c('-', '5-15', '5-9', '5-100')#,'5-17m', 

# Rounds of PEV mass vaccination
PEVrounds <- c('-','single','3yrs') #, '5yrs'

# Status of fifth RTSS dose- for simplicity, can only be 0 
# fifth <- c(0)

# EPI booster timing 
EPIbooster <- c('12m', '12m boost', '18m', 'seasonal', '-')

# Extra boosters for EPI 
EPIextra <- c('5y', '10y', '-')

# adding vaccine boosters: 0 - no fifth dose, fifth (2 booster doses), annual (every year), 6mo (every 6 months), 2yrs (every 2 years)
massbooster_rep <- c('-', '4 annual', 'annual', '4 annual, no drop', 'annual no drop')


interventions <- crossing(treatment, SMC, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, "0", drawID, sep = "_")) 

# remove non-applicable scenarios 
combo <- combo |>
  filter(!(PEV == 'none' & PEVstrategy %in% c('AB', 'hybrid','catch-up', 'mass'))) |>
  filter(!(PEV %in% c('R21', 'RTSS') & PEVstrategy == 'none')) |> 
  filter(!(PEV == 'R21' & EPIbooster %in% c('12m boost', '18m'))) |> # only 12 month booster
  filter(!(PEVstrategy == 'none' & EPIbooster %in% c('12m', '12m boost', '18m', 'seasonal'))) |> # If no vaccination, then no EPI boosters
  filter(!(PEVstrategy %in% c('none', 'mass') & EPIextra %in% c('5y', '10y'))) |> #no extra epi boosters if no vax or if only mass
  filter(!(PEVstrategy == 'none' & PEVrounds %in% c('3yrs','single'))) |> # no rounds of PEVs if no vaccination
  filter(!(PEVstrategy == 'none' & PEVage %in% c('5-15', '5-9', '5-100'))) |> # if no vaccine, then age groups are irrelevant
  filter(!(PEVstrategy %in% c('AB', 'hybrid') & PEVage %in% c('5-15', '5-9', '5-100'))) |> # with EPI vax, no age groups specified
  filter(!(PEVstrategy == 'catch-up' & PEVage %in% c('-', '5-100'))) |> # catch-up vax to only 5-19 or 5-15
  filter(!(PEVstrategy == 'mass' & PEVage %in% c('5-15', '5-9', '-'))) |> # mass vax only to 5-100
  filter(!(PEVstrategy %in% c('AB', 'hybrid', 'catch-up') & PEVrounds %in% c('single','3yrs', '5yrs'))) |> # AB and hybrid would not have repeated mass rounds 
  filter(!(PEVstrategy == 'AB' & EPIbooster %in% c('seasonal','-'))) |> # no seasonal boosters for AB
  filter(!((PEV == 'RTSS' & PEVstrategy == 'catch-up') & EPIbooster %in% c('-', '12m', '18m','seasonal'))) |> # always has an EPI booster timing of 12 m boost
  filter(!((PEV == 'R21' & PEVstrategy == 'catch-up') & EPIbooster %in% c('-', '12m boost', '18m','seasonal'))) |> # always has an EPI booster timing of 12 m boost
  filter(!(PEVstrategy == 'hybrid' & EPIbooster %in% c('12m', '12m boost', '18m','-'))) |> # hybrid vaccination will not have these boosters
  filter(!(PEVstrategy == 'mass' & EPIbooster %in% c('12m', '12m boost', '18m','seasonal'))) |> # mass vaccination doesn't have epi boosters
  filter(!(PEVstrategy %in% c('AB', 'hybrid','mass', 'none') & EPIextra %in% c('5y','10y'))) |> # only catch-up vax gets extra boosters
  filter(!(seas_name == 'perennial' & PEVstrategy == 'hybrid')) |> # no hybrid vaccination in perennial settings 
  filter(!(PEVstrategy %in% c('AB', 'hybrid','catch-up') & MDA == 1)) |># MDA is only to PEVstrategy == none and PEVstrategy ==mass 
  filter(!(PEVstrategy %in% c('AB', 'hybrid','catch-up', 'mass') & PEVcov == 0)) |> # if vaccination, then coverage is not 0
  filter(!(PEVstrategy == 'none' & PEVcov == 0.8)) |>
  filter(!(PEVstrategy == 'mass' & PEVrounds == '-')) |># this means the same thing as single 
  filter(!(PEVstrategy %in% c('none', 'catch-up', 'hybrid', 'AB') & massbooster_rep %in% c('4 annual', 'annual', '4 annual, no drop', 'annual no drop'))) # only mass booster repetition in mass scenarios 
  
# put variables into the same order as function arguments
combo <- combo |> 
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         PEV,               # Which PEV
         PEVstrategy,       # PEV strategy
         PEVcov,            # PEV coverage
         PEVage,            # PEV age groups
         PEVrounds,         # PEV rounds of mass vax
         EPIbooster,        # Timing of EPI booster 
         EPIextra,          # Extra EPI boosters for catch-up campaigns 
         massbooster_rep,   # how many mass boosters
         MDA,               # MDA coverage
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()# |>
  # filter(PEV == 'RTSS' | PEV == 'none') 

# Check that it looks right 
check <- combo %>% group_by(PEV, PEVstrategy, PEVage, EPIextra, EPIbooster, massbooster_rep, PEVrounds, MDA) %>% summarize(n = n())
# new <- combo%>% filter(PEV == 'RTSS' & massbooster_rep%in%c('4 annual, no drop', 'annual no drop'))
  # filter(!(massbooster_rep %in% c('4 annual, no drop', 'annual no drop'))) 
# combo <- rbind(combo, new)
# new <- combo |> filter(massbooster_rep=='annual')
# combo <- rbind(combo, new)

saveRDS(combo, paste0(path, '03_output/scenarios_torun_RTSS.rds'))
# combo <- readRDS(paste0(path, '03_output/scenarios_torun_R21.rds'))
# generate parameter list in malariasimulation format
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/scenarios_torun_RTSS.rds'), # file path to pull
                paste0(HPCpath, "03_output/parameters_torun_RTSS.rds"))      # file path to push

## Setting up cluster ------------------------------------------------------------------------------------------------
setwd(HPCpath)

# didehpc::web_login()

# to edit HPC username and password below
# usethis::edit_r_environ()
share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")#'kem22','Q:','//qdrive.dide.ic.ac.uk/homes','Q:')

config <- didehpc::didehpc_config(credentials = list(
  username = 'kem22',#Sys.getenv("DIDE_USERNAME"),
  password = 'pwdRease1449!'),#Sys.getenv("DIDE_PASSWORD")),
  workdir = HPCpath,
  shares = share,
  use_rrq = FALSE,
  cores = 1,
  cluster =  "fi--didemrchnb",#"fi--dideclusthn",
  template = "GeneralNodes",#'8Core',"32Core", "12Core", "16Core", "12and16Core", "20Core", "24Core", "32Core"
  parallel = FALSE) 

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, '02_code/Functions/run_simulation.R')),
                             packages = c("dplyr", "malariasimulation", "purrr", "tidyr"),
                             package_sources = src)

obj <- didehpc::queue_didehpc(ctx, config = config)
# obj$install_packages("github::mrc-ide/malariasimulation@dev")

# Run tasks -------------------------------------------------------------------------------------------------------------
combo <- readRDS(paste0(HPCpath, "03_output/parameters_torun_test.rds"))
combo <- readRDS(paste0(HPCpath, "03_output/parameters_torun_RTSS.rds"))
combo <- readRDS(paste0(HPCpath, "03_output/parameters_torun_R21.rds"))
combo <- readRDS(paste0(HPCpath, "03_output/parameters_torun_elimination.rds"))
# x <- c(10,27,37,38, 43,44,79,80,85,86)
x = c(1:nrow(combo)) # runs

x <- which(combo$massbooster_rep == '4 annual'|combo$massbooster_rep == '4 annual, no drop')

# index <- tibble(x=x)
ii <- index %>%
 mutate(f = paste0(HPCpath, "HPC_RTSS_wrongorder/raw_modelrun_", index$x, ".rds"))
catchupfiles <- lapply(ii$f, function(m) readRDS(m)) %>% rbindlist(fill = TRUE, idcol = 'identifier')
catchupsimpl <- catchupfiles %>% group_by(identifier) %>% slice(1) %>% select(c(PEV, PEVstrategy, PEVage, EPIextra, EPIbooster))

catchupsimpl$catchup <- ifelse(catchupsimpl$PEVstrategy == 'catch-up', 1, 0)

# x <- x[900:960]
# x <- which((combo$PEVstrategy == 'mass' & combo$massbooster_rep != '-' & combo$PEVrounds !='single') |
#              (combo$PEVstrategy == 'catch-up' & combo$EPIextra != '-') |
#              combo$PEVstrategy == 'hybrid')

# define all combinations of scenarios and draws
index <- tibble(x = x)

# remove ones that have already been run
index <- index |>
  mutate(f = paste0(HPCpath, "HPC_RTSS/raw_modelrun_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with test scenario
t <- obj$enqueue_bulk(27, runsim) 
t$status()

# break tasks to run in groups of 1,000
# y <- c(seq(6401, 9500 - 100, 100), 9401)
# z <- seq(6500, 9500, 100)
# data <- tibble(y, z)
# t <- obj$enqueue_bulk(data, runsim)


# index <- seq(1, nrow(combo), 1)
# 
# # break tasks to run in groups of 100
# x <- c(seq(1, length(index) - 100, 100), 5001)
# y <- seq(100, length(index), 100)
# 
# data <- tibble(x, y)

# submit jobs, 100 as a time
sjob <- function(x, y){

  t <- obj$enqueue_bulk(index[x:y,], runsim)
  print(paste0(x, " to ", y))

}

map2_dfr(seq(0, nrow(index) - 100, 100),
         seq(99, nrow(index), 100),
         sjob)
 
# # submit all remaining tasks
t <- obj$enqueue_bulk(index, runsim)
t$status()
