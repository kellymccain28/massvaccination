# Run scenarios in malariasimulation on HPC
### This script first sets up the HPC. 
### Next, creates a set of scenarios that we want to test and then uses the `generate_params` 
### function to create a list of parameters for each scenario and outputs this as a file. 
### The function `run_simulation` is run for each scenario in this parameter list.

source('02_code/packages_data.R')

## Model set up ------------------------------------------------------------------------------------------------
# year
year <- 365

# population
population <- 100000 # increased from 50,000 

# run time
warmup <- 15 * year       # 
sim_length <- 21 * year   # value > 0 

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
MDAcov <- c(0, 1)

# MDA timing
# MDAtiming <- c('none', 'during') 

# SMC: 0, 1
SMC <- c(0) 

# RTS,S
RTSS <- c('none', 'AB', 'hybrid', 'catch-up', 'mass')#SVmass+AB', 'SVmass+hybrid', 'mass+AB') #, 'SV'

# RTS,S coverage
# RTSScov <- c(0, 0.8) 

# RTS,S boost in immunogenicity for 4th dose (for AB with booster returning to same leve and for hybrid)
# RTSSboost <- c(0, 1)
# don't need to boost here - will do in parameterization 

# RTS,S age group
RTSSage <- c('-', '5-15', '5-9', '5-100')#,'5-17m', 

# Rounds of RTSS mass vaccination
RTSSrounds <- c('-','single','3yrs') #, '5yrs'

# Status of fifth RTSS dose- for simplicity, can only be 0 
# fifth <- c(0)

# EPI booster timing 
EPIbooster <- c('12m', '18m', 'seasonal', '-')

# Extra boosters for EPI 
EPIextra <- c('5y', '10y', '-')

# adding vaccine boosters: 0 - no fifth dose, fifth (2 booster doses), annual (every year), 6mo (every 6 months), 2yrs (every 2 years)
# massbooster_rep <- c(0, 'annual')

interventions <- crossing(treatment, SMC, RTSS, RTSSage, RTSSrounds, EPIbooster, EPIextra, MDAcov)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, drawID, sep = "_")) 

# remove non-applicable scenarios 
combo <- combo |>
  filter(!(RTSS == 'none' & EPIbooster %in% c('12m', '18m', 'seasonal'))) |> # If no vaccination, then no EPI boosters
  filter(!(RTSS %in% c('none', 'mass') & EPIextra %in% c('5y', '10y'))) |> #no extra epi boosters if no vax or if only mass
  filter(!(RTSS == 'none' & RTSSrounds %in% c('3yrs','single'))) |> # no rounds of RTSs if no vaccination
  # filter(!(RTSS == 'none' & RTSSboost == 1)) |> # if no vaccine, then no boost in immunity 
  filter(!(RTSS == 'none' & RTSSage %in% c('5-15', '5-9', '5-100'))) |> # if no vaccine, then age groups are irrelevant
  filter(!(RTSS %in% c('AB', 'hybrid') & RTSSage %in% c('5-15', '5-9', '5-100'))) |> # with EPI vax, no age groups specified
  filter(!(RTSS == 'catch-up' & RTSSage %in% c('-', '5-100'))) |> # catch-up vax to only 5-19 or 5-15
  filter(!(RTSS == 'mass' & RTSSage %in% c('5-15', '5-9', '-'))) |> # mass vax only to 5-100
  filter(!(RTSS %in% c('AB', 'hybrid', 'catch-up') & RTSSrounds %in% c('single','3yrs', '5yrs'))) |> # AB and hybrid would not have repeated mass rounds 
  #filter(!(RTSS %in% c('catch-up', 'mass') & RTSSrounds == '-')) |> # need to specify number of mass rounds 
  filter(!(RTSS == 'AB' & EPIbooster %in% c('seasonal','-'))) |> # no seasonal boosters for AB
  filter(!(RTSS == 'catch-up' & EPIbooster %in% c('-', '18m','seasonal'))) |> # always has an EPI booster timing of 12 m
  filter(!(RTSS == 'hybrid' & EPIbooster %in% c('12m', '18m','-'))) |> # hybrid vaccination will not have these boosters
  filter(!(RTSS == 'mass' & EPIbooster %in% c('12m', '18m','seasonal'))) |> # mass vaccination doesn't have epi boosters
  filter(!(RTSS %in% c('AB', 'hybrid','mass', 'none') & EPIextra %in% c('5y','10y'))) |> # only catch-up vax gets extra boosters
  # filter(!(RTSS == 'none' & massbooster_rep == 'annual')) |> # if no vaccine, then no annual mass bosoters 
  filter(!(seas_name == 'perennial' & RTSS == 'hybrid')) |> # no hybrid vaccination in perennial settings 
  # filter(!(RTSSboost == 1 & RTSSrounds == '-')) |> # no boost if no rounds of RTSS
  # filter(!(RTSSboost == 1 & RTSSage == '-')) |> # no boost if no age group 
  # filter(!(RTSSboost == 1 & ))
  filter(!(RTSS %in% c('AB', 'hybrid','catch-up') & MDAcov == 1)) # MDA is only to rtss == none and rtss ==mass 
  
  
  # filter(!(RTSS=='none' & RTSScov > 0)) |>
  # filter(!(RTSS %in% c('mass','mass','mass','AB','hybrid','catch-up') & RTSScov == 0)) |>
  # filter(!(RTSScov == 0 & RTSSage %in% c('5-9','5-100','5-15'))) |> #,'5-17m'
  # filter(!(RTSScov > 0 & RTSSage == 'none')) |>
  # filter(!(RTSS %in% c('AB','hybrid') & RTSSage %in% c('5-9', '5-100','5-15','none'))) |> # age and hybrid and SV only given to young children 
  # filter(!(RTSS %in% c('AB','hybrid') & RTSSrounds %in% c('single','every 3 years'))) |> # they don't get mass vaccination rounds
  # filter(!(RTSSage == 'none' & RTSSrounds %in% c('single','every 3 years'))) |>
  # filter(!((RTSS =='none' | RTSS == 'AB') & fifth == 1)) |>
  # filter(!(seas_name == "perennial" & (RTSS == "SVmass+AB" |RTSS =='SVmass+hybrid' | RTSS == "hybrid" | RTSS == 'SV'))) |>
  # filter(!((RTSS == 'SVmass+AB' |RTSS =='SVmass+hybrid'| RTSS =='mass+AB') & RTSSrounds == 'none')) |>
  # filter(!((RTSS == 'SVmass+AB'|RTSS =='SVmass+hybrid'|RTSS == 'mass+AB') & RTSSage %in% c('5-17m'))) |> # when routine vaccination is to children 5-17 months, wouldn't do mass to them too
  # filter(!((MDAcov > 0 | MDAtiming == 'during') & RTSS == 'none')) |> # only testing MDA combined with RTSS or RTSS alone, not MDA on its own
  # filter(!(MDAcov > 0 & MDAtiming == 'none')) |> # can't have coverage of MDA and no timing of MDA
  # filter(!(MDAcov == 0 & MDAtiming == 'during')) |> # can't have coverage of MDA at 0% and MDA during vax
  # filter(!(seas_name == 'seasonal' & RTSS == 'mass+AB')) |> # no non-seasonal mass vaccination in seasonal areas
  # filter(!(booster_rep == 'annual' & RTSS %in% c('hybrid', 'AB','none'))) |> # no annual boosters for AB and hybrid vaccination or when no vax
  # filter(!(booster_rep %in% c('annual', '6mo', '2yrs') & RTSSrounds == 'every 3 years')) |># for now, there cannot be a repeated mass campaign
  # filter(!(booster_rep == 'annual' & RTSScov == 0)) |> # no annual boosters when no vaccination 
  # filter(!(booster_rep == 'annual' & RTSS == 'none')) |> # no annual boosters when no vaccination 
  # filter(!(RTSSboost == 1 & (RTSS == 'none' | RTSSrounds == 'none' | RTSSage == 'none' | RTSScov == 0)))  # can't have boosted immunogenicity if wasn't vaccinated 


  
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
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         RTSSage,           # RTS,S age groups
         RTSSrounds,        # RTS,S rounds of mass vax
         RTSSboost,         # RTS,S boosted immunogenicity for 4th dose
         booster_rep,       # how many boosters 
         MDAcov,            # MDA coverage
         MDAtiming,         # timing of MDA round
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()#|>
  # filter(RTSSrounds!='every 3 years')
# # rearrange so new ones are at the end 
# new <- combo |>
#   filter(RTSSrounds!='every 3 years' & booster_rep == 'annual')
# combo <- combo |>
#   filter(!(RTSSrounds!='every 3 years' & booster_rep == 'annual')) |>
#   rbind(new)


saveRDS(combo, paste0(path, '03_output/scenarios_torun.rds'))

# generate parameter list in malariasimulation format
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/scenarios_torun.rds'), # file path to pull
                paste0(HPCpath, "03_output/parameters_torun.rds"))      # file path to push

## Setting up cluster ------------------------------------------------------------------------------------------------
setwd(HPCpath)

# didehpc::web_login()

# to edit HPC username and password below
# usethis::edit_r_environ()
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

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation@dev"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, '02_code/Functions/run_simulation.R')),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

obj <- didehpc::queue_didehpc(ctx, config = config)
# obj$install_packages("github::mrc-ide/malariasimulation@dev")

# Run tasks -------------------------------------------------------------------------------------------------------------
x = c(1:nrow(combo)) # runs
# x <- which(combo$MDAcov>0 ) # to get only indices where there is MDA because we only need to rerun those
# x <- 11049:(nrow(combo)) # just the new ones (with annual boosters and 0.01)
# define all combinations of scenarios and draws
index <- tibble(x = x)

# remove ones that have already been run
index <- index |>
  mutate(f = paste0(HPCpath, "HPC/raw_modelrun_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with the first scenario
t <- obj$enqueue_bulk(4020, runsim) #936-svmass+hybrid #545
t$status()
# t$wait(1000)
#t$results()

# submit jobs, 100 as a time
sjob <- function(x, y){
  
  t <- obj$enqueue_bulk(index[x:y,], runsim)
  return(1)
  
}

map2_dfr(seq(0, nrow(index)- 100, 100),
         seq(99, nrow(index), 100),
         sjob)
# map2_dfr(3998, 4080, sjob)
# submit all remaining tasks
t <- obj$enqueue_bulk(index, runsim)
# t$status()
