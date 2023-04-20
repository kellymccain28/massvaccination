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
population <- 50000

# run time
warmup <- 15 * year       # needs to be multiple of 3 for ITN distribution
sim_length <- 21 * year   # value > 0 

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- c(0, 1:50)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.03, 0.05, 0.25, 0.45, 0.65)#c(0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65) # 

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
# ITN type: pyr, pbo
ITN <- c('pyr') 

# ITN usage
ITNuse <- c(0)#0.25, 0.50, 

# boosting ITN use by 10%: 0, 1
ITNboost <- c(0) 

# insecticide resistance level: 0, 0.4, 0.8
resistance <- c(0) 

# IRS coverage
IRS <-  c(0)

# treatment coverage (baseline)
treatment <- c(0.45) 

# MDAcoverage
MDAcov <- c(0, 0.8)

# MDA timing
MDAtiming <- c('none', 'during') #'prior', , 'after'

# SMC: 0, 1
SMC <- c(0) 

# RTS,S: none, EPI, SVmass+EPI, SVmass+hybrid, mass+EPI, hybrid, catch-up
RTSS <- c('none', 'SVmass+EPI', 'SVmass+hybrid', 'EPI', 'hybrid', 'mass+EPI', 'SV') #,'catch-up'

# RTS,S coverage
RTSScov <- c(0, 0.8) 

# RTS,S age group
RTSSage <- c('none', 'school-aged', 'young children','everyone')#,'all children','under 5s'

# Rounds of RTSS mass vaccination
RTSSrounds <- c('none','single','every 3 years')

# adding a fifth RTS,S dose: 0, 1
fifth <- c(0)

interventions <- crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, RTSSrounds, fifth, MDAcov, MDAtiming)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse, drawID, sep = "_")) 

# remove non-applicable scenarios 
combo <- combo |>
  filter(!(RTSS=='none' & RTSScov > 0)) |>
  filter(!(RTSS %in% c('SVmass+EPI','SVmass+hybrid','mass+EPI','EPI','hybrid','catch-up') & RTSScov == 0)) |>
  filter(!(RTSScov == 0 & RTSSage %in% c('all children','everyone','school-aged','under 5s','young children','over 5s'))) |>
  filter(!(RTSScov > 0 & RTSSage == 'none')) |>
  filter(!(RTSS %in% c('EPI','hybrid', 'SV') & RTSSage %in% c('all children', 'everyone','school-aged','under 5s','over 5s','none'))) |> # age and hybrid and SV only given to young children 
  filter(!(RTSS %in% c('EPI','hybrid', 'SV') & RTSSrounds %in% c('single','every 3 years'))) |> # they don't get mass vaccination rounds
  filter(!(RTSSage == 'none' & RTSSrounds %in% c('single','every 3 years'))) |>
  filter(!((RTSS =='none' | RTSS == 'EPI') & fifth == 1)) |>
  filter(!(seas_name == "perennial" & (RTSS == "SVmass+EPI" |RTSS =='SVmass+hybrid' | RTSS == "hybrid" | RTSS == 'SV'))) |>
  filter(!((RTSS == 'SVmass+EPI' |RTSS =='SVmass+hybrid'| RTSS =='mass+EPI') & RTSSrounds == 'none')) |>
  filter(!((RTSS == 'SVmass+EPI'|RTSS =='SVmass+hybrid'|RTSS == 'mass+EPI') & RTSSage %in% c('young children'))) |> # when routine vaccination is to children 5-17 months, wouldn't do mass to them too
  filter(!((MDAcov > 0 | MDAtiming == 'during') & RTSS == 'none')) |> # only testing MDA combined with RTSS or RTSS alone, not MDA on its own
  filter(!(MDAcov > 0 & MDAtiming == 'none')) |> # can't have coverage of MDA and no timing of MDA
  filter(!(MDAcov == 0 & MDAtiming == 'during')) # can't have coverage of MDA at 0% and MDA during vax

# put variables into the same order as function arguments
combo <- combo |> 
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         ITN,               # ITN type
         ITNuse,            # ITN usage
         ITNboost,          # if ITN usage is boosted by 10%
         resistance,        # resistance level
         IRS,               # IRS coverage
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         RTSSage,           # RTS,S age groups
         RTSSrounds,        # RTS,S rounds of mass vax
         fifth,             # status of 5th dose for mass or hybrid strategies
         MDAcov,            # MDA coverage
         MDAtiming,         # timing of MDA round
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()

# rearrange so new ones are at the end 
# new <- combo |>
#   filter(RTSSrounds=='every 3 years' | RTSS == 'mass+EPI' | pfpr == 0.03)
# combo <- combo |>
#   filter(!(RTSSrounds=='every 3 years' | RTSS == 'mass+EPI' | pfpr == 0.03)) |>
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

# define all combinations of scenarios and draws
index <- tibble(x = x)
# index <- index[1837:5355,]
# remove ones that have already been run
index <- index |>
  mutate(f = paste0(HPCpath, "HPC/raw_modelrun_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with the first scenario
t <- obj$enqueue_bulk(897, runsim) #936-svmass+hybrid #545
t$status()
# t$wait(1000)
#t$results()

# submit jobs, 10 as a time
sjob <- function(x, y){
  
  t <- obj$enqueue_bulk(index[x:y,], runsim)
  return(1)
  
}

map2_dfr(seq(0, nrow(index)- 100, 100),
         seq(99, nrow(index), 100),
         sjob)


# submit all remaining tasks
t <- obj$enqueue_bulk(index, runsim)
# t$status()
