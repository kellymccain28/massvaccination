# Run scenarios in malariasimulation on HPC
### This script first sets up the HPC. 
### Next, creates a set of scenarios that we want to test and then uses the `generate_params` 
### function to create a list of parameters for each scenario and outputs this as a file. 
### The function `run_simulation` is run for each scenario in this parameter list.

## Setting up cluster ------------------------------------------------------------------------------------------------
source('02_code/data_libraries.R')

setwd(HPCpath)

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "kem22")
didehpc::web_login()

# to edit HPC username and password below
# usethis::edit_r_environ

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, 'Functions/run_simulation.R')),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")
config <- didehpc::didehpc_config(credentials = list(
  username = Sys.getenv("DIDE_USERNAME"),
  password = Sys.getenv("DIDE_PASSWORD")),
  shares = share,
  use_rrq = FALSE,
  cores = 1,
  cluster = "fi--didemrchnb",
  template = "32Core", # "GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core", "32Core"
  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)


## Model set up ------------------------------------------------------------------------------------------------
# year
year <- 365

# population
population <- 10000

# run time
warmup <- 3 * year       # needs to be multiple of 3 for ITN distribution
sim_length <- 6 * year   # value > 0 

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- c(0)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.03, seq(.05, .65, .05)) # start at 10% (min rec for RTSS use)

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

stable <- bind_rows(s1, s2, s3)

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

# treatment coverage
treatment <- c(0.45) 

# SMC: 0, 1
SMC <- c(0) 

# RTS,S: none, EPI, SV, hybrid
RTSS <- c('none','SV') 

# RTS,S coverage
RTSScov <- c(0, 0.8) 

# RTS,S age group
RTSSage <- c('all children','everyone','school-aged','under 5s','young children')

# adding a fifth RTS,S dose: 0, 1
fifth <- c(0, 1)  

interventions <- crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, fifth)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse, drawID, sep = "_")) 

# remove non-applicable scenarios -- we are not assuming SMC or RTSS so not applicable
combo <- combo |>
  filter(!(seas_name == 'highly seasonal' & SMC == 0)) |>
  filter(!(seas_name %in% c('perennial') & SMC != 0)) |>
  filter(!(RTSS == 'none' & RTSScov > 0)) |>
  filter(!(RTSS == 'SV' & RTSScov == 0))

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
         fifth,             # status of 5th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()

saveRDS(combo, paste0(path, '03_output/scenarios_torun.rds'))

# generate parameter list in malariasimulation format
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/scenarios_torun.rds'), # file path to pull
                paste0(path, "03_output/parameters_torun.rds"))      # file path to push

# Run tasks -------------------------------------------------------------------------------------------------------------
x = c(1:nrow(combo)) # runs

# define all combinations of scenarios and draws
index <- tibble(x = x)

# remove ones that have already been run
index <- index |> 
  mutate(f = paste0(HPCpath, "HPC/general_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with the first scenario
# t <- obj$enqueue_bulk(index[1,], runsim)

# submit jobs, 10 as a time
sjob <- function(x, y){
  
  t <- obj$enqueue_bulk(index[x:y,], run_simulation)
  return(1)
  
}

map2_dfr(seq(0, nrow(index) - 10, 10),
         seq(9, nrow(index), 10),
         sjob)


# submit all remaining tasks
# t <- obj$enqueue_bulk(index, runsim)
# t$status()
