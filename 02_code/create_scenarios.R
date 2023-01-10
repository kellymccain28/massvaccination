source('02_code/data_libraries.R')

# MODEL set-up ----
# year
year <- 365

# population
population <- 10000

# run time
warmup <- 9 * year       # needs to be multiple of 3 for ITN distribution
sim_length <- 3 * year   # value > 0 

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- c(0)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- seq(.10, .75, .05) # start at 10% (min rec for RTSS use)

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
ITNuse <- c(0, 0.75)#0.25, 0.50, 

# boosting ITN use by 10%: 0, 1
ITNboost <- c(0) 

# insecticide resistance level: 0, 0.4, 0.8
resistance <- c(0) 

# IRS coverage
IRS <-  c(0)

# treatment coverage
treatment <- c(0.45) 

# SMC: 0, 1
SMC <- c(0, 1) 

# RTS,S: none, EPI, SV, hybrid
RTSS <- c('SV') 

# RTS,S coverage
RTSScov <- c(0.8) 

# adding a fifth RTS,S dose: 0, 1
fifth <- c(0)  


interventions <- crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse, drawID, sep = "_")) 

# remove non-applicable scenarios
combo <- combo |>
  filter(!(seas_name == 'highly seasonal' & SMC == 0)) |>
  filter(!(seas_name %in% c('seasonal', 'perennial') & SMC != 0)) |>
  filter(!(RTSS == 'none' & RTSScov == 0.9))

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
         fifth,             # status of 5th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()

saveRDS(combo, paste0(path, '03_output/runscenarios.rds'))

# generate parameter list in malariasimulation format
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/runscenarios.rds'), # file path to pull
                paste0(path, "03_output/run_parameters.rds"))      # file path to push


# Run tasks --------------------------------------------------------------------
x = c(1:nrow(combo)) # runs

index <- crossing(x)

### unsure how to run these***