source('02_code/packages_data.R')

# MODEL set-up ----
# year
year <- 365

# population
population <- 50000

# run time
warmup <- 12 * year       # needs to be multiple of 3 for ITN distribution
sim_length <- 12 * year   # value > 0 

# number of parameter draws
# 0 = use mean values, 1 to 50 = draws
drawID <- c(0)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65) # same profiles as Penny et al. 

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
RTSS <- c('none') 

# RTS,S age group
RTSSage <- c(0)#c('all children','everyone','school-aged','under 5s','young children')

# RTS,S coverage
RTSScov <- c(0) 

# adding a fifth RTS,S dose: 0, 1
fifth <- c(0)  


interventions <- crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, RTSSage, fifth)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse, RTSSage, drawID, sep = "_")) 

# remove non-applicable scenarios -- we are not assuming SMC or RTSS so not applicable
# combo <- combo |>
#   filter(!(seas_name == 'highly seasonal' & SMC == 0)) |>
#   filter(!(seas_name %in% c('perennial') & SMC != 0)) |>
#   filter(!(RTSS == 'none' & RTSScov == 0.9))

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

saveRDS(combo, paste0(path, '03_output/baseline_scenarios.rds'))

# generate parameter list in malariasimulation format -----
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/baseline_scenarios.rds'), # file path to pull
                paste0(path, "03_output/baseline_parameters.rds"))      # file path to push

### Do PfPR/EIR matching -----
source('02_code/Functions/eir_prev_matching.R')

lapply(1:nrow(combo), pr_match)

# Results ----------------------------------------------------------------------
# read in results
files <- list.files(path = paste0("03_output/PrEIR"), pattern = "PRmatch_draws_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("rbind", dat_list) |> as_tibble()

summary(match$starting_EIR)

# take a look at failed jobs
anti_join(combo, match, by = "ID") 


# if needed, use the code to set these IDs to a starting EIR of XX (the upper limit in the PRmatch.R function)
# remainder <- anti_join(combo |> mutate(scenarioID = row_number()), match, by = "ID") |> 
#   mutate(starting_EIR = 1000) |>
#   select(scenarioID, drawID, starting_EIR, ID)
# 
# match <- full_join(match, remainder)


# save EIR estimates, both on local machine and on shared drive
saveRDS(match, paste0(path, "03_output/PrEIR/EIRestimates.rds"))
saveRDS(match, paste0(HPCpath, "03_output/EIRestimates.rds"))


