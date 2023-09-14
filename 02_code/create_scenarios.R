source(paste0(HPCpath,'02_code/packages_data.R'))

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
drawID <- c(942,  40, 541, 497, 877, 697, 400, 450, 806, #0, 
            600, 670, 363, 838, 478, 403, 375, 335, 598, 142,
            919, 444, 986, 659,  71, 457, 891, 188, 432, 975, 
            488, 867, 538, 912, 534, 215, 540, 866, 613, 973, 
            917, 937, 931, 296, 835, 328, 147, 701, 889, 708, 888)# c(0, 1:50)

# SITE set-up ----
# parasite prevalence 2-10 year olds
pfpr <- c(0.01, 0.03, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65) # same profiles as Penny et al. 

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

stable <- bind_rows(s2, s3)#s1, 

# vectors
# list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

# INTERVENTIONS ----
# treatment coverage
treatment <- c(0.45) 

# SMC: 0, 1
SMC <- c(0) 

# MDA y/n
MDA <- c(0)

# Which vaccine
PEV <- c('none') 

# Vaccination strategy
PEVstrategy <- c('none')

# Age group
PEVage <- c('-')

# PEV coverage
PEVcov <- c(0) 

# Rounds of PEV mass vaccination
PEVrounds <- c('-')

# EPI booster timing 
EPIbooster  <- c('-')  

# Extra boosters for EPI 
EPIextra <- c('-')
  
# adding vaccine boosters: 0 - no fifth dose, fifth (2 booster doses), annual (every year), 6mo (every 6 months), 2yrs (every 2 years)
massbooster_rep <- c('-')

interventions <- crossing(treatment, SMC, PEV, PEVstrategy, PEVcov, PEVage, PEVrounds, EPIbooster, EPIextra, massbooster_rep, MDA)

# create combination of all runs 
combo <- crossing(population, pfpr, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, "0", drawID, sep = "_")) #, RTSSage

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
  ) |> as.data.frame()

saveRDS(combo, paste0(path, '03_output/baseline_scenarios.rds'))

# generate parameter list in malariasimulation format -----
source(paste0(path, '02_code/Functions/generate_params.R'))

generate_params(paste0(path, '03_output/baseline_scenarios.rds'), # file path to pull
                paste0(HPCpath, "03_output/baseline_parameters.rds"))      # file path to push

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
  cluster = "fi--dideclusthn", #"fi--didemrchnb",
  template = '8Core',#"32Core","GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core", "32Core"
  parallel = FALSE) 

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation@dev", "github::mrc-ide/cali"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, '02_code/Functions/eir_prev_matching.R')),
                             packages = c("dplyr", "malariasimulation","cali"),
                             package_sources = src)

obj <- didehpc::queue_didehpc(ctx, config = config)


# Run tasks -------------------------------------------------------------------------------------------------------------
x <- c(1:nrow(combo)) # baseline scenarios

# define all combinations of scenarios and draws
index <- tibble(x=x, y = combo$drawID)#, ID = combo$ID)

# remove ones that have already been run
index <- index |>
  mutate(f = paste0(HPCpath, "03_output/PrEIR/PRmatch_draws_", index$ID, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist, -ID)

# run a test with the first scenario
# t <- obj$enqueue_bulk(index[1530,], pr_match)
# t$wait(1000)

# submit jobs, 30 as a time
sjob <- function(x, y){
  t <- obj$enqueue_bulk(index[x:y,], pr_match)
  print(paste0(x, ' to ', y))
}

map2_dfr(seq(1, nrow(index), 100),
         seq(100, nrow(index), 100),
         sjob)
# map2_dfr(900,918, sjob)
### Do PfPR/EIR matching without cluster -----
# source('02_code/Functions/eir_prev_matching.R')
# 
# lapply(1:nrow(combo), pr_match)

# Results ----------------------------------------------------------------------
# read in results
files <- list.files(path = paste0("03_output/PrEIR"), pattern = "PRmatch_draws_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("rbind", dat_list) |> as_tibble()

summary(match$starting_EIR)

# take a look at failed jobs
failed <- anti_join(combo, match, by = "ID") 


# if needed, use the code to set these IDs to a starting EIR of XX (the upper limit in the PRmatch.R function)
# remainder <- anti_join(combo |> mutate(scenarioID = row_number()), match, by = "ID") |>
#   mutate(starting_EIR = 1000) |>
#   select(scenarioID, drawID, starting_EIR, ID)
# 
# match <- full_join(match, remainder)


# save EIR estimates, both on local machine and on shared drive
saveRDS(match, paste0(path, "03_output/PrEIR/EIRestimates.rds"))
saveRDS(match, paste0(HPCpath, "03_output/PrEIR/EIRestimates.rds"))


