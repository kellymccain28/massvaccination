# Generate parameters for malariasimulation ------------------------------------

generate_params <- function(inputpath,   # path to input scenarios
                            outputpath){ # path where output file will be stored
  #setwd(path)
  # read in dataframe of all scenario combinations
  scenarios <- readRDS(inputpath)
  
  # generate parameters
  generate_params2 <- function(x){ # x = scenario number
    
    # pull one scenario at a time
    data <- scenarios[x, ]
    
    # assign values
    population = data$population
    seasonality = data$seasonality
    seas_name = data$seas_name
    pfpr = data$pfpr
    warmup = data$warmup
    sim_length = data$sim_length
    speciesprop = data$speciesprop
    treatment = data$treatment
    SMC = data$SMC
    PEV = data$PEV
    PEVstrategy = data$PEVstrategy
    PEVcov = data$PEVcov
    PEVage = data$PEVage
    PEVrounds = data$PEVrounds
    EPIbooster = data$EPIbooster
    EPIextra = data$EPIextra
    # fifth = data$fifth
    massbooster_rep = data$massbooster_rep
    MDA = data$MDA
    # MDAtiming = data$MDAtiming
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
    params$clinical_incidence_rendering_min_ages = c(c(0, 0.25, seq(1, 20, by = 1))*year, seq(0, 95, by = 5)*year, 0, 5*year, 5*year) #
    params$clinical_incidence_rendering_max_ages = c(c(0.25, seq(1, 21, by = 1))*year, seq(5, 100, by = 5)*year, 100*year, 15*year, 100*year) #
    
    # Set severe incidence rendering 
    params$severe_incidence_rendering_min_ages = c(c(0, 0.25, seq(1, 20, by = 1))*year, seq(0, 95, by = 5)*year, 0, 5*year, 5*year) #
    params$severe_incidence_rendering_max_ages = c(c(0.25, seq(1, 21, by = 1))*year, seq(5, 100, by = 5)*year, 100*year, 15*year, 100*year) #
    
    # Set age group rendering 
    params$age_group_rendering_min_ages = c(c(0, 0.25, seq(1, 20, by = 1))*year, seq(0, 95, by = 5)*year, 0, 5*year, 5*year) #
    params$age_group_rendering_max_ages = c(c(0.25, seq(1, 21, by = 1))*year, seq(5, 100, by = 5)*year, 100*year, 15*year, 100*year) #
    
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
    
    # vectors ----------
    # params <- set_species(
    #   parameters = params,
    #   species = list(arab_params, fun_params, gamb_params),
    #   proportions = unlist(speciesprop))
    # 
    # # proportion of bites taken in bed for each species
    # # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    # params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
    # # proportion of bites taken indoors for each species
    # params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020
    
    # Set initial carrying capacity ----
    # carryingcapacity <- get_init_carrying_capacity(params)
    # 
    # params <- set_carrying_capacity(carrying_capacity = carryingcapacity, 
    #                                 timesteps = )
    
    # ITNs ----------
    # find values in S.I. of 10.1038/s41467-018-07357-w
    # or in Table S1.3 of 10.1016/S2542-5196(21)00296-5
    # same value for all species
    # bednet_timesteps <- c(0)
    # 
    # # no resistance
    # dn0_1 <- 0.387 # pyr, 0 resistance
    # 
    # # resistance
    # dn0_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.387,
    #                    ITN=='pyr' & resistance==0.4 ~ 0.352,
    #                    ITN=='pyr' & resistance==0.8 ~ 0.270,
    #                    ITN=='pbo' & resistance==0 ~ 0.517,
    #                    ITN=='pbo' & resistance==0.4 ~ 0.494,
    #                    ITN=='pbo' & resistance==0.8 ~ 0.419)
    # # no resistance
    # rn_1 <- 0.563 # pyr, 0 resistance
    # 
    # # resistance
    # rn_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.563,
    #                   ITN=='pyr' & resistance==0.4 ~ 0.568,
    #                   ITN=='pyr' & resistance==0.8 ~ 0.626,
    #                   ITN=='pbo' & resistance==0 ~ 0.474,
    #                   ITN=='pbo' & resistance==0.4 ~ 0.493,
    #                   ITN=='pbo' & resistance==0.8 ~ 0.525)
    # 
    # # no resistance
    # gamman_1 <- 2.64 # pyr, 0 resistance
    # 
    # # resistance
    # gamman_2 <- case_when(ITN=='pyr' & resistance==0 ~ 2.64,
    #                       ITN=='pyr' & resistance==0.4 ~ 2.226,
    #                       ITN=='pyr' & resistance==0.8 ~ 1.616,
    #                       ITN=='pbo' & resistance==0 ~ 2.64,
    #                       ITN=='pbo' & resistance==0.4 ~ 2.160,
    #                       ITN=='pbo' & resistance==0.8 ~ 1.311)
    # 
    # ITNuse1 = ITNuse
    # ITNuse2 = ITNuse1 + .10   # ITN boost by 10%
    # 
    # if (ITNboost == 0) {      # if ITNs are not boosted, keep ITN use constant
    #   ITNuse2 = ITNuse1
    # }
    # 
    # npre <- ceiling(warmup / (3 * year))      # number of distributions during warmup
    # npost <- ceiling(sim_length / (3 * year)) # number of distributions during sim_length
    # 
    # params <- set_bednets(
    #   parameters = params,
    #   timesteps = c(
    #     # baseline coverage starts at timestep 1
    #     seq(1, (warmup), 3 * year),  
    #     # intervention coverage starts at sim_length
    #     seq(warmup + 1, (warmup + sim_length), 3 * year)), 
    #   
    #   coverages = c(rep(ITNuse1, npre),      # set baseline coverage
    #                 rep(ITNuse2, npost)),    # set intervention coverage
    #   
    #   retention = 5 * year,
    #   
    #   dn0 = matrix(c(rep(dn0_1, npre), rep(dn0_2, npost),
    #                  rep(dn0_1, npre), rep(dn0_2, npost),
    #                  rep(dn0_1, npre), rep(dn0_2, npost)),
    #                nrow = npre + npost, ncol = 3),
    #   rn = matrix(c(rep(rn_1, npre), rep(rn_2, npost),
    #                 rep(rn_1, npre), rep(rn_2, npost),
    #                 rep(rn_1, npre), rep(rn_2, npost)),
    #               nrow = npre + npost, ncol = 3),
    #   rnm = matrix(c(rep(.24, npre + npost),
    #                  rep(.24, npre + npost),
    #                  rep(.24, npre + npost)),
    #                nrow = npre + npost, ncol = 3),
    #   gamman = c(rep(gamman_1 * year, npre), rep(gamman_2 * year, npost))
    # )
    # 
    # # var for outputting to check ITN timings are correct
    # bednet_timesteps <- params$bednet_timesteps - warmup
    # 
    # 
    # # IRS ----------
    # # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    # if (IRS > 0) {
    #   params <- set_spraying(
    #     parameters = params,
    #     timesteps = seq(1, sim_length, year),
    #     coverages = rep(IRS, sim_length/year),
    #     ls_theta = matrix(c(rep(2.025, sim_length/year),
    #                         rep(2.025, sim_length/year),
    #                         rep(2.025, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3),
    #     ls_gamma = matrix(c(rep(-0.009, sim_length/year),
    #                         rep(-0.009, sim_length/year),
    #                         rep(-0.009, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3),
    #     ks_theta = matrix(c(rep(-2.222, sim_length/year),
    #                         rep(-2.222, sim_length/year),
    #                         rep(-2.222, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3),
    #     ks_gamma = matrix(c(rep(0.008, sim_length/year),
    #                         rep(0.008, sim_length/year),
    #                         rep(0.008, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3),
    #     ms_theta = matrix(c(rep(-1.232, sim_length/year),
    #                         rep(-1.232, sim_length/year),
    #                         rep(-1.232, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3),
    #     ms_gamma = matrix(c(rep(-0.009, sim_length/year),
    #                         rep(-0.009, sim_length/year),
    #                         rep(-0.009, sim_length/year)),
    #                       nrow=sim_length/year, ncol=3)
    #   )  }
    # 
    # # treatment ----------
    if (treatment > 0) {
      params <- set_drugs(
        parameters = params,
        drugs = list(AL_params, SP_AQ_params))
      
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_clinical_treatment(
        parameters = params,
        drug = 1,
        timesteps = c(1),
        coverages = c(treatment)
      )  
    }
    
    
    # SMC ----------
    # smc_timesteps <- 0
    # 
    # if (SMC > 0 & seas_name == 'seasonal') {
    #   peak <- peak_season_offset(params)
    #   # 5 doses, centered around peak
    #   first <- round(warmup + c(peak + c(-2, -1, 0, 1, 2) * month), 0)
    #   firststeps <- sort(rep(first, sim_length/year))
    #   yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
    #   timesteps <- yearsteps + firststeps
    #   
    #   params <- set_drugs(
    #     parameters = params,
    #     list(AL_params, SP_AQ_params))
    #   
    #   params$drug_prophylaxis_scale <- c(10.6, 39.34)
    #   params$drug_prophylaxis_shape <- c(11.3, 3.40)
    #   
    #   params <- set_smc(
    #     parameters = params,
    #     drug = 2,
    #     timesteps = sort(timesteps),
    #     coverages = rep(SMC, length(timesteps)),
    #     min_age = round(0.25 * year),
    #     max_age = round(5 * year))
    #   
    #   smc_timesteps <- params$smc_timesteps - warmup
    # }
    # 
    # if (SMC > 0 & seas_name == 'highly seasonal') {
    #   peak <- peak_season_offset(params)
    #   # 4 doses, centered around peak
    #   first <- round(c(peak + c(-1, 0, 1, 2) * month), 0)
    #   firststeps <- sort(rep(first, (warmup + sim_length)/year))
    #   yearsteps <- rep(c(0, seq(year, (warmup + sim_length) - year, year)), length(first))
    #   timesteps <- yearsteps + firststeps
    #   
    #   params <- set_drugs(
    #     parameters = params,
    #     list(AL_params, SP_AQ_params))
    #   
    #   params$drug_prophylaxis_scale <- c(10.6, 39.34)
    #   params$drug_prophylaxis_shape <- c(11.3, 3.40)
    #   
    #   params <- set_smc(
    #     parameters = params,
    #     drug = 2,
    #     timesteps = sort(timesteps),
    #     coverages = rep(SMC, length(timesteps)),
    #     min_ages = rep((0.25 * year), length(timesteps)),
    #     max_ages = rep((5 * year), length(timesteps)))
    #   
    #   # var for outputting to check SMC timings are correct
    #   smc_timesteps <- params$smc_timesteps - warmup
    # }
    
    # Vaccination with a PEV ----
    if (PEV == 'RTSS'){
      # RTSS ----------
      if (PEVcov > 0) {
        program_start <- 5 * year
        
        if (PEVage == '5-9'){
          min_ages = 5 * year
          max_ages = 9 * year
        } else if (PEVage == '5-15'){
          min_ages = 5 * year
          max_ages = 15 * year
        } else if (PEVage == '5-100'){
          min_ages = 5 * year
          max_ages = 100 * year
        }
        
        # AB 
        if (PEVstrategy == "AB") {
          params$pev_doses <- round(c(0, 1.5 * month, 3 * month)) # spacing from Hillary's RTSS work 
          # params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
          
          if(EPIbooster == '12m' | EPIbooster == '12m boost') {
            epiboosters <- round(12 * month)
          } else if (EPIbooster == '18m') {
            epiboosters <- round(18 * month) # from phase III R21 trial
          } 
          
          boost_cov <- 0.8
          
          pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
          
          if (EPIbooster == '12m boost'){
            # update booster to have same effect as dose 3 per Thompson et al. 2022 (when time between 3rd and 4th dose is 12 mo)
            rtss_booster_profile$cs <- c(6.37008, 0.35)
          }
          
          params <- set_pev_epi(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps, 
            coverages = PEVcov,
            age = round(5 * month),
            min_wait = 0,
            booster_timestep = epiboosters,
            booster_coverage = boost_cov,
            booster_profile = list(rtss_booster_profile),
            seasonal_boosters = FALSE
          )  
          
          print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
          print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
          
        }
        
        # hybrid ----------
        if (PEVstrategy == "hybrid") {
          params$pev_doses <- round(c(0, 1.5 * month, 3 * month)) # spacing from the RTSS work 
          # params$pev_doses <- round(c(0, 1 * month, 2 * month)) # spacing from phase iii trial R21
          
          peak <- peak_season_offset(params)
          
          hybridbooster <- round(c(peak - 3.5 * month), 0)
          
          boost_cov <- 0.8  # coverage from 10.1016/S2214-109X(22)00416-8 #else c(PEVcov*0.8, PEVcov*0.8*0.9)
          pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
          
          # update booster to have same effect as dose 3 per Thompson et al. 2022 (when time between 3rd and 4th dose is 12 mo)
          rtss_booster_profile$cs <- c(6.37008, 0.35)
          
          params <- set_pev_epi(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps,
            coverages = PEVcov,
            age = round(5 * month),
            min_wait = 0,
            booster_timestep = hybridbooster,
            booster_coverage = boost_cov,
            booster_profile = list(rtss_booster_profile),
            seasonal_boosters = TRUE
          ) 
          
          print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
          print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
        }
        
        # mass ----------
        if (PEVstrategy == 'mass'){
          
          params$pev_doses <- round(c(0, 1.5 * month, 3 * month)) # spacing from the RTSS work 
          # params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
          
          # First set the EPI strategy 
          EPIboosters <- round(c(12 * month))  # from phase III R21 trial
          pevtimesteps <- warmup + program_start # starting when warmup ends + program start
          EPIboost_cov <- 0.8 
          
          # update 4th booster to have same effect as dose 3 per Thompson et al. 2022 (when time between 3rd and 4th dose is 12 mo)
          rtss_booster_4th <- rtss_booster_profile
          rtss_booster_4th$cs <- c(6.37008, 0.35)
          
          params <- set_pev_epi(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps,
            coverages = PEVcov,
            age = round(5 * month),
            min_wait = 0,
            booster_timestep = EPIboosters,
            booster_coverage = EPIboost_cov,
            booster_profile = list(rtss_booster_4th),
            seasonal_boosters = FALSE)
          
          # Get timing for mass vaccination rounds
          peak <- peak_season_offset(params)
          
          if (massbooster_rep == '-') {
            massboosters <- round(12 * month) # 1 year after 3rd dose 
          } else if (massbooster_rep == '4 annual' | massbooster_rep == '4 annual, no drop') {
            massboosters <- round(seq(12, 48) * month) # 4 annual boosters after 3rd dose
          } else if (massbooster_rep == 'annual' | massbooster_rep == 'annual no drop'){
            massboosters <- round(seq(12, (sim_length - program_start)/month, by = 12) * month)
          }
          
          if (massbooster_rep %in% c('-', '4 annual', 'annual')){
            massbooster_cov <- c(0.8, rep(0.9, length(massboosters)-1)) # coverage from 10.1016/S2214-109X(22)00416-8
          } else if (massbooster_rep %in% c('4 annual, no drop', 'annual no drop')){
            massbooster_cov <- c(0.8, rep(1, length(massboosters)-1))
          }
          
          if(seas_name == 'seasonal'){
            first <- round(warmup + program_start + (peak - month * 3.5), 0) # after program start, 3.5 months prior to peak for seasonal (RTSS-CE repo)
          } else if(seas_name == 'perennial'){
            first <- round(warmup + program_start) # non seasonally distributed mass campaign
          }
          
          if (PEVrounds == 'single'){
            pevtimesteps <- c(first)
            min_wait <- 0
          } else if (PEVrounds == '1yr'){
            pevtimesteps <- c(first, first + seq(1 * year, sim_length, 1 * year))
            min_wait <- 20 * year
          } else if (PEVrounds == '3yrs'){
            pevtimesteps <- c(first, first + seq(3 * year, sim_length, 3 * year))
            min_wait <- 20 * year
          } else if (PEVrounds == '5yrs'){
            pevtimesteps <- c(first, first + seq(5 * year, sim_length, 5 * year))
            min_wait <- 20 * year
          }
          
          # Set mass vaccination 
          proppop_notpregnant <- 1 - 0.078/2 # from DHS data - see Get_pregnancy_rate.R
          
          # add mass vaccination taking into account pregnant women
          params <- set_mass_pev(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps, 
            coverages = rep(PEVcov * proppop_notpregnant, length(pevtimesteps)), 
            min_ages = min_ages,
            max_ages = max_ages,
            min_wait = min_wait,
            booster_timestep = massboosters, # timesteps following initial vaccination 
            booster_profile = append(replicate(length(massboosters)-1, rtss_booster_profile, simplify = FALSE), list(rtss_booster_4th), after = 0),
            booster_coverage = massbooster_cov # prop of vaccinated pop who will receive booster vaccine
          )
          
          
          # var for outputting to check RTS,S timings are correct
          print(paste0("Mass timesteps: ", mass_pev_timesteps <- params$mass_pev_timesteps - warmup))
          print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
          print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
          print(paste0('Mass booster: ', params$mass_pev_booster_timestep))
        }
        
        if (PEVstrategy == 'catch-up'){
          
          # Get timing for EPI boosters in catch-up campaigns
          params$pev_doses <- round(c(0, 1.5 * month, 3 * month)) # spacing from the RTSS work 
          # params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
          
          if (EPIbooster == '12m boost' & EPIextra == '5y') {
            epiboosters <- round(c(12 * month, 5 * year))
          } else if (EPIbooster == '12m boost' & EPIextra == '10y') {
            epiboosters <- round(c(12 * month, 10 * year))
          } else if (EPIbooster == '12m boost' & EPIextra == '-') {
            epiboosters <- round(12 * month)
          }
          
          epiboost_cov <- c(0.8, rep(0.9, length(epiboosters)-1))#0.8 * (0.9 ^ (0:(length(epiboosters)-1)))
          
          massboosters <- round(c(1 * year)) # 1 year after 3rd dose 
          
          massboost_cov <- 0.8
          
          pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
          
          # update 4th booster to have same effect as dose 3 per Thompson et al. 2022 (when time between 3rd and 4th dose is 12 mo)
          rtss_booster_4th <- rtss_booster_profile
          rtss_booster_4th$cs <- c(6.37008, 0.35)
          # any boosters afterward would be normal boosters
          epiboosterprofiles <- if (length(epiboosters) == 1){ 
            list(rtss_booster_4th) 
          } else if (length(epiboosters) == 2){
            list(rtss_booster_4th, rtss_booster_profile)}
          
          # Set EPI strategy for young children 
          params <- set_pev_epi(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps, 
            coverages = PEVcov,
            age = round(5 * month),
            min_wait = 0,
            booster_timestep = epiboosters,
            booster_coverage = epiboost_cov,
            booster_profile = rep(list(rtss_booster_profile), length(epiboosters)),
            seasonal_boosters = FALSE
          )  
          
          # Set catch-up mass campaigns
          params <- set_mass_pev(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps,
            coverages = rep(PEVcov, length(pevtimesteps)),
            min_ages = min_ages,
            max_ages = max_ages,
            min_wait = 0,
            booster_timestep = massboosters,
            booster_coverage = massboost_cov,
            booster_profile = rep(list(rtss_booster_profile), length(massboosters))
          )
          
          print(paste0("Mass timesteps: ", mass_pev_timesteps <- params$mass_pev_timesteps - warmup))
          print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
          print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
          print(paste0('Mass booster: ', params$mass_pev_booster_timestep))
        }
      }
    } else if (PEV == 'R21'){
    
      # R21 ----------
      # r21_efficacy <- read.csv(paste0(path, '01_data/efficacy_parameters.csv')) %>%
      #   summarise(v_max = median(v_max),
      #             alpha = median(alpha),
      #             beta = median(beta))
      # r21_params <- read.csv(paste0(path,'01_data/r21_malarisimulation_parameter_draws.csv')) %>%
      #   group_by(par) %>%
      #   summarize(median_mu = median(mu),
      #             median_sd = median(sd))
      if (PEVcov >0){
          r21_profile <- rtss_profile
          r21_profile$vmax <- 0.84#0.8441944
          r21_profile$alpha <- 0.99#0.9879701
          r21_profile$beta <- 455.69#455.6919
          r21_profile$cs <- c(9.32, 0.839)
          r21_profile$rho <- c(0.791, 0.607)
          r21_profile$ds <- c(3.79, 0.165)
          r21_profile$dl <- c(6.28, 0.433)
          
          r21_booster_profile <- r21_profile
          r21_booster_profile$rho <- c(0.0821, 0.547)
          r21_booster_profile$cs <- c(9.24, 0.719)
          
          r21_booster_profile2 <- r21_booster_profile
          r21_booster_profile2$cs <- c(9.02, 0.845)
          
          program_start <- 5 * year
          
          if (PEVage == '5-9'){
            min_ages = 5 * year
            max_ages = 9 * year
          } else if (PEVage == '5-15'){
            min_ages = 5 * year
            max_ages = 15 * year
          } else if (PEVage == '5-100'){
            min_ages = 5 * year
            max_ages = 100 * year
          }
          
          # AB
          if (PEVstrategy == "AB") {
            params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
            
            epiboosters <- round(12 * month)
            
            boost_cov <- 0.8
            
            pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
            
            params <- set_pev_epi(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = PEVcov,
              age = round(5 * month),
              min_wait = 0,
              booster_timestep = epiboosters,
              booster_coverage = boost_cov,
              booster_profile = list(r21_booster_profile),
              seasonal_boosters = FALSE
            )
            
            print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
            print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
            
          }
          
          # hybrid ----------
          if (PEVstrategy == "hybrid") {
            params$pev_doses <- round(c(0, 1 * month, 2 * month)) # spacing from phase iii trial R21
            
            peak <- peak_season_offset(params)
            
            hybridbooster <- round(c(peak - 3.5 * month), 0)
            
            boost_cov <- 0.8  # coverage from 10.1016/S2214-109X(22)00416-8 #else c(R21cov*0.8, R21cov*0.8*0.9)
            pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
            
            params <- set_pev_epi(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = PEVcov,
              age = round(5 * month),
              min_wait = 0,
              booster_timestep = hybridbooster,
              booster_coverage = boost_cov,
              booster_profile = list(r21_booster_profile),
              seasonal_boosters = TRUE
            )
            
            print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
            print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
            print('hybrid parameterized)')
          }
          
          # mass ----------
          if (PEVstrategy == 'mass'){
            
            params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
            
            # First set the EPI strategy
            EPIboosters <- round(c(12 * month))  # from phase III R21 trial
            pevtimesteps <- warmup + program_start # starting when warmup ends + program start
            EPIboost_cov <- 0.8
            
            params <- set_pev_epi(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = PEVcov,
              age = round(5 * month),
              min_wait = 0,
              booster_timestep = EPIboosters,
              booster_coverage = EPIboost_cov,
              booster_profile = list(r21_booster_profile),
              seasonal_boosters = FALSE)
            
            # Get timing for mass vaccination rounds
            peak <- peak_season_offset(params)
            
            if (massbooster_rep == '-') {
              massboosters <- round(12 * month) # 1 year after 3rd dose 
            } else if (massbooster_rep == '4 annual' | massbooster_rep == '4 annual, no drop') {
              massboosters <- round(seq(12, 48) * month) # 4 annual boosters after 3rd dose
            } else if (massbooster_rep == 'annual' | massbooster_rep == 'annual no drop'){
              massboosters <- round(seq(12, (sim_length - program_start)/month, by = 12) * month)
            }
            
            if (massbooster_rep %in% c('-', '4 annual', 'annual')){
              massbooster_cov <- c(0.8, rep(0.9, length(massboosters)-1)) # coverage from 10.1016/S2214-109X(22)00416-8
            } else if (massbooster_rep %in% c('4 annual, no drop', 'annual no drop')){
              massbooster_cov <- c(0.8, rep(1, length(massboosters)-1))
            }
            
            if (length(massboosters) == 1){ 
              massboosterprofiles <- list(r21_booster_profile) 
            } else if (length(massboosters) >=4){
              massboosterprofiles <- append(replicate(length(massboosters)-1, r21_booster_profile2, simplify = FALSE), list(r21_booster_profile), after = 0)
            } 
            
            if(seas_name == 'seasonal'){
              first <- round(warmup + program_start + (peak - month * 3.5), 0) # after program start, 3.5 months prior to peak for seasonal (R21-CE repo)
            } else if(seas_name == 'perennial'){
              first <- round(warmup + program_start) # non seasonally distributed mass campaign
            }
            
            if (PEVrounds == 'single'){
              pevtimesteps <- c(first)
              min_wait <- 0
            } else if (PEVrounds == '1yr'){
              pevtimesteps <- c(first, first + seq(1 * year, sim_length, 1 * year))
              min_wait <- 20 * year
            } else if (PEVrounds == '3yrs'){
              pevtimesteps <- c(first, first + seq(3 * year, sim_length, 3 * year))
              min_wait <- 20 * year
            } else if (PEVrounds == '5yrs'){
              pevtimesteps <- c(first, first + seq(5 * year, sim_length, 5 * year))
              min_wait <- 20 * year
            }
            
            # Set mass vaccination
            proppop_notpregnant <- 1 - 0.078/2 # from DHS data - see Get_pregnancy_rate.R
            
            # add mass vaccination taking into account pregnant women
            params <- set_mass_pev(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = rep(PEVcov * proppop_notpregnant, length(pevtimesteps)),
              min_ages = min_ages,
              max_ages = max_ages,
              min_wait = min_wait,
              booster_timestep = massboosters, # timesteps following initial vaccination
              booster_profile = massboosterprofiles,#list(r21_booster_profile, rep(r21_booster_profile2, length(massboosters) -1)),
              booster_coverage = massbooster_cov # prop of vaccinated pop who will receive booster vaccine
            )
            
            # var for outputting to check RTS,S timings are correct
            print(paste0("Mass timesteps: ", mass_pev_timesteps <- params$mass_pev_timesteps - warmup))
            print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
            print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
            print(paste0('Mass booster: ', params$mass_pev_booster_timestep))
          }
          
          if (PEVstrategy == 'catch-up'){
            
            # Get timing for EPI boosters in catch-up campaigns
            params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
            
            if (EPIextra == '5y') {
              epiboosters <- round(c(12 * month, 5 * year))
            } else if (EPIextra == '10y') {
              epiboosters <- round(c(12 * month, 10 * year))
            } else if (EPIextra == '-') {
              epiboosters <- round(12 * month)
            }
            
            epiboost_cov <- c(0.8, rep(0.9, length(epiboosters)-1))#0.8 * (0.9 ^ (0:(length(epiboosters)-1)))
            
            massboosters <- round(c(1 * year)) # 1 year after 3rd dose
            
            massboost_cov <- 0.8
            
            pevtimesteps <- warmup + program_start # starting 5 years after warmup ends
            
            epiboosterprofiles <- if (length(epiboosters) == 1){ 
              list(r21_booster_profile) 
            } else if (length(epiboosters) == 2){
              list(r21_booster_profile, r21_booster_profile2)}
            
            # Set EPI strategy for young children
            params <- set_pev_epi(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = PEVcov,
              age = round(5 * month),# Hillary did 6 months of age 
              min_wait = 0,
              booster_timestep = epiboosters,
              booster_coverage = epiboost_cov,
              booster_profile = epiboosterprofiles, # first booster is one thing, then any others are different
              seasonal_boosters = FALSE
            )
            
            # Set catch-up mass campaigns
            params <- set_mass_pev(
              parameters = params,
              profile = r21_profile,
              timesteps = pevtimesteps,
              coverages = rep(PEVcov, length(pevtimesteps)),
              min_ages = min_ages,
              max_ages = max_ages,
              min_wait = 0,
              booster_timestep = massboosters,
              booster_coverage = massboost_cov,
              booster_profile = list(r21_booster_profile) 
            )
            
            print(paste0("Mass timesteps: ", mass_pev_timesteps <- params$mass_pev_timesteps - warmup))
            print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
            print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
            print(paste0('Mass booster: ', params$mass_pev_booster_timestep))
          }
      }
    }
    
    # MDA  ----------
    if (MDA > 0) { # this will only run if there is MDA, which will only be the case if there is mass vaccination 
      MDAcov <- 0.8
      program_start <- 5 * year
      
      if(PEV !='none'){
        mdatimesteps <- pevtimesteps # round of MDA occurs at same time as first dose
      } else if(PEV == 'none'){
        mdatimesteps <- warmup + program_start
      }
      
      params <- set_drugs(
        parameters = params, 
        list(AL_params, SP_AQ_params)
      )
      
      proppop_notpregnant <- 1 - 0.078/2 # from DHS data - see Get_pregnancy_rate.R
      
      # https://www.who.int/publications/i/item/9789241513104 - for drug, min_ages, and coverage (pregnancy)
      params <- set_mda(
        parameters = params,
        drug = 1, # AL for MDA 
        timesteps = mdatimesteps,
        coverages = rep(MDAcov * proppop_notpregnant, length(mdatimesteps)), # excluding pregnant women
        min_ages = rep(6 * month, length(mdatimesteps)), # starting from 6 months of age 
        max_ages = rep(100 * 365, length(mdatimesteps))
      ) 
    }
    
    # synergy SMC & RTS,S ----------
    # if (SMC > 0 & RTSS %in% c("EPI", "mass", "hybrid")) {
    #   
    #   params$rtss_beta <- 70.9
    #   params$rtss_alpha <- 0.868
    #   params$rtss_vmax <- 0.843
    #   params$rtss_cs_boost <- c(6.37008, 0.35)
    #   
    #   params$drug_prophylaxis_scale <- c(10.6, 45.76)
    #   params$drug_prophylaxis_shape <- c(11.3, 2.87)
    #   
    # }
  

    # correlate interventions  ----------
    # correlations <- get_correlation_parameters(params)
    
    # if (RTSScov == 0.77 & pfpr == 0.40) {
    #   correlations$inter_intervention_rho('rtss', 'bednets', 0.04)
    # }
    #
    # if (RTSScov == 0.72 & pfpr == 0.18) {
    #   correlations$inter_intervention_rho('rtss', 'bednets', 0.07)
    # }
    
    # parameter draws  ----------
    # choose a parameter draw
    
    if(drawID > 0){
      params <- set_parameter_draw(params, drawID)
    }
    
    # save as data.frame
    data$params <- list(params)
    data$scenarioID <- x
    
    # print count
    print(paste(x,'pfpr=',data$pfpr))
    
    return(data)
}
  
  # loop through function to generate parameters one by one
  output <- map_dfr(1:nrow(scenarios), generate_params2)
  
  # save output ----------
  saveRDS(output, outputpath)
  
}
