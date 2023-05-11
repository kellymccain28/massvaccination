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
    ITN = data$ITN
    ITNuse = data$ITNuse
    ITNboost = data$ITNboost
    resistance = data$resistance
    IRS = data$IRS
    treatment = data$treatment
    SMC = data$SMC
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
    params$clinical_incidence_rendering_min_ages = c(c(0, 0.25, seq(1, 20, by = 1))*year, seq(0, 95, by = 5)*year, 0)
    params$clinical_incidence_rendering_max_ages = c(c(0.25, seq(1, 21, by = 1))*year, seq(5, 100, by = 5)*year, 100*year) 
    
    # Set severe incidence rendering 
    params$severe_incidence_rendering_min_ages = c(c(0, 0.25, seq(1, 20, by = 1))*year, seq(0, 95, by = 5)*year, 0) 
    params$severe_incidence_rendering_max_ages = c(c(0.25, seq(1, 21, by = 1))*year, seq(5, 100, by = 5)*year, 100*year) 
    
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
    params <- set_species(
      parameters = params,
      species = list(arab_params, fun_params, gamb_params),
      proportions = unlist(speciesprop))
    
    # proportion of bites taken in bed for each species
    # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
    # proportion of bites taken indoors for each species
    params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020
    
    # ITNs ----------
    # find values in S.I. of 10.1038/s41467-018-07357-w
    # or in Table S1.3 of 10.1016/S2542-5196(21)00296-5
    # same value for all species
    bednet_timesteps <- c(0)
    
    # no resistance
    dn0_1 <- 0.387 # pyr, 0 resistance
    
    # resistance
    dn0_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.387,
                       ITN=='pyr' & resistance==0.4 ~ 0.352,
                       ITN=='pyr' & resistance==0.8 ~ 0.270,
                       ITN=='pbo' & resistance==0 ~ 0.517,
                       ITN=='pbo' & resistance==0.4 ~ 0.494,
                       ITN=='pbo' & resistance==0.8 ~ 0.419)
    # no resistance
    rn_1 <- 0.563 # pyr, 0 resistance
    
    # resistance
    rn_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.563,
                      ITN=='pyr' & resistance==0.4 ~ 0.568,
                      ITN=='pyr' & resistance==0.8 ~ 0.626,
                      ITN=='pbo' & resistance==0 ~ 0.474,
                      ITN=='pbo' & resistance==0.4 ~ 0.493,
                      ITN=='pbo' & resistance==0.8 ~ 0.525)
    
    # no resistance
    gamman_1 <- 2.64 # pyr, 0 resistance
    
    # resistance
    gamman_2 <- case_when(ITN=='pyr' & resistance==0 ~ 2.64,
                          ITN=='pyr' & resistance==0.4 ~ 2.226,
                          ITN=='pyr' & resistance==0.8 ~ 1.616,
                          ITN=='pbo' & resistance==0 ~ 2.64,
                          ITN=='pbo' & resistance==0.4 ~ 2.160,
                          ITN=='pbo' & resistance==0.8 ~ 1.311)
    
    ITNuse1 = ITNuse
    ITNuse2 = ITNuse1 + .10   # ITN boost by 10%
    
    if (ITNboost == 0) {      # if ITNs are not boosted, keep ITN use constant
      ITNuse2 = ITNuse1
    }
    
    npre <- ceiling(warmup / (3 * year))      # number of distributions during warmup
    npost <- ceiling(sim_length / (3 * year)) # number of distributions during sim_length
    
    params <- set_bednets(
      parameters = params,
      timesteps = c(
        # baseline coverage starts at timestep 1
        seq(1, (warmup), 3 * year),  
        # intervention coverage starts at sim_length
        seq(warmup + 1, (warmup + sim_length), 3 * year)), 
      
      coverages = c(rep(ITNuse1, npre),      # set baseline coverage
                    rep(ITNuse2, npost)),    # set intervention coverage
      
      retention = 5 * year,
      
      dn0 = matrix(c(rep(dn0_1, npre), rep(dn0_2, npost),
                     rep(dn0_1, npre), rep(dn0_2, npost),
                     rep(dn0_1, npre), rep(dn0_2, npost)),
                   nrow = npre + npost, ncol = 3),
      rn = matrix(c(rep(rn_1, npre), rep(rn_2, npost),
                    rep(rn_1, npre), rep(rn_2, npost),
                    rep(rn_1, npre), rep(rn_2, npost)),
                  nrow = npre + npost, ncol = 3),
      rnm = matrix(c(rep(.24, npre + npost),
                     rep(.24, npre + npost),
                     rep(.24, npre + npost)),
                   nrow = npre + npost, ncol = 3),
      gamman = c(rep(gamman_1 * year, npre), rep(gamman_2 * year, npost))
    )
    
    # var for outputting to check ITN timings are correct
    bednet_timesteps <- params$bednet_timesteps - warmup
    
    
    # IRS ----------
    # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    if (IRS > 0) {
      params <- set_spraying(
        parameters = params,
        timesteps = seq(1, sim_length, year),
        coverages = rep(IRS, sim_length/year),
        ls_theta = matrix(c(rep(2.025, sim_length/year),
                            rep(2.025, sim_length/year),
                            rep(2.025, sim_length/year)),
                          nrow=sim_length/year, ncol=3),
        ls_gamma = matrix(c(rep(-0.009, sim_length/year),
                            rep(-0.009, sim_length/year),
                            rep(-0.009, sim_length/year)),
                          nrow=sim_length/year, ncol=3),
        ks_theta = matrix(c(rep(-2.222, sim_length/year),
                            rep(-2.222, sim_length/year),
                            rep(-2.222, sim_length/year)),
                          nrow=sim_length/year, ncol=3),
        ks_gamma = matrix(c(rep(0.008, sim_length/year),
                            rep(0.008, sim_length/year),
                            rep(0.008, sim_length/year)),
                          nrow=sim_length/year, ncol=3),
        ms_theta = matrix(c(rep(-1.232, sim_length/year),
                            rep(-1.232, sim_length/year),
                            rep(-1.232, sim_length/year)),
                          nrow=sim_length/year, ncol=3),
        ms_gamma = matrix(c(rep(-0.009, sim_length/year),
                            rep(-0.009, sim_length/year),
                            rep(-0.009, sim_length/year)),
                          nrow=sim_length/year, ncol=3)
      )  }
    
    # treatment ----------
    if (treatment > 0) {
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_clinical_treatment(
        parameters = params,
        drug = 1,
        timesteps = c(1),
        coverages = c(treatment)
      )  }
    
    
    # SMC ----------
    smc_timesteps <- 0
    
    if (SMC > 0 & seas_name == 'seasonal') {
      peak <- peak_season_offset(params)
      # 5 doses, centered around peak
      first <- round(warmup + c(peak + c(-2, -1, 0, 1, 2) * month), 0)
      firststeps <- sort(rep(first, sim_length/year))
      yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
      timesteps <- yearsteps + firststeps
      
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_smc(
        parameters = params,
        drug = 2,
        timesteps = sort(timesteps),
        coverages = rep(SMC, length(timesteps)),
        min_age = round(0.25 * year),
        max_age = round(5 * year))
      
      smc_timesteps <- params$smc_timesteps - warmup
    }
    
    if (SMC > 0 & seas_name == 'highly seasonal') {
      peak <- peak_season_offset(params)
      # 4 doses, centered around peak
      first <- round(c(peak + c(-1, 0, 1, 2) * month), 0)
      firststeps <- sort(rep(first, (warmup + sim_length)/year))
      yearsteps <- rep(c(0, seq(year, (warmup + sim_length) - year, year)), length(first))
      timesteps <- yearsteps + firststeps
      
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))
      
      params$drug_prophylaxis_scale <- c(10.6, 39.34)
      params$drug_prophylaxis_shape <- c(11.3, 3.40)
      
      params <- set_smc(
        parameters = params,
        drug = 2,
        timesteps = sort(timesteps),
        coverages = rep(SMC, length(timesteps)),
        min_ages = rep((0.25 * year), length(timesteps)),
        max_ages = rep((5 * year), length(timesteps)))
      
      # var for outputting to check SMC timings are correct
      smc_timesteps <- params$smc_timesteps - warmup
    }
    
    # RTSS ----------
    if (RTSScov > 0) {
      program_start <- 1 * year
      
      boost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) # coverage from 10.1016/S2214-109X(22)00416-8
      if (RTSSage == 'young children'){
        min_ages = round(5*month)
        max_ages = round(17*month)
      } else if (RTSSage == 'all children'){
        min_ages = round(5*month)
        max_ages = 15*year
      } else if (RTSSage == 'under 5s'){
        min_ages = round(5*month)
        max_ages = 5*year
      } else if (RTSSage == 'school-aged'){
        min_ages = 5*year
        max_ages = 15*year
      } else if (RTSSage == 'everyone'){
        min_ages = 5*year
        max_ages = 100*year
      }
      
      # EPI 
      if (RTSS == "EPI") {
        params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
        epiboosters <- if(fifth == 0) round(c(12 * month)) else round(c(12 * month, 24 * month)) # from phase III R21 trial
        pevtimesteps <- warmup + program_start # starting when warmup ends
        
        params <- set_pev_epi(
          parameters = params,
          profile = rtss_profile,
          timesteps = pevtimesteps, 
          coverages = RTSScov,
          age = round(5 * month),
          min_wait = 0,
          booster_timestep = epiboosters,
          booster_coverage = boost_cov,
          booster_profile = list(rtss_booster_profile),
          seasonal_boosters = FALSE
        )  
      }
    
      # mass ----------
      
      if (RTSS == "SVmass+EPI" | RTSS == 'SVmass+hybrid' | RTSS == 'mass+EPI') {
        
        params$pev_doses <- round(c(0, 1 * month, 2 * month)) # monthly spacing from phase III R21 trial
        peak <- peak_season_offset(params)
        massboosters <- if(booster_rep == 0){ # only fourth dose
          round(c(1 * year + 2 * month)) 
          } else if (booster_rep == 'fifth'){ # fifth dose
            round(c(1 * year + 2 * month, 2 * year + 2 * month))
          } else if (booster_rep == 'annual'){ # annual doses 
            round(c(1 * year + 2 * month, seq(2 * year + 2 * month, (sim_length / year) * year + 2 * month, by = year)))
          } else if (booster_rep == '2yr'){ # doses every 2 years
            round(c(1 * year + 2 * month, seq(2 * year + 2 * month, (sim_length / year) * year + 2 * month, by = 2 * year)))
          }
        boost_cov <- c(RTSScov * 0.8, rep(RTSScov*0.8*0.9, length(massboosters)-1)) # coverage from 10.1016/S2214-109X(22)00416-8
        # if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) 
        
        if(RTSS == 'SVmass+EPI' | RTSS == "mass+EPI"){
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
            min_wait = 0,
            booster_timestep = EPIboosters,
            booster_coverage = EPIboost_cov,
            booster_profile = list(rtss_booster_profile),
            seasonal_boosters = FALSE)
        } else if (RTSS == 'SVmass+hybrid'){
          # Set the hybrid strategy
          hybridboosters <- if(fifth == 0) round(c(peak - 0.5 * month), 0) else round(((peak - 0.5 * month) + c(0, year)), 0)
          pevtimesteps <- warmup + program_start # starting when warmup ends
          hybridboost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) 
          
          params <- set_pev_epi(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps,
            coverages = RTSScov,
            age = round(5 * month),
            min_wait = 0,
            booster_timestep = hybridboosters,
            booster_coverage = hybridboost_cov,
            booster_profile = list(rtss_booster_profile),
            seasonal_boosters = TRUE) 
        }
      
      
        # Get timing for mass vaccination rounds
        if(RTSS == 'SVmass+EPI' | RTSS == 'SVmass+hybrid'){
          first <- round(warmup + program_start + (peak - month * 3.5), 0) # start a year later
        } else if(RTSS == 'mass+EPI'){
          first <- round(warmup + program_start)
        }
        
        if (RTSSrounds == 'single'){
          pevtimesteps <- c(first)
        } else if (RTSSrounds == 'every 3 years'){
          pevtimesteps <- c(first, first + seq(3 * year, sim_length, 3 * year))
        }
      
        # Next, set the mass vaccination strategy in addition to routine vaccination
        # if(RTSSage %in% c('young children', 'all children', 'under 5s')){
        #   min_wait = 3 * year  # I have no idea what the min wait would actually be - if vaccinated around 6 m, then next would be around 3.5-4 yrs old.
        # } else {
        #   min_wait = 1 * month
        # }
        
        if(RTSSage != 'everyone'){
          params <- set_mass_pev(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps,
            coverages = rep(RTSScov, length(pevtimesteps)),
            min_ages = min_ages,
            max_ages = max_ages,
            min_wait = 0,
            booster_timestep = massboosters,
            booster_coverage = boost_cov,#rep(boost_cov, length(massboosters)),
            booster_profile = rep(list(rtss_booster_profile), length(massboosters)))
        } else if(RTSSage == 'everyone'){
          proppop_notpregnant <- 1 - 0.078/2 # from DHS data - see Get_pregnancy_rate.R
          
          # add mass vaccination for taking into account pregnant women
          params <- set_mass_pev(
            parameters = params,
            profile = rtss_profile,
            timesteps = pevtimesteps, 
            coverages = rep(RTSScov * proppop_notpregnant, length(pevtimesteps)), 
            min_ages = min_ages,
            max_ages = max_ages,
            min_wait = 0,
            booster_timestep = massboosters, # timesteps following initial vaccination 
            booster_profile = rep(list(rtss_booster_profile), length(massboosters)),
            booster_coverage = boost_cov#rep(boost_cov, length(massboosters)) # prop of vaccinated pop who will receive booster vaccine
          )
        }
        
        # var for outputting to check RTS,S timings are correct
        print(paste0("Mass timesteps: ", mass_pev_timesteps <- params$mass_pev_timesteps - warmup))
        print(paste0('EPI timesteps: ', pev_epi_timesteps <- params$pev_epi_timesteps - warmup))
        print(paste0('EPI booster: ', params$pev_epi_booster_timestep))
        print(paste0('Mass booster: ', params$mass_pev_booster_timestep))
      }
      
      # hybrid ----------
      if (RTSS == "hybrid") {
        params$pev_doses <- round(c(0, 1 * month, 2 * month)) # spacing from phase iii trial R21
        
        peak <- peak_season_offset(params)
        first <- round(warmup + program_start + (peak - month * 3.5), 0)
        boosters <- if(fifth == 0) round(c(first - warmup + 3 * month), 0) else round(((first - warmup + 3 * month) + c(0, year)), 0)
        boost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) # coverage from 10.1016/S2214-109X(22)00416-8
        pevtimesteps <- warmup + program_start # starting when warmup ends
        
        params <- set_pev_epi(
          parameters = params,
          profile = rtss_profile,
          timesteps = pevtimesteps,
          coverages = RTSScov,
          age = round(5 * month),
          min_wait = 0,
          booster_timestep = boosters,
          booster_coverage = boost_cov,
          booster_profile = list(rtss_booster_profile),
          seasonal_boosters = TRUE
        ) 
      }
      
      if (RTSS == 'SV'){ # only to children 5-17 months of age
        params$pev_doses <- round(c(0, 1 * month, 2 * month)) # spacing from phase iii trial R21
        
        peak <- peak_season_offset(params)
        first <- round(warmup + program_start + (peak - month * 3.5), 0)
        pevtimesteps <- c(first)
        boosters <- if(fifth == 0) round(c(first - warmup + 3 * month), 0) else round(((first - warmup + 3 * month) + c(0, year)), 0)
        boost_cov <- if(fifth == 0) RTSScov * 0.8 else c(RTSScov*0.8, RTSScov*0.8*0.9) # coverage from 10.1016/S2214-109X(22)00416-8
        
        params <- set_mass_pev(
          params,
          profile = rtss_profile,
          timesteps = pevtimesteps, 
          coverages = rep(RTSScov, length(pevtimesteps)), 
          min_ages = min_ages,
          max_ages = max_ages,
          min_wait = 0,
          booster_timestep = boosters, # timesteps following initial vaccination 
          booster_profile = list(rtss_booster_profile),
          booster_coverage = rep(boost_cov, length(boosters)) # prop of vaccinated pop who will receive booster vaccine
        )
      }
      
      if (RTSS == 'catch-up'){
        
      }
    }
    
    # MDA  ----------
    if (MDAcov > 0 & MDAtiming != 'none') { # this will only run if there is MDA, which will not be the case if there is no vaccination 
      if (MDAtiming == 'during'){
        mdatimesteps <- pevtimesteps # round of MDA occurs at same time as first dose
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
    if (SMC > 0 & RTSS %in% c("EPI", "mass", "hybrid")) {
      
      params$rtss_beta <- 70.9
      params$rtss_alpha <- 0.868
      params$rtss_vmax <- 0.843
      params$rtss_cs_boost <- c(6.37008, 0.35)
      
      params$drug_prophylaxis_scale <- c(10.6, 45.76)
      params$drug_prophylaxis_shape <- c(11.3, 2.87)
      
    }
    
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
    #   #with malariasimulation@dev 
    #   # params <- params |> set_parameter_draw(sample(1:1000,1))
    #        
    #   d <- readRDS(paste0(path,'03_output/parameter_draws.rds'))[drawID,]
    #   
    #   # over-write malariasimulation parameters to match the parameter draw
    #   params$dd = d$dur_D
    #   params$dt = d$dur_T
    #   params$da = d$dur_A # value 195 from the old model is the right one! New says 200
    #   params$du = d$dur_U
    #   params$sigma_squared = d$sigma2
    #   params$rm = d$dm
    #   params$rvm = d$dvm
    #   params$rb = d$db
    #   params$rc = d$dc
    #   params$rva = d$dv
    #   params$rid =	d$dd
    #   params$b0 = d$bh
    #   params$b1 = d$bmin
    #   params$ib0 = d$IB0
    #   params$kb = d$kb
    #   params$ub = d$ub
    #   params$uc = d$uc
    #   params$uv = d$uv
    #   params$ud = d$ud
    #   params$cd = d$cD
    #   params$gamma1 = d$gamma_inf
    #   params$cu = d$cU
    #   params$ct = d$cT # values in malariasim are more precise
    #   params$a0 = d$a0
    #   params$rho = d$rho
    #   params$phi0 = d$phi0
    #   params$phi1 = d$phi1
    #   params$ic0 = d$IC0
    #   params$kc = d$kc
    #   params$theta0 = d$theta0
    #   params$theta1 = d$theta1
    #   params$kv = d$kv
    #   params$fv0 = d$fv0
    #   params$av = d$av0
    #   params$gammav = d$gammav
    #   params$iv0 = d$IV0
    #   params$de = d$dur_E
    #   params$delay_gam = d$latgam
    #   params$dem = d$latmosq
    #   params$fd0 = d$fd0
    #   params$ad = d$ad0
    #   params$gammad = d$gammad
    #   params$d1 = d$dmin / 1000
    #   params$id0 = d$ID0
    #   params$kd = d$kd
    #   params$average_age = round(1 / d$eta)
    #   params$pcm = d$P_IC_M
    #   params$pvm = d$P_IV_M
    #   
    # }
    
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
