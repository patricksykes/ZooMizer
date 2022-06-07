# Optimise ZooMizer parameters

source("uncoupledmodel.R")
source("ZooMizerResourceFunctions.R")

ZooMizer_coupled_gg <- function(ID, tmax = 1000, effort = 0, GGEvect) {
  groups <- read.csv("data/TestGroups_mizer.csv") #load in ZooMSS groups
  
  # Run the ZooMSS model inside Mizer:
  #lower "starting proportion" of salps and chaetognaths
  # groups$Prop[c(7,8)] <- 0.1*groups$Prop[c(7,8)]
  # #lower "starting proportion" of euphausiids
  # groups$Prop[6] <- 0.1*groups$Prop[c(6)]
  
  # change start of senescence for salps and euphs
  # groups$w_mat[c(6,8)] <- groups$w_mat[c(6,8)] - 2  #note these are in log10(w)
  
  zoo_groups <- groups[which(groups$Type == "Zooplankton"),] #just keep the zooplankton
  
  input <- readRDS("data/enviro_grid20.RDS")[ID,]
  
  stable_zoomizer <- readRDS(paste0("Output/", zoomssjobname, "_ZooMSS_", ID_char,".RDS"))
  
  times <- length(getTimes(stable_zoomizer))
  stable_zoo <- colMeans(stable_zoomizer@n[ceiling(times/2+1):times,,]) 
  
  # stable_zoomizer@n <- sweep(stable_zoomizer@n, 2, 2.5 * stable_zoomizer@params@species_params$Carbon, "*") #convert to carbon
  intercept <- mean(getCommunitySlope(stable_zoomizer,
                                      species = which(groups$Type == "Zooplankton"),
                                      max_w = 1e-5,
  )$intercept[ceiling(times/2+1):times])
  slope <- mean(getCommunitySlope(stable_zoomizer,
                                  species = which(groups$Type == "Zooplankton"),
                                  max_w = 1e-5,
  )$slope[ceiling(times/2+1):times])
  
  
  new <- stable_zoo[which(groups$Type == "Zooplankton"),]
  
  temp_eff <- 2.^((input$sst - 30)/10)
  R_factor = 4
  fish_params <- newTraitParams(no_sp = 13,
                                min_w = 10^(-3),
                                min_w_inf = 10^1,
                                max_w_inf = 10^7,
                                no_w = 10*(7-(-3))+1, # sets dx = 0.1
                                min_w_pp = 10^(-14.5),
                                alpha = 0.25, # * temp_eff, #takes care of assimilation efficiency & temp effect on Encounter
                                kappa = exp(intercept), #10^(input$phyto_int) * 10^-1,
                                lambda = 1 - slope, # 1 - input$phyto_slope,
                                gamma = 640 * temp_eff, #takes care of temp effect on PredRate and Encounter
                                # f0 = 0.6,
                                h = 40,
                                beta = 100,
                                sigma = 1.3,
                                # R_factor = 1.01, #RMax = RFactor * RDI, default 4. Note RDI 
                                reproduction_level = 1/R_factor,
                                fc = 0,
                                ks = 0, #set metabolism to zero
                                ext_mort_prop = 0, #currently zeroed since this is fitted. No need for temp effect here, calculates from PredMort
                                knife_edge_size = 10  # min size for fishing
  )
  
  # fish_params@species_params$gamma <- fish_params@species_params$gamma * temp_eff
  # fish_params <- setSearchVolume(fish_params)


  zoo_params <- newZooMizerParams(groups = zoo_groups, input = input, fish_params = fish_params)
  
  zoo_params@initial_n <- new[which(groups$Type == "Zooplankton"),]
  
  zoo_params@species_params$interaction_resource[9] <- 1 # make jellies eat resource (i.e. fish)
  zoo_params@species_params$GrossGEscale <- GGEvect
  ## Zooplankton Mortality
  fish_params@other_params$temp_eff <- 2.^((input$sst - 30)/10)
  
  
  #U-shaped mortality
  M <- getExtMort(zoo_params) # pulls out senescence
  z0pre <- rep(0, nrow(M))
  z0pre[c(6,8)] <- 1 * (zoo_params@species_params$w_inf[c(6,8)])^(1/4) * fish_params@other_params$temp_eff #choose groups to apply to
  allo_mort <- outer(z0pre, w(zoo_params)^(-1/4)) # mortality of z0pre * w^(-1/4)
  zoo_params <- setExtMort(zoo_params, ext_mort = M + allo_mort)
  
  # hockey stick mortality (mizer background mortality + senescence)
  # M <- getExtMort(zoo_params) # pulls out senescence
  # z0pre <- rep(0, nrow(M))
  # z0pre[c(6,8)] <- 0.6 #choose groups to apply to; default Mizer value 0.6
  # mu_b <- z0pre * zoo_params@species_params$w_inf^(-1/4) * fish_params@other_params$temp_eff
  # zoo_params <- setExtMort(zoo_params, ext_mort = M + mu_b) # adds 0.6 * w_inf^(-1/4)
  
  #Both
  # M <- getExtMort(zoo_params) # pulls out senescence
  # z0pre <- rep(0, nrow(M))
  # z0pre[c(6,8)] <- 1 * (zoo_params@species_params$w_inf[c(6,8)])^(1/4) * fish_params@other_params$temp_eff #choose groups to apply to
  # larv_mort <- outer(z0pre, w(zoo_params)^(-1/4)) # mortality of z0pre * w^(-1/4)
  # 
  # z0pre1 <- rep(0, nrow(M))
  # z0pre1[c(6,8)] <- 6 * fish_params@other_params$temp_eff #choose groups to apply to; default Mizer value 0.6
  # back_mort <- z0pre * zoo_params@species_params$w_inf^(-1/4)
  # zoo_params <- setExtMort(zoo_params, ext_mort = M + larv_mort + back_mort) # adds 0.6 * w_inf^(-1/4)!
  
  # #type II feeding for zoo
  # zoo_params <- setRateFunction(zoo_params, "PredRate", "mizerPredRate")
  # zoo_params <- setRateFunction(zoo_params, "FeedingLevel", "mizerFeedingLevel")
  # zoo_params@species_params$h <- 40
  # zoo_params <- setMaxIntakeRate(zoo_params)
  
  zoo_params@species_params$mu0DD <- rep(80, 9) * zoo_params@species_params$w_inf
  zoo_params@species_params$mu0DD[6] <- zoo_params@species_params$mu0DD[6] * 50
  zoo_params <- setRateFunction(zoo_params, "Mort", "totalMortDD")
  
  fish_params <- fish_params %>%
    setComponent(component = "zoo",
                 initial_value = zoo_params@initial_n,
                 dynamics_fun = "zoo_dynamics",
                 mort_fun = "zoo_predation",
                 component_params = list(params = zoo_params)
    ) %>%
    setResource(resource_dynamics = "resource_zooMizer")
  
  
  
  # this line not needed if gamma is set to be temperature-dependent
  # fish_params <- setRateFunction(fish_params, "PredRate", "PredRate_temp")
  
  
  
  initialNResource(fish_params) <- resource_zooMizer(fish_params, fish_params@initial_n_other) #TODO: set this to be stable state of a regular ZooMizer run
  # plotSpectra(fish_params, wlim = c(10^-14.5, NA), total = FALSE)
  
  # initialN(fish_params) <- get_initial_n(fish_params) * 10^3 # adjust initial n for fish
  initialNResource(fish_params@other_params$zoo$params) <- initialNResource(fish_params@other_params$zoo$params)
  fish_idx <- (length(fish_params@w_full)-length(fish_params@w) + 1):length(fish_params@w_full)
  initialNResource(fish_params@other_params$zoo$params)[fish_idx] <- initialNResource(fish_params@other_params$zoo$params)[fish_idx] + colSums(initialN(fish_params))
  # fish_params <- setParams(fish_params)
  
  # fish_params@species_params$R_max <- readRDS("data/rmaxs.RDS")/10
  fish_params@species_params$erepro <- 1
  # fish_params@species_params$R_max <- fish_params@resource_params$kappa * fish_params@species_params$w_inf ^(-1.5)
  # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10^(5+3*log10(input$chlo))
  fish_params <- setBevertonHolt(fish_params, erepro = 1)
  # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10000
  # fish_params <- setRateFunction(fish_params, "FeedingLevel", "FeedingLevel_type3")
  # fish_params <- setExtMort(fish_params, ext_mort = array(data=0.1, dim=dim(fish_params@mu_b)))
  
  return(project(fish_params, t_max = tmax, dt = 0.01, effort = effort))
}  

optimGGEfun <- function(ID, GGEvect) {
  # Groups <- read.csv("data/TestGroups_mizer.csv") # Load in functional group information
  # 
  # Groups$GrossGEscale <- c(GGEvect,2.5,2.5,2.5)
  # 
  #   # ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX'))
  # 
  
  # Choose environmental data to use
  enviro <- readRDS("data/enviro_grid20.RDS")[ID,]
  enviro$tmax <- 10
  
  
  out <- ZooMizer_coupled_gg(ID, tmax = 10, effort = 0, GGEvect = GGEvect)
  
  params <- setInitialValues(out@params, out)
  times <- max(getTimes(out))
  params@initial_n[] <- colMeans(out@n[ceiling(times/2+1):times,,])
  nothers <- array(0,c(times,dim(params@initial_n_other[[1]])))
  for (i in 1:times) {
    nothers[i,,]  <- out@n_other[[i]]
  }
  params@initial_n_other$zoo[] <- colMeans(nothers[ceiling(times/2+1):times,,])
  params@other_params$zoo$params@initial_n_other[] <- colMeans(nothers[ceiling(times/2+1):times,,])
  params@initial_n_pp <- colMeans(out@n_pp[ceiling(times/2+1):times,])
  
  growth <- getEGrowth(params@other_params$zoo$params)
  for(i in 1:9) growth[i,1:params@other_params$zoo$params@w_min_idx[i]] <- 0
  w <- params@other_params$zoo$params@w
  
  
  gg_all <- growth*(2^((15-30)/10)/2^((as.numeric(enviro['sst']) - 30)/10))
  
  gg <- gg_all/365 # Growth rates in years from model, convert to days
  gg[c(3,8),] <- 15*gg[c(3,8),] # Multiply larvs and salps by 15 ()
  gg <- sweep(gg, 2, w, "/")
  return(gg)
}

compflag <- function(gg){
  flag_gs <- read.table('data/FlagHirst.txt', header = TRUE)
  flag_gs$BodyWeight <- (flag_gs$BodyWeight/1e6)*(1/0.15)
  flag_gs$Growth <- flag_gs$Growth*24

  gg <- melt(gg) %>% rename(BodyWeight = w, GROUP = sp, Growth = value) %>% filter(GROUP == "Flagellates" & Growth > 0)
  
  ggfun <- stats::approxfun(gg$BodyWeight, gg$Growth)
  
  flag_gs$pred <- ggfun(flag_gs$BodyWeight)
  err <- (flag_gs$pred-flag_gs$Growth)^2
  err <- sum(err, na.rm = TRUE)
  return(err)
}

optimgg_flag <- function(GGEflag){
  GGEvect <- c(GGEflag, rep(2.5,8))
  gg <- optimGGEfun(ID=10, GGEvect=GGEvect)
  return(compflag(gg))
}

zoomssjobname <- '20220118_alldt_pt01' #job name used on queue
ID <- 10
ID_char <- sprintf("%04d",ID)

# k <- optimgg_flag(ID)
# gg <- optimGGEfun(ID)

errs <- foreach(g=seq(from=2, to=3, by=0.1), .combine = rbind) %do% {
  phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^-1.584554, lambda=1.972042, ... ) {
    npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
    npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    return(npp)
  }
  
  err <- optimgg_flag(g)
  c(g, err)
}

flagval <- optim(par=2.5,optimgg_flag,method ="L-BFGS-B",lower=2,upper =3)



ggs <- optimGGEfun(10, c(flagval$par, rep(2.5,8)))
