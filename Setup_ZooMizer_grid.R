## Run the ZooMizer coupled model over multiple cores
##
## Adapted from Jason's ZooMSS code for UNSWs Katana
##
## Updated Tuesday 4th of January 2022

library(mizer)
library(tidyverse)
library(assertthat)

#job specifics
Groups <- read.csv("data/TestGroups_mizer.csv") # Load in functional group information

jobname <- '20220104_grid' #job name used on queue

ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX'))
ID_char <- sprintf("%04d",ID)

# Choose environmental data to use
enviro <- readRDS("data/enviro_grid20210705.RDS")[ID,]

source("uncoupledmodel.R")
source("ZooMizerResourceFunctions.R")

environment(new_project_simple) <- asNamespace('mizer') 
assignInNamespace("project_simple", new_project_simple, ns = "mizer")

# Run the ZooMSS model inside Mizer:
phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
  npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}

zoomss <- fZooMizer_run(groups = Groups, input = enviro, no_w = 177+1) # dx=0.1
saveRDS(zoomss, file = paste0("Output/", jobname, "_ZooMSS_", ID_char,".RDS"))



environment(project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", project_simple, ns = "mizer")

ZooMizer_coupled <- function(ID, tmax = 1000, effort = 0, zoomssjobname) {
  groups <- read.csv("data/TestGroups_mizer.csv") #load in ZooMSS groups
  ID_char <- sprintf("%04d",ID)
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
  fish_params <- newTraitParams(no_sp = 13,
                                min_w = 10^(-3),
                                min_w_inf = 10^1,
                                max_w_inf = 10^7,
                                no_w = 10*(7-(-3))+1, # sets dx = 0.1
                                min_w_pp = 10^(-14.5),
                                alpha = 0.6, # * temp_eff, #takes care of assimilation efficiency & temp effect on Encounter
                                kappa = exp(intercept), #10^(input$phyto_int) * 10^-1,
                                lambda = 1 - slope, # 1 - input$phyto_slope,
                                gamma = 640 * temp_eff, #takes care of temp effect on PredRate and Encounter
                                # f0 = 0.6,
                                h = 40,
                                beta = 100,
                                sigma = 1.3,
                                # R_factor = 1.01, #RMax = RFactor * RDI, default 4. Note RDI 
                                reproduction_level = 0,
                                fc = 0,
                                ks = 0, #set metabolism to zero
                                ext_mort_prop = 0, #currently zeroed since this is fitted. No need for temp effect here, calculates from PredMort
                                knife_edge_size = 10  # min size for fishing
  )
  
  # fish_params@species_params$gamma <- fish_params@species_params$gamma * temp_eff
  # fish_params <- setSearchVolume(fish_params)
  
  fish_params@other_params$temp_eff <- 2.^((input$sst - 30)/10)
  
  zoo_params <- newZooMizerParams(groups = zoo_groups, input = input, fish_params = fish_params)
  
  zoo_params@initial_n <- new[which(groups$Type == "Zooplankton"),]
  
  zoo_params@species_params$interaction_resource[9] <- 1 # make jellies eat resource (i.e. fish)
  
  ## Zooplankton Mortality
  
  #U-shaped mortality
  M <- getExtMort(zoo_params) # pulls out senescence
  z0pre <- rep(0, nrow(M))
  z0pre[c(6,8)] <- 1 * (zoo_params@species_params$w_inf[c(6,8)])^(1/4) * fish_params@other_params$temp_eff #choose groups to apply to
  allo_mort <- outer(z0pre, w(zoo_params)^(-1/4)) # mortality of z0pre * w^(-1/4)
  zoo_params <- setExtMort(zoo_params, ext_mort = M + allo_mort)
  
  # add fixed level of background mortality (so in total have mizer fixed background mortality + senescence)
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
  
  #Density-dependent mortality
  zoo_params@species_params$mu0DD <- rep(80, 9) * zoo_params@species_params$w_inf
  zoo_params@species_params$mu0DD[6] <- zoo_params@species_params$mu0DD[6] * 50 #turn up density-dependent mortality on krill
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
  
  
  
  # initialN(fish_params) <- get_initial_n(fish_params) * 10^3 # adjust initial n for fish
  initialNResource(fish_params@other_params$zoo$params) <- kappa*params@w_full^(1-lambda) / params@dw_full #returns the fixed spectrum at every time step
  initialNResource(fish_params@other_params$zoo$params)[fish_params@w_full > fish_params@other_params$zoo$params@resource_params$w_pp_cutoff* (1 - 1e-06)] <- 0
  fish_idx <- (length(fish_params@w_full)-length(fish_params@w) + 1):length(fish_params@w_full)
  initialNResource(fish_params@other_params$zoo$params)[fish_idx] <- initialNResource(fish_params@other_params$zoo$params)[fish_idx] + colSums(initialN(fish_params))
  # fish_params <- setParams(fish_params)
  
  initialNResource(fish_params) <- resource_zooMizer(params = fish_params, n_other = fish_params@initial_n_other) #TODO: set this to be stable state of a regular ZooMizer run
  # plotSpectra(fish_params, wlim = c(10^-14.5, NA), total = FALSE)
  
  
  # fish_params@species_params$R_max <- readRDS("data/rmaxs.RDS")/10
  fish_params@species_params$erepro <- 1
  # fish_params@species_params$R_max <- fish_params@resource_params$kappa * fish_params@species_params$w_inf ^(-1.5)
  # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10^(5+3*log10(input$chlo))
  # fish_params <- setBevertonHolt(fish_params, erepro = 1)
  # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10000
  # fish_params <- setRateFunction(fish_params, "FeedingLevel", "FeedingLevel_type3")
  # fish_params <- setExtMort(fish_params, ext_mort = array(data=0.1, dim=dim(fish_params@mu_b)))
  
  return(project(fish_params, t_max = tmax, dt = 0.01, effort = effort))
} 

out <- ZooMizer_coupled(ID, tmax = 1000, effort = 0)
saveRDS(out, file = paste0("Output/", jobname, "_ZooMizer_", ID_char,".RDS"))
