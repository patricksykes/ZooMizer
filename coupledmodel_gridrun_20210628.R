# Script to run ZooMizer coupled model over a range of environmental conditions
# Patrick Sykes
# 27 June 2021

require(mizer)
require(tidyverse)
require(foreach)
require(assertthat)

ZooMizer_coupled <- function(ID, tmax = 25, effort = 0) {
  groups <- read.csv("data/TestGroups_mizer.csv") #load in ZooMSS groups
  zoo_groups <- groups[which(groups$Type == "Zooplankton"),] #just keep the zooplankton
  
  input <- read_csv("data/enviro_grid20210628.csv")[ID,]
  
  stable_zoomizer <- readRDS("test_grid_20210627.RDS")[[ID]]
  
  times <- length(getTimes(stable_zoomizer))
  stable_zoo <- colSums(stable_zoomizer@n[ceiling(times/2+1):times,,])
  
  stable_zoomizer@n <- sweep(stable_zoomizer@n, 2, 2.5 * stable_zoomizer@params@species_params$Carbon, "*") #convert to carbon
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
  
  fish_params <- newTraitParams(no_sp = 5,
                                min_w = 10^(-3),
                                min_w_inf = 10^3,
                                max_w_inf = 10^7,
                                no_w = (7-(-3))*10+1, # sets dx = 0.1 basically
                                min_w_pp = 10^(-14.5),
                                alpha = 0.25, # * temp_eff, #takes care of assimilation efficiency & temp effect on Encounter
                                kappa = exp(intercept), #10^(input$phyto_int) * 10^-1,
                                lambda = 1 - slope, # 1 - input$phyto_slope,
                                gamma = 1280 * temp_eff, #takes care of temp effect on PredRate and Encounter
                                # f0 = 0.6,
                                # h = 10^50,
                                R_factor = 1.01, #RMax = RFactor * RDI, default 4. Note RDI 
                                ks = 0, #set metabolism to zero
                                ext_mort_prop = 0, #currently zeroed since this is fitted. No need for temp effect here, calculates from PredMort
                                knife_edge_size = 10  # min size for fishing
  )
  
  # fish_params@species_params$gamma <- fish_params@species_params$gamma * temp_eff
  # fish_params <- setSearchVolume(fish_params)
  
  zoo_params <- newZooMizerParams(groups = zoo_groups, input = input, fish_params = fish_params)
  
  zoo_params@initial_n <- new[which(groups$Type == "Zooplankton"),]
  
  fish_params <- fish_params %>%
    setComponent(component = "zoo",
                 initial_value = zoo_params@initial_n,
                 dynamics_fun = "zoo_dynamics",
                 component_params = list(params = zoo_params)
    ) %>%
    setResource(resource_dynamics = "resource_zooMizer")
  
  # this line not needed if gamma is set to be temperature-dependent
  # fish_params <- setRateFunction(fish_params, "PredRate", "PredRate_temp")
  
  fish_params@other_params$temp_eff <- 2.^((input$sst - 30)/10)
  
  initialNResource(fish_params) <- resource_zooMizer(fish_params, fish_params@initial_n_other) #TODO: set this to be stable state of a regular ZooMizer run
  
  plotSpectra(fish_params, wlim = c(10^-14.5, NA), total = FALSE)
  
  initialN(fish_params) <- get_initial_n(fish_params)
  fish_params <- setParams(fish_params)
  
  return(project(fish_params, t_max = tmax, dt = 0.1, effort = effort))
}  



library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, lapply(c("mizer", "assertthat", "tidyverse"), require, character.only = TRUE)) %>% invisible()
clusterExport(cl, c("ZooMizer_coupled", "PredRate_temp"))



sims <- foreach(ID = 1:24,
                              .packages = c("mizer", "assertthat", "tidyverse")
                              #.export = c("ZooMizer_coupled", "PredRate_temp")
) %dopar% {
  source("ZooMizerResourceFunctions.R")
  ZooMizer_coupled(ID, tmax = 150, effort = 0)
}

saveRDS(sims, "coupledmodel_24grid_20210628.RDS")

biomasses <- foreach(i=1:length(sims), .combine = rbind) %dopar% sum(colMeans(tail(getBiomass(sims[[i]]),25)))


