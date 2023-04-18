# Script to run ZooMizer coupled model over a range of environmental conditions
# Patrick Sykes
# 1 July 2021
# Updated 10 March 2023

require(mizer)
require(tidyverse)
require(foreach)
require(assertthat)



library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, lapply(c("mizer", "assertthat", "tidyverse"), require, character.only = TRUE)) %>% invisible()
clusterExport(cl, c("ZooMizer_coupled"))

jobname <- "20230418_erepro_1_Rmax_Julia_sst15"

sims <- foreach(ID = 1:21,
                .packages = c("mizer", "assertthat", "tidyverse")
                # .export = c("ZooMizer_coupled")
) %dopar% {
  start <- Sys.time()
  zoomssjobname <- "20230404_sst15"
  ID_char <- sprintf("%04d",ID)
  source("uncoupledmodel.R")
  source("ZooMizerResourceFunctions.R")

  # environment(new_project_simple) <- asNamespace('mizer')
  # assignInNamespace("project_simple", new_project_simple, ns = "mizer")

  # Run the ZooMSS model inside Mizer:
  # phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
  #   npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  #   npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  #   return(npp)
  # }

  # density-dependent mortality
  # totalMortDD <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...){
  #   mort <- pred_mort + params@mu_b + f_mort
  #   # Add contributions from other components
  #   # for (i in seq_along(params@other_mort)) {
  #   #   mort <- mort +
  #   #     do.call(params@other_mort[[i]],
  #   #             list(params = params,
  #   #                  n = n, n_pp = n_pp, n_other = n_other, t = t,
  #   #                  component = names(params@other_mort)[[i]], ...))
  #   # }
  #   ddmort_raw <- sweep(n, "w", params@w^params@species_params$q[1] * params@dw, "*")
  #   ddmort <- ddmort_raw * 0
  #   sp_sel <- 1:nrow(ddmort) #c(6,7,8,9)
  #   mu0 <- 80 # Benoit & Rochet value = 80, but no density-independent mortality in that model
  #   ddmort[sp_sel,] <- mu0 * ddmort_raw[sp_sel,] #* params@species_params$w_inf
  #   ddmort[6,] <- ddmort[6,] * 50
  #   return(mort + ddmort)
  # }

  # zoomss <- fZooMizer_run(groups = read.csv("data/TestGroups_mizer.csv"), input =  readRDS("data/enviro_sst_15_chlox25.RDS")[ID,], no_w = (1*177+1)) # dx=0.1
  # saveRDS(zoomss, file = paste0("Output/", zoomssjobname, "_ZooMSS_", ID_char,".RDS"))

  mid <- Sys.time()
  
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
    
    input <- readRDS("data/enviro_sst_15_chlox25.RDS")[ID,]
    
    stable_zoomizer <- readRDS(paste0("Output/", zoomssjobname, "_ZooMSS_", ID_char,".RDS"))
    
    
    times <- length(getTimes(stable_zoomizer))
    stable_zoomizer@n <- sweep(stable_zoomizer@n, 2, 2.5 * stable_zoomizer@params@species_params$Carbon, "*") #convert to carbon
    intercept <- mean(getCommunitySlope(stable_zoomizer,
                                        species = which(groups$Type == "Zooplankton"),
                                        max_w = 1e-5
    )$intercept[ceiling(times/2+1):times])
    slope <- mean(getCommunitySlope(stable_zoomizer,
                                    species = which(groups$Type == "Zooplankton"),
                                    max_w = 1e-5
    )$slope[ceiling(times/2+1):times])
    new <- colMeans(stable_zoomizer@n[ceiling(times/2+1):times,,])[which(groups$Type == "Zooplankton"),]
    
    
    ## Alternative method using existing ZooMizer run
    # stable_zoomizer <- readRDS("zooplankton_sst15.RDS")[[ID]]
    
    # stable_zoo <- colMeans(stable_zoomizer@params@initial_n_other$zoo) #@n[ceiling(times/2+1):times,,]) 
    # sel <- which(stable_zoomizer@params@other_params$zoo$params@w <= 1e-5)
    # fit <- lm(log(stable_zoomizer@params@other_params$zoo$params@w[sel]) ~ log(stable_zoo[sel]))
    # intercept <- fit$coefficients[1]
    # slope <- fit$coefficients[2]
    
    # new <- stable_zoomizer@params@initial_n_other$zoo[which(groups$Type == "Zooplankton"),]
    
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
    initialNResource(fish_params@other_params$zoo$params) <- exp(intercept)*fish_params@w_full^(slope) / fish_params@dw_full #returns the fixed spectrum at every time step
    initialNResource(fish_params@other_params$zoo$params)[fish_params@w_full > fish_params@other_params$zoo$params@resource_params$w_pp_cutoff* (1 - 1e-06)] <- 0
    fish_idx <- (length(fish_params@w_full)-length(fish_params@w) + 1):length(fish_params@w_full)
    initialNResource(fish_params@other_params$zoo$params)[fish_idx] <- initialNResource(fish_params@other_params$zoo$params)[fish_idx] + colSums(initialN(fish_params))
    # fish_params <- setParams(fish_params)
    
    initialNResource(fish_params) <- resource_zooMizer(params = fish_params, n_other = fish_params@initial_n_other) #TODO: set this to be stable state of a regular ZooMizer run
    # plotSpectra(fish_params, wlim = c(10^-14.5, NA), total = FALSE)
    
    
    # fish_params@species_params$R_max <- readRDS("data/rmaxs.RDS")/10
    fish_params@species_params$erepro <- 1
    fish_params@species_params$R_max <- fish_params@resource_params$kappa * fish_params@species_params$w_inf^(-1.5) #relationship from Julia's model
    # fish_params@species_params$erepro <- readRDS("data/hpcmeanerepros.RDS") # mean of big HPC run - benefit of same across all runs
    # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10^(5+3*log10(input$chlo)) #scale R_max with chlorophyll
    # fish_params <- setBevertonHolt(fish_params, erepro = 1)
    # fish_params@species_params$R_max <- fish_params@species_params$R_max * 10000get
    # fish_params <- setRateFunction(fish_params, "FeedingLevel", "FeedingLevel_type3")
    # fish_params <- setExtMort(fish_params, ext_mort = array(data=0.1, dim=dim(fish_params@mu_b)))
    
    return(project(fish_params, t_max = tmax, dt = 0.01, effort = effort))
  }  
  
  environment(project_simple) <- asNamespace('mizer')
  assignInNamespace("project_simple", project_simple, ns = "mizer")
  
  sim <- ZooMizer_coupled(ID, tmax = 1000, effort = 0, zoomssjobname)
  end <- Sys.time()
  sim@params@other_params$time <- data.frame(zoomss = mid-start, zoomizer = end-mid)
  sim
}

saveRDS(sims, paste0("Output/",jobname,"_ZooMizer.RDS"))

enviro <- readRDS("data/enviro_sst_15_chlox25.RDS")
zoo_params <- readRDS("Output/20220118_alldt_pt01_ZooMSS_0001.RDS")@params %>% removeSpecies(c("Fish_Small", "Fish_Med", "Fish_Large"))

# sims <- foreach(ID = 1:20,
#                 .packages = c("mizer", "assertthat", "tidyverse")
#                 #.export = c("ZooMizer_coupled", "PredRate_temp")
# ) %dopar% {
#   ID_char <- sprintf("%04d",ID)
#   readRDS(paste0("Output/20220212_type1approx_ZooMizer_", ID_char,".RDS"))
# }

# check biomass time series plots
timeseries <- foreach(ID = 1:21) %dopar% {
  source("ZooMizerPlots.R")
  source("ZooMizerSummaryFunctions.R")
  plotBiomass_ZooMizer(sims[[ID]], zoo_params)+
    labs(title = paste0("SST = ", enviro$sst[ID], ", chlo = ", round(enviro$chlo[ID],3)))
}

library(patchwork)
timeseriesgrid <- wrap_plots(timeseries, nrow =  7, ncol = 3) + plot_layout(guides = "collect")
timeseriesgrid
ggsave(paste0(jobname,"_timeseriesplots.png"), timeseriesgrid, width = 15, height = 28)
# biomasses <- foreach(i=1:length(sims), .combine = rbind) %dopar% sum(colMeans(tail(getBiomass(sims[[i]]),25)))

#averaged size spectra plots
spectra <- foreach(i = 1:21) %dopar% {
  plotSpectra_ZooMizer(sims[[i]], zoo_params, time_range = c(501,1000), wlim = c(1e-14, NA))+
    labs(title = paste0("SST = ", enviro$sst[i], ", chlo = ", round(enviro$chlo[i],3)))
}

spectplots <- wrap_plots(spectra, nrow =  7, ncol = 3) + plot_layout(guides = "collect")
ggsave(paste0(jobname,"spectraplots.png"), spectplots, width = 15, height = 28)


#get biomass data
b <- foreach(ID=1:21, .combine = rbind) %dopar% {
  bioms <- getBiomass_ZooMizer(sims[[ID]], sims[[ID]]@params@other_params$zoo$params)
  rows <- ceiling(nrow(bioms)/2):nrow(bioms)
  colMeans(bioms[rows,])
}
saveRDS(b, paste0(jobname,"_biomdf.RDS"))

df <- cbind(enviro[1:21,], b)
propdf <- df
propdf[,c(3:5,11:12)] <- propdf[,c(3:5,11:12)] / rowSums(propdf[,c(3:5,11:12)])
propdf[,13:19] <- propdf[,13:19] / rowSums(propdf[,13:19])
propdf[,20:32] <- propdf[,20:32] / rowSums(propdf[,20:32])

# plot group biomass vs chlo
# Microzoo/phyto
df2 <- gather(propdf, Group, Proportion, c(pico_biom:micro_biom,Flagellates:Ciliates))
gg1 <- ggplot(df2, aes(x=log10(chlo), y = Proportion, fill = Group)) + 
  geom_area() +
  labs(title = "Microzoo proportion")

#zoo
df3 <- gather(propdf, Group, Proportion, Larvaceans:Jellyfish)
gg2 <- ggplot(df3, aes(x=log10(chlo), y = Proportion, fill = Group)) + 
  geom_area() +
  labs(title = "Zooplankton proportion")

# fish
df4 <- gather(propdf, Group, Proportion, "1":"13")
gg3 <- ggplot(df4, aes(x=log10(chlo), y = Proportion, fill = Group)) + 
  geom_area() +
  labs(title = "Fish proportion")

#combine and save plots
biomchlo <- wrap_plots(gg1,gg2,gg3, nrow =  3, ncol = 1) + plot_layout(guides = "collect")
ggsave(paste0(jobname, "_biomvschlo.png"), biomchlo, width = 18, height = 15)
