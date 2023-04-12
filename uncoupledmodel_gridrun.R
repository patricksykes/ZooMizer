source("fZooMSS_CalculatePhytoParam.R")
enviro <- readRDS("data/enviro_grid_sst15.RDS")
enviro$dt <- 0.01
enviro$tmax <- 1000

source("uncoupledmodel.R")

require(mizer)
require(tidyverse)
require(foreach)
require(assertthat)

Groups <- read_csv("data/TestGroups_mizer.csv")

library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, source("uncoupledmodel.R")) %>% invisible()

zoomssgrid <- list()

zoomssgrid <- foreach(i=1:nrow(enviro),
                      .packages = c("assertthat", "mizer")
                      .export = c("Groups", "enviro")
                      ) %dopar% {

                        input <- enviro[i,]

                        phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
                          npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
                          npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
                          return(npp)
                        }
                        
                        fZooMizer_run(groups = Groups, input = input, no_w = 177+1)
                      }

saveRDS(zoomssgrid, file = "20230314_sst15_ZooMSS.RDS")