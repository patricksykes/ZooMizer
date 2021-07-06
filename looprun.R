require(mizer)
require(assertthat)
require(tidyverse)
source("fZooMizer_run.R")
source("fZooMSS_CalculatePhytoParam.R")
#source("ZooMizer_setup.R")

environment(new_project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", new_project_simple, ns = "mizer")

# environment(new_newMultispeciesParams) <- asNamespace('mizer')
# assignInNamespace("newMultispeciesParams", new_newMultispeciesParams, ns = "mizer")
# 
# environment(new_emptyParams) <- asNamespace('mizer')
# assignInNamespace("emptyParams", new_emptyParams, ns = "mizer")


Groups <- read.csv("data/TestGroups_mizer.csv")

# enviro <- readRDS("data/enviro_test20.RDS")
enviro <- expand.grid(sst = seq(0,30,10), chlo = 10^(seq(-1.5,1,.5)))
enviro <- fZooMSS_CalculatePhytoParam(enviro)
enviro$dt <- 0.01
enviro$tmax <- 1000

assim_eff <- matrix(Groups$GrossGEscale * groups$Carbon, nrow = nrow(Groups), ncol = nrow(Groups))
phyto_cc <- 0.1
assim_phyto <- Groups$GrossGEscale * phyto_cc
# get_filterfeeders <- which(Groups$FeedType == "FilterFeeder")

# for (i in get_filterfeeders) {
#   assim_eff[,i] <- assim_eff[,i] / Groups$Carbon[i]
#   assim_phyto[i] <- assim_phyto[i] / Groups$Carbon[i]
# }

zoomizergrid <- list()

library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, lapply(c("mizer", "assertthat", "tidyverse"), require, character.only = TRUE)) %>% invisible()
clusterEvalQ(cl, lapply(c("fZooMizer_run.R","fZooMSS_CalculatePhytoParam.R"), source)) %>% invisible()
clusterExport(cl, c("Groups", "enviro"))

zoomizergrid <- foreach(i =1:nrow(enviro),
                        .packages = c("mizer", "assertthat", "tidyverse")
                        # .export = c("Groups", "enviro", "phyto_fixed")
                        ) %dopar% {
      
# load required functions on the cluster instance
  environment(new_project_simple) <- asNamespace('mizer') 
  assignInNamespace("project_simple", new_project_simple, ns = "mizer")
  
  phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
    npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
    npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    return(npp)
  }
  
  #run the function
  
  fZooMizer_run(groups = Groups, input = enviro[i,])
  
}

saveRDS(zoomizergrid, file="test_grid_20210628.RDS", version = 2)

# #apples to apples comparison:
# 
# zoomizeraves <- list()
# 
# for (i in 1:nrow(enviro)) {
#   zoomizeraves[[i]] <- apply(as.array(zoomizergrid[[i]]@n[501:1001,,]),c(2,3),'mean')
# }
# 
# plotSpectra(hack2, time_range = 1)
# 

# load("Output/full_ZooMizer.RData")
# 
#  zgrid <- list()
#  j <- 0
#  for (i in 1:20) {
#    j <- j+1
#    zgrid[[j]] <- res[[i]]
#  }
#  
# saveRDS(zgrid, "zoomssgrid_redo_yearly.rds")
