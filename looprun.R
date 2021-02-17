require(mizer)
require(assertthat)
source("fZooMizer_run.R")
#source("ZooMizer_setup.R")

environment(new_project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", new_project_simple, ns = "mizer")


Groups <- read.csv("data/TestGroups_mizer.csv")

enviro <- readRDS("data/enviro_test20.RDS")
enviro$dt <- 0.01
enviro$tmaxx <- 1000

assim_eff <- matrix(Groups$GrossGEscale * groups$Carbon, nrow = nrow(Groups), ncol = nrow(Groups))
phyto_cc <- 0.1
assim_phyto <- Groups$GrossGEscale * phyto_cc
get_filterfeeders <- which(Groups$FeedType == "FilterFeeder")

for (i in get_filterfeeders) {
  assim_eff[,i] <- assim_eff[,i] / Groups$Carbon[i]
  assim_phyto[i] <- assim_phyto[i] / Groups$Carbon[i]
}

zoomizergrid <- list()


for (i in 1:nrow(enviro)) {
  phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
    npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
    npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    return(npp)
  }
  
  sim <- fZooMizer_run(groups = Groups, input = enviro[i,])
  zoomizergrid[[i]] <- sim
  rm(sim)
}

saveRDS(zoomizergrid, file="test_grid.RDS", version = 2)

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
# zgrid <- list()
# j <- 0
# for (i in enviro$cellID) {
#   j <- j+1
#   zgrid[[j]] <- zoomssgrid[[i]]
# }
# 
# saveRDS(zgrid, "zoomssgrid.rds")
