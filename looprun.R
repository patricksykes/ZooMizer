require(mizer)
source("fZooMizer_run.R")
#source("ZooMizer_setup.R")

environment(new_project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", new_project_simple, ns = "mizer")


Groups <- read.csv("data/TestGroups_mizer.csv")
enviro <- readRDS("data/envirofull_20200317.RDS")

enviro <- enviro[enviro$sst >= 1,]
enviro$dt <- 0.01
enviro$tmaxx <- 1000

assim_eff <- matrix(groups$GrossGEscale * groups$Carbon, nrow = nrow(groups), ncol = nrow(groups))
get_filterfeeders <- which(groups$FeedType == "FilterFeeder")

for (i in get_filterfeeders) {
  assim_eff[,i] <- assim_eff[,i] / groups$Carbon[i]
}

zoomizergrid <- list()

enviro_subset=enviro[seq(71,1428,71),]

for (i in 1:nrow(enviro_subset)) {
  phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
    n_pp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
    n_pp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
    return(n_pp)
  }
  
  sim <- fZooMizer_run(groups = Groups, input = enviro_subset[i,])
  zoomizergrid[i] <- sim
  rm(sim)
}

saveRDS(zoomizergrid, file="test_grid.RDS", version = 2)