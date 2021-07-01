# Script to generate tile plots of the ZooMizer output
# Patrick Sykes
# 29 June 2021

library(mizer)
library(assertthat)
library(tidyverse)


today <- format(Sys.Date(), format = "%Y%m%d")

sims <- read_rds("coupledmodel_24grid_20210628.RDS")
enviro <- read_csv("data/enviro_grid20210628.csv")

zoo_groups <- read_csv(file = "data/TestGroups_mizer.csv")[1:9,]
zoo_params <- newZooMizerParams(groups = zoo_groups, input = enviro[1,], fish_params = sims[[1]]@params)
  
all_bioms <- foreach(i=1:24, .combine = rbind) %dopar% colMeans(tail(getBiomass_ZooMizer(sims[[i]], zoo_params),25))

biom.df <- cbind(enviro, all_bioms) %>% 
  rename(fish1=`1`,
         fish2=`2`,
         fish3=`3`,
         fish4=`4`,
         fish5=`5`) %>% 
  mutate(AllFish = fish1 + fish2 + fish3 + fish4 +fish5)

library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, lapply(c("mizer", "assertthat", "tidyverse"), require, character.only = TRUE)) %>% invisible()

tiles <- foreach(i = 9:23) %dopar% {
  column <- sym(colnames(biom.df)[i])
  ggplot(biom.df, aes(x=sst, y = log10(chlo)))+
    geom_tile(aes(fill=log10(!!column)))+
    scale_fill_viridis_c()+
    labs(title = paste(colnames(biom.df)[i], "Biomass vs. SST and chlorophyll-a"))
  ggsave(filename = paste0("tileplot_",column,"_",today,".png"))
}

slopes <- foreach(i=1:24, .combine = rbind) %dopar% colMeans(tail(getCommunitySlope(sims[[i]]),25))[1]
biom.df$fish_slope <- as.numeric(slopes)

ggplot(biom.df, aes(x = phyto_slope, y = fish_slope))+
  geom_point(aes(colour = sst))+
  geom_smooth()+
  geom_abline(slope = 1, intercept = 0)+
  labs(title = "Fish community slope vs. phyto slope")
ggsave(filename = paste0("phytovsfishslopes_",column,"_",today,".png"))


ggplot(biom.df, aes(x = sst, y = fish_slope/phyto_slope))+
  geom_point(aes(colour = log10(chlo)))+
  geom_smooth()

source("ZooMizerPlots.R")


# problem here where fish#5 doesn't have the line shown...
timeseries <- foreach(i = 1:24) %dopar% {
  plotBiomass_ZooMizer(sims[[i]], zoo_params)+
    labs(title = paste0("Species biomass over time, SST = ", enviro$sst[i], ", chlo = ", enviro$chlo[i]))
}

  