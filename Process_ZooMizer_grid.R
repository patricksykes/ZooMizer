# Script to collect and process ZooMizer runs from HPC
# Patrick Sykes
# 8th of July 2021

library(mizer)
library(tidyverse)
library(assertthat)

#job specifics
Groups <- read.csv("data/TestGroups_mizer.csv") # Load in functional group information
zoo_groups <- Groups[1:9,]
zoo_params <- newZooMizerParams(groups = zoo_groups, input = enviro[1,], fish_params = sims[[1]]@params)



jobname <- '20210705_grid' #job name used on queue

enviro <- readRDS("data/enviro_grid20210705.RDS")


library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, lapply(c("mizer", "assertthat", "tidyverse"), require, character.only = TRUE)) %>% invisible()


biomasses <- foreach(ID=1:nrow(enviro), .combine = rbind) %dopar% {
  source("ZooMizerSummaryFunctions.R")
  ID_char <- sprintf("%04d",ID)
  sim <- readRDS(paste0("Output/",jobname,"_ZooMizer_",ID_char,".RDS"))
  colMeans(tail(getBiomass_ZooMizer(sim, zoo_params),250))
}

biom.df <- cbind(enviro, biomasses) %>% 
  rename(fish1=`1`,
         fish2=`2`,
         fish3=`3`,
         fish4=`4`,
         fish5=`5`) %>% 
  mutate(AllFish = fish1 + fish2 + fish3 + fish4 +fish5)

tiles <- foreach(i = 11:25) %dopar% {
  column <- sym(colnames(biom.df)[i])
  ggplot(biom.df, aes(x=sst, y = log10(chlo)))+
    geom_tile(aes(fill=log10(!!column)))+
    scale_fill_viridis_c()+
    labs(title = paste(colnames(biom.df)[i], "Biomass vs. SST and chlorophyll-a"))
  ggsave(filename = paste0("tileplot_",column,"_biomass_20210705.png"))
}

df$totalzoobiom <- rowSums(biom.df[,13:19])
tiles <- foreach(i = 13:19) %dopar% {
  column <- sym(colnames(df)[i])
  ggplot(biom.df, aes(x=sst, y = log10(chlo)))+
    geom_tile(aes(fill=(!!column)/totalzoobiom))+
    scale_fill_viridis_c()+
    labs(title = paste(colnames(biom.df)[i], "Biomass vs. SST and chlorophyll-a"))
  ggsave(filename = paste0("tileplot_",column,"_relbiom_20210705.png"))
}

abunds <- foreach(ID=1:nrow(enviro), .combine = rbind) %dopar% {
  source("ZooMizerSummaryFunctions.R")
  ID_char <- sprintf("%04d",ID)
  sim <- readRDS(paste0("Output/",jobname,"_ZooMizer_",ID_char,".RDS"))
  colMeans(tail(getAbundance_ZooMizer(sim, zoo_params),250))
}

abund.df <- cbind(enviro, abunds) %>% 
  rename(fish1=`1`,
         fish2=`2`,
         fish3=`3`,
         fish4=`4`,
         fish5=`5`) %>% 
  mutate(AllFish = fish1 + fish2 + fish3 + fish4 +fish5)

abund.df$totalzooabund <- rowSums(abund.df[, 13:19])

abundtiles <- foreach(i = 13:19) %dopar% {
  column <- sym(colnames(abund.df)[i])
  ggplot(abund.df, aes(x=sst, y = log10(chlo)))+
    geom_tile(aes(fill=(!!column)/totalzooabund))+
    scale_fill_viridis_c()+
    labs(title = paste(colnames(abund.df)[i], "Biomass vs. SST and chlorophyll-a"))
  ggsave(filename = paste0("tileplot_",column,"_relabund_20210705.png"))
}

df$totalzoobiom <- rowSums(biom.df[,13:19])
tiles <- foreach(i = 13:19) %dopar% {
  column <- sym(colnames(df)[i])
  ggplot(biom.df, aes(x=sst, y = log10(chlo)))+
    geom_tile(aes(fill=(!!column)/totalzoobiom))+
    scale_fill_viridis_c()+
    labs(title = paste(colnames(biom.df)[i], "Biomass vs. SST and chlorophyll-a"))
  ggsave(filename = paste0("tileplot_",column,"_relbiom_20210705.png"))
}
