# Code for analysis of slope sensitivity to choice of endpoints
# Patrick Sykes
# May 2021


library(mizer)
library(tidyverse)

# data
enviro <- readRDS("data/enviro_test20.RDS")
zoomizergrid <- readRDS("test_grid_20210317.RDS")

for(ID in 1:length(zoomizergrid)) zoomizergrid[[ID]]@n <- sweep(zoomizergrid[[ID]]@n, 2, 2.5 * zoomizergrid[[ID]]@params@species_params$Carbon, "*")


# choices for wmin/wmax
wmins <- c(10^(-10.7), 10^seq(from = -10.5, to = -6, by = 0.5))
wmaxs <- 10^seq(from = -5, to = 2, by = 0.5)


#set up to run parallel
library(doParallel)
cl <- makePSOCKcluster(max(1, detectCores()-1))     ## set up cores-1 machines
registerDoParallel(cl, cores = (max(1, detectCores()-1)))
clusterEvalQ(cl, library(mizer)) %>% invisible()
clusterExport(cl, c("df", "enviro", "zoomizergrid"))



Res3 <- foreach(ID = 1:length(zoomizergrid), .combine = rbind) %:%
  foreach(wmin = wmins, .combine = rbind) %dopar% {
    sst <- enviro$sst[ID]
    chlo <- enviro$chlo[ID]
    slopes <- getCommunitySlope(zoomizergrid[[ID]], min_w = wmin, max_w = 1e-3, species = 1:9)
    slopes <- slopes[501:1000,]
    slopes <- colMeans(slopes)
    c("ID" = ID,"chlo" = chlo, "sst" = sst,"min_w" = wmin,"max_w" = wmax, slopes)
  }
saveRDS(Res3, "zooplanktonslopes-3.rds")  

Res5 <- foreach(ID = 1:length(zoomizergrid), .combine = rbind) %:%
  foreach(wmin = wmins, .combine = rbind) %dopar% {
    sst <- enviro$sst[ID]
    chlo <- enviro$chlo[ID]
    slopes <- getCommunitySlope(zoomizergrid[[ID]], min_w = wmin, max_w = 1e-5, species = 1:9)
    slopes <- slopes[501:1000,]
    slopes <- colMeans(slopes)
    c("ID" = ID,"chlo" = chlo, "sst" = sst,"min_w" = wmin,"max_w" = wmax, slopes)
  }
saveRDS(Res5, "zooplanktonslopes-5.rds")  

stopCluster(cl)
rm(cl)


Res3.df <- as.data.frame(Res3)
Res5.df <- as.data.frame(Res5)
library(ggplot2)

# ggplot(Res3.df, aes(group = ID,x=ID, y = slope))+geom_boxplot()+stat_summary(fun=mean, colour="red", geom="point")
# ggplot(Res5.df, aes(group = ID,x=ID, y = slope))+geom_boxplot()+stat_summary(fun=mean, colour="red", geom="point")

ggplot(Res3.df, aes(x = log10(min_w), group = ID, colour = ID))+geom_point(aes(y=slope))
ggplot(Res3.df, aes(x = log10(min_w), group = ID, colour = ID))+geom_point(aes(y=intercept))

ggplot(Res5.df, aes(x = log10(min_w), group = ID, colour = ID))+geom_point(aes(y=slope))
ggplot(Res5.df, aes(x = log10(min_w), group = ID, colour = ID))+geom_point(aes(y=intercept))

Res3.df <- Res3.df %>% filter(min_w == 10^(-10.7))
Res5.df <- Res5.df %>% filter(min_w == 10^(-10.7))

Res5.df
