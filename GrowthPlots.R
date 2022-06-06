library(ggplot2)
library(patchwork)
# library(egg)
# library(rgdal)
# library(raster)
# library(reshape2)

# ## IMPORT MODEL SIMULATIONS, AND ENVIRONMENTAL DATA
# load('17_09_2019_01/res.RData')
# N <- res
# load('17_09_2019_01/growth.RData')
# load('17_09_2019_01/add_objects.RData')
# enviro_data <- readRDS('17_09_2019_01/params.RDS')
# w = 10^(seq(from = -10.7, to =  7, 0.1))

### MULTIPLY LARVS AND SALPS BY 15 (TO INCLUDE HIGH CARBON EXPORT BY THESE GROUPS),
### STANDARDISE TO 15C, THEN CALCULATE MEAN
#gg_all <- 0

#for(i in 1:dim(enviro_data)[1]){ # There would be a much smarter way to do this...
#  gg_all <- gg_all + growth[[i]]*(2^((15-30)/10)/2^((enviro_data[i,'sst'] - 30)/10))
#}

#gg_all <- gg_all/dim(enviro_data)[1]

## SST = 26, Chlo = 3.34 mg m-3. I used the cell with the highest chlo and sst values, to
## make sure zoo growth rates did not exceed empirical observations in the absence of phytoplankton
## dynamics

plotZooMizerGrowthsComparison <- function(out){

params <- setInitialValues(out@params, out)
times <- max(getTimes(out))
params@initial_n[] <- colMeans(out@n[ceiling(times/2+1):times,,])
nothers <- array(0,c(times,dim(params@initial_n_other[[1]])))
for (i in 1:times) {
  nothers[i,,]  <- out@n_other[[i]]
}
params@initial_n_other$zoo[] <- colMeans(nothers[ceiling(times/2+1):times,,])
params@other_params$zoo$params@initial_n_other[] <- colMeans(nothers[ceiling(times/2+1):times,,])
params@initial_n_pp <- colMeans(out@n_pp[ceiling(times/2+1):times,])

growth <- getEGrowth(params@other_params$zoo$params)
for(i in 1:9) growth[i,1:params@other_params$zoo$params@w_min_idx[i]] <- 0
w <- params@other_params$zoo$params@w
# for(i in 1:9) growth[i,1:params[[2]]@other_params$zoo$params@w_min_idx[i]] <- 0
# w <- params[[2]]@other_params$zoo$params@w

gg_all <- growth*(2^((15-30)/10)/2^((enviro[2,'sst'] - 30)/10))

gg <- gg_all/365 # Growth rates in years from model, convert to days
gg[c(3,8),] <- 15*gg[c(3,8),] # Multiply larvs and salps by 15 ()

# par(mfrow = c(2,2))
# cols = rainbow(9)
### FLAGELLATES AND CILIATES
flag_gs <- read.table('data/FlagHirst.txt', header = TRUE)
flag_gs$BodyWeight <- (flag_gs$BodyWeight/1e6)*(1/0.15)
flag_gs$Growth <- flag_gs$Growth*24

cil_gs <- read.table('data/CilHirst.txt', header = TRUE)
cil_gs$BodyWeight <- (cil_gs$BodyWeight/1e6)*(1/0.15)
cil_gs$Growth <- cil_gs$Growth*24

# plot(log10(flag_gs$BodyWeight), log10(flag_gs$Growth), xlim = c(-10.8, -6), ylim = c(-0.7, 0.8), col = 'blue',
#      ylab = expression(paste('Growth (d'^{-1},')', sep = '')),
#      xlab = expression(paste('log'[10],'(Body Size, g)', sep = '')), main = 'Microzoo (Cils = blue, Flag = red)')
# points(log10(cil_gs$BodyWeight), log10(cil_gs$Growth), col = 'red')
# 
# lines(log10(w), log10(gg[1,]/w), col = 'blue', lwd = 2)
# lines(log10(w), log10(gg[2,]/w), col = 'red', lwd = 2)

## Microzoo ggplot
model_microzoo <- data.frame('w' = log10(w), 'cils' = log10(gg[1,]/w), 'flags' = log10(gg[2,]/w))
model_microzoo[which(model_microzoo$cils == -Inf), 'cils'] <- NA
model_microzoo[which(model_microzoo$flags == -Inf), 'flags'] <- NA

cil_gs$BodyWeight <- log10(cil_gs$BodyWeight)
cil_gs$Growth <- log10(cil_gs$Growth)
flag_gs$BodyWeight <- log10(flag_gs$BodyWeight)
flag_gs$Growth <- log10(flag_gs$Growth)


microzoo_plot <- ggplot() + 
  geom_line(data = model_microzoo, aes(x = w, y = cils), color = 'red', size = 1.2)+
  geom_line(data = model_microzoo, aes(x = w, y = flags), color = 'blue', size = 1.2)+
  geom_point(data = flag_gs, aes(x = BodyWeight, y = Growth), color = 'red')+
  geom_point(data = cil_gs, aes(x = BodyWeight, y = Growth), color = 'blue')+
  theme_bw() + xlim(-10.5,-6) + theme(plot.subtitle=element_text(size=12, color="black"),
                          axis.text = element_text(size = 12),
                          axis.title.x = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                          axis.title.y = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                          legend.position = 'top')+xlab(expression(paste("log"[10], "(Body Size, g)")))+ 
  ylab(expression(paste("Growth (d"^{-1}, ')')))+scale_linetype_manual(values = c("solid","dashed"),
name="Equations",labels=c("Equation 1","Equation 2"),guide="legend")+ 
  ggtitle(label = "",subtitle = 'a) Microzooplankton - flagellates blue, ciliates red') 
  

### SALPS AND LARVACEANS
larv_gs <- read.table('data/LarvHirst.txt', header = TRUE)
larv_gs$BodyWeight <- (larv_gs$BodyWeight/1e6)*100
larv_gs[is.na(larv_gs$Temp), 'Temp'] <- 15
larv_gs$Growth_TempCorrect <- larv_gs$Growth*exp(-0.65/(8.62e-5*(273.15+15)))/exp(-0.65/(8.62e-5*(273.15+larv_gs$Temp)))

salp_gs <- read.table('data/SalpHirst.txt', header = TRUE)
salp_gs$BodyWeight <- (salp_gs$BodyWeight/1e6)*100
salp_gs[is.na(salp_gs$Temp), 'Temp'] <- 15
salp_gs$Growth_TempCorrect <- salp_gs$Growth*exp(-0.65/(8.62e-5*(273.15+15)))/exp(-0.65/(8.62e-5*(273.15+salp_gs$Temp)))

# plot(log10(larv_gs$BodyWeight), log10(larv_gs$Growth), xlim = c(-7,1), ylim = c(-2,0.5), col = 'blue',
#      ylab = expression(paste('Growth (d'^{-1},')', sep = '')),
#      xlab = expression(paste('log'[10],'(Body Size, g)', sep = '')),
#      main = 'Filter Feeders (Larvs = blue, \n Salps = red')
# points(log10(salp_gs$BodyWeight), log10(salp_gs$Growth), col = 'red')
# 
# lines(log10(w), log10(gg[3,]/w), col = 'blue', lwd = 2)
# lines(log10(w), log10(gg[8,]/w), col = 'red', lwd = 2)

## Filter feeder ggplot
model_filters <- data.frame('w' = log10(w), 'larvs' = log10(gg[3,]/w), 'salps' = log10(gg[8,]/w))
model_filters[which(model_filters$larvs == -Inf), 'larvs'] <- NA
model_filters[which(model_filters$salps == -Inf), 'salps'] <- NA

salp_gs$BodyWeight <- log10(salp_gs$BodyWeight)
salp_gs$Growth_TempCorrect <- log10(salp_gs$Growth_TempCorrect)
larv_gs$BodyWeight <- log10(larv_gs$BodyWeight)
larv_gs$Growth_TempCorrect <- log10(larv_gs$Growth_TempCorrect)

filter_feeder_plot <- ggplot() + 
  geom_line(data = model_filters, aes(x = w, y = salps), color = 'green', size = 1.2)+
  geom_line(data = model_filters, aes(x = w, y = larvs), color = 'purple', size = 1.2)+
  geom_point(data = salp_gs, aes(x = BodyWeight, y = Growth_TempCorrect), color = 'green')+
  geom_point(data = larv_gs, aes(x = BodyWeight, y = Growth_TempCorrect), color = 'purple')+
  theme_bw() + xlim(-6.5, 1) + theme(plot.subtitle=element_text(size=12, color="black"),
                                      axis.text = element_text(size = 12),
                                      axis.title.x = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                      axis.title.y = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                      legend.position = 'top')+xlab(expression(paste("log"[10], "(Body Size, g)")))+ 
  ylab(expression(paste("Growth (d"^{-1}, ')')))+scale_linetype_manual(values = c("solid","dashed"),
                                                                       name="Equations",labels=c("Equation 1","Equation 2"),guide="legend")+ 
  ggtitle(label = "",subtitle = 'b) Filter Feeders - salps green, larvaceans purple') 

### OMNIVOROUS COPEPODS AND EUPHAUSIIDS
cops_gs <- read.table('data/CCopepodHirst.txt', header = TRUE)
ocops_gs <- cops_gs[cops_gs$Diet == 'O',]
ocops_gs$BodyWeight <- (as.numeric(as.character(ocops_gs$BodyWeight))/1e6)*(1/0.12)
ocops_gs[is.na(ocops_gs$Temp), 'Temp'] <- 15
ocops_gs$Growth_TempCorrect <- as.numeric(as.character(ocops_gs$g))*exp(-0.65/(8.62e-5*(273.15+15)))/exp(-0.65/(8.62e-5*(273.15+ocops_gs$Temp)))

euphs_gs <- read.table('data/KrillHirst.txt', header = TRUE)
euphs_gs$BodyWeight <- (euphs_gs$BodyWeight/1e6)*(1/0.12)
euphs_gs[is.na(euphs_gs$Temp), 'Temp'] <- 15
euphs_gs$Growth_TempCorrect <- euphs_gs$Growth*exp(-0.65/(8.62e-5*(273.15+15)))/exp(-0.65/(8.62e-5*(273.15+euphs_gs$Temp)))

# plot(log10(ocops_gs$BodyWeight), log10(ocops_gs$Growth_TempCorrect), xlim = c(-7.5, 0), ylim = c(-3,0.2), col = 'blue', 
#      ylab = expression(paste('Growth (d'^{-1},')', sep = '')),
#      xlab = expression(paste('log'[10],'(Body Size, g)', sep = '')), main = 'Omnivores (Cops = blue, Euphs = red)')
# points(log10(euphs_gs$BodyWeight), log10(euphs_gs$Growth_TempCorrect), col = 'red')
# 
# lines(log10(w), log10(gg[4,]/w), col = 'black', lwd = 2)
# lines(log10(w), log10(gg[6,]/w), col = 'red', lwd = 2)

## Omnivore ggplot
model_omnis <- data.frame('w' = log10(w), 'omnis' = log10(gg[4,]/w), 'euphs' = log10(gg[6,]/w))
model_omnis[which(model_omnis$omnis == -Inf), 'omnis'] <- NA
model_omnis[which(model_omnis$euphs == -Inf), 'euphs'] <- NA

ocops_gs$BodyWeight <- log10(ocops_gs$BodyWeight)
ocops_gs$Growth_TempCorrect <- log10(ocops_gs$Growth_TempCorrect)
ocops_gs[which(ocops_gs$Growth_TempCorrect < -3), 'Growth_TempCorrect'] <- NA

euphs_gs$BodyWeight <- log10(euphs_gs$BodyWeight)
euphs_gs$Growth_TempCorrect <- log10(euphs_gs$Growth_TempCorrect)
euphs_gs[which(euphs_gs$Growth_TempCorrect < -3), 'Growth_TempCorrect'] <- NA

omni_plot <- ggplot() + 
  geom_point(data = ocops_gs, aes(x = BodyWeight, y = Growth_TempCorrect), color = 'cyan3')+
  geom_line(data = model_omnis, aes(x = w, y = omnis), color = 'cyan2', size = 1.2)+
  geom_line(data = model_omnis, aes(x = w, y = euphs), color = 'cornflowerblue', size = 1.2)+
  geom_point(data = euphs_gs, aes(x = BodyWeight, y = Growth_TempCorrect), color = 'cornflowerblue')+
  theme_bw() + xlim(-6.4, 0.1) + ylim(-3.8, 0.1) + theme(plot.subtitle=element_text(size=12, color="black"),
                                     axis.text = element_text(size = 12),
                                     axis.title.x = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                     axis.title.y = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                     legend.position = 'top')+xlab(expression(paste("log"[10], "(Body Size, g)")))+ 
  ylab(expression(paste("Growth (d"^{-1}, ')')))+scale_linetype_manual(values = c("solid","dashed"),
                                                                       name="Equations",labels=c("Equation 1","Equation 2"),guide="legend")+ 
  ggtitle(label = "",subtitle = 'c) Omnivores - Omni. Cops cyan, Krill blue') 

### CARN COPS, CHAETOGNATHS AND JELLYFISH
cops_gs <- read.table('data/CCopepodHirst.txt', header = TRUE)
ccops_gs <- cops_gs[cops_gs$Diet == 'C',]
ccops_gs$BodyWeight <- (as.numeric(as.character(ccops_gs$BodyWeight))/1e6)*(1/0.12)
ccops_gs[is.na(ccops_gs$Temp), 'Temp'] <- 15
ccops_gs$Growth_TempCorrect <- as.numeric(as.character(ccops_gs$g))*exp(-0.65/(8.62e-5*(273.15+15)))/exp(-0.65/(8.62e-5*(273.15+ccops_gs$Temp)))

jelly_gs <- read.table('data/CnidariaHirst.txt', header = TRUE)
jelly_gs$BodyWeight <- (jelly_gs$BodyWeight/1e6)*(1/0.005)
jelly_gs$Growth <- jelly_gs$Growth

chaet_gs <- read.table('data/ChaetHirst.txt', header = TRUE)
chaet_gs$BodyWeight <- (chaet_gs$BodyWeight/1e6)*(1/0.02)
chaet_gs$Growth <- chaet_gs$Growth

# plot(log10(ccops_gs$BodyWeight), log10(ccops_gs$Growth), xlim = c(-7.5, 2), ylim = c(-3.1,0.5),col = 'blue',
#      ylab = expression(paste('Growth (d'^{-1},')', sep = '')),
#      xlab = expression(paste('log'[10],'(Body Size, g)', sep = '')),
#      main = 'Carnivores (Cops = blue, \n Chaets = red, Jellys = green)')
# points(log10(chaet_gs$BodyWeight), log10(chaet_gs$Growth), col = 'red')
# points(log10(jelly_gs$BodyWeight), log10(jelly_gs$Growth), col = 'green')
# 
# lines(log10(w), log10(gg[5,]/w), col = 'blue', lwd = 2)
# lines(log10(w), log10(gg[7,]/w), col = 'red', lwd = 2)
# lines(log10(w), log10(gg[9,]/w), col = 'green', lwd = 2)

## GGPLOT
model_carns <- data.frame('w' = log10(w), 'carns' = log10(gg[5,]/w), 'chaets' = log10(gg[7,]/w), 'jellys' = log10(gg[9,]/w))
model_carns[which(model_carns$carns == -Inf), 'carns'] <- NA
model_carns[which(model_carns$chaets == -Inf), 'chaets'] <- NA
model_carns[which(model_carns$jellys == -Inf), 'jellys'] <- NA

ccops_gs$BodyWeight <- log10(ccops_gs$BodyWeight)
ccops_gs$Growth_TempCorrect <- log10(ccops_gs$Growth_TempCorrect)
ccops_gs[which(ccops_gs$Growth_TempCorrect < -3), 'Growth_TempCorrect'] <- NA

chaet_gs$BodyWeight <- log10(chaet_gs$BodyWeight)
chaet_gs$Growth <- log10(chaet_gs$Growth)

jelly_gs$BodyWeight <- log10(jelly_gs$BodyWeight)
jelly_gs$Growth <- log10(jelly_gs$Growth)

carn_plot <- ggplot() + 
  geom_line(data = model_carns, aes(x = w, y = carns), color = 'orange', size = 1.2)+
  geom_point(data = ccops_gs, aes(x = BodyWeight, y = Growth_TempCorrect), color = 'orange')+
  geom_line(data = model_carns, aes(x = w, y = chaets), color = 'darkgreen', size = 1.2)+
  geom_point(data = chaet_gs, aes(x = BodyWeight, y = Growth), color = 'darkgreen')+
  geom_line(data = model_carns, aes(x = w, y = jellys), color = 'black', size = 1.2)+
  geom_point(data = jelly_gs, aes(x = BodyWeight, y = Growth), color = 'black')+
  theme_bw() + xlim(-7, 2) + ylim(-3.5, 0) + theme(plot.subtitle=element_text(size=12, color="black"),
                                                         axis.text = element_text(size = 12),
                                                         axis.title.x = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                                         axis.title.y = element_text(size =12, margin = unit(c(0.2,0.2,0,0), unit = "cm")),
                                                         legend.position = 'top')+xlab(expression(paste("log"[10], "(Body Size, g)")))+ 
  ylab(expression(paste("Growth (d"^{-1}, ')')))+
  ggtitle(label = "",subtitle = 'd) Carnivores - Carn.Cops orange, Chaet. blue, Jellyfish black') 



plotgg <- (microzoo_plot + filter_feeder_plot) / (omni_plot + carn_plot) + plot_layout(guides = "collect")

return(plotgg)
}
# plot_storage <- list(microzoo_plot, filter_feeder_plot, omni_plot, carn_plot)
# ggsave(filename = "growth_figures.png", plot = ggarrange(plots = plot_storage, nrow = 2), width = 7.5, height = 5.5)
# 
# 

