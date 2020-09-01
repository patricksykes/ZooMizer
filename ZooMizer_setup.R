# Setup and helper functions for ZooMizer

library(devtools)
#most up to date master branch of mizer
#install_github("sizespectrum/mizer")
#install_github("astaaudzi/mizer-rewiring", ref = "temp-model-comp")
#documentation here:
#https://sizespectrum.org/mizer/dev/index.html
library(mizer)
require(tidyverse)

#remotes::install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)

setZooMizerConstants <- function(params, Groups, sst <- ){
  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  SearchVol <- getSearchVolume(params)
  M_sb <- getExtMort(mf.params)
  ZSpre <- 1 # senescence mortality prefactor
  ZSexp <- 0.3 # senescence mortality exponent
  
  pred_kernel <- getPredKernel(params)
  prey_weight_matrix <- matrix(params@w_full, nrow = length(params@w), ncol = length(params@w_full), byrow = TRUE)
  pred_weight_matrix <- matrix(params@w, nrow = length(params@w), ncol = length(params@w_full))
  
  for(i in 1:nrow(params@species_params)){
    ## Senescence mortality
    if(params@species_params$Type[i] == "Zooplankton"){
      M_sb[i,] <- ZSpre*(params@w/(10^(params@species_params$w_mat[i])))^ZSexp
      M_sb[i, 10^(params@species_params$w_max[i]) < params@w] <- 0
      M_sb[i, 10^(params@species_params$w_mat[i]) > params@w] <- 0
    }
    
    if(params@species_params$Type[i] == "Fish"){
      M_sb[i,] <- 0.1*ZSpre*(params@w/(10^(params@species_params$w_mat[i])))^ZSexp
      M_sb[i, 10^(params@species_params$w_max[i]) < par <- <- <- <- 
             ### Search volume
             SearchVol[i,] <- (params@species_params$gamma[i])*(params@w^(params@species_params$gamma[i]))
           SearchVol[i, 10^(params@species_params$w_max[i]) < params@w] <- 0
           SearchVol[i, 10^(params@species_params$w_min[i]) > params@w] <- 0
           
           ### Predation Kernels
           if(is.na(params@species_params$PPMRscale[i]) == FALSE){ # If group has an m-value (zooplankton)
             # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
             D.z <- 2*(3*params@w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
             betas <- (exp(0.02*log(D.z)^2 - params@species_params$PPMRscale[i] + 1.832))^3 # Wirtz's equation
             beta_mat <- matrix(betas, nrow = length(params@w), ncol = length(mf.params@w_full))
             
             # Calculate feeding kernels
             pred_kernel[i, , ] <- exp(-0.5*(log((beta_mat*prey_weight_matrix)/
                                                   pred_weight_matrix)/params@species_params$FeedWidth[i])^2)/
               sqrt(2*pi*params@species_params$FeedWidth[i]^2)
             # The feeding kernel of filter feeders is not expected to change  with increasing size so we fix it here
             
             # if (param$fixed_filterPPMR == TRUE){
             if(i == 3){
               pred_kernel[i, , ] <- matrix(pred_kernel[i,44,], nrow = length(params@w), ncol = length(mf.params@w_full), byrow = TRUE)
             }
             if(i == 8){
               pred_kernel[i, , ] <- matrix(pred_kernel[i,61,], nrow = length(params@w), ncol = length(mf.params@w_full), byrow = TRUE)
             }
             # }
             
           } else { # If group does not have an m-value (fish)
             beta_mat <- matrix(params@species_params$PPMR[i], nrow = length(params@w), ncol = length(params@w_full))
             
             # Calculate feeding kernels
             pred_kernel[i, , ] <- exp(-0.5*(log((beta_mat*prey_weight_matrix)/
                                                   pred_weight_matrix)/params@species_params$FeedWidth[i])^2)/
               sqrt(2*pi*params@species_params$FeedWidth[i]^2)
           }
           
    }
    
    #temperature effect 
    temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(params@species_params$species), ncol = length(params@w))
    M_sb <- temp_eff * M_sb # Incorporate temp effect on senscence mortality
    
    params <- setExtMort(params, z0 = M_sb)
    params <- setSearchVolume(params, SearchVol)
    params <- setPredKernel(params, pred_kernel)
    
    return(params)
  }
  
  
  new_project_simple <- function(params, n, n_pp, n_other, t, dt, steps, 
                                 effort, resource_dynamics_fn, other_dynamics_fns,
                                 rates_fns, ...) {    
    # Handy things
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    idx <- 2:no_w
    w_max_idx <- params@w_min_idx
    for (i in 1:length(w_max_idx)) {
      w_max_idx[i] <- which(round(log10(params@w),2)==round(log10(params@species_params$w_inf[i]),2)) 
    }
    
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    # This references the egg size bracket for all species, so for example
    # n[w_min_idx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    # Matrices for solver
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)
  
  for (i_time in 1:steps) {
    r <- rates_fns$Rates(
      params, n = n, n_pp = n_pp, n_other = n_other,
      t = t, effort = effort, rates_fns = rates_fns, ...)
    
    # Update time
    t <- t + dt
    
    # Update other components
    n_other_current <- n_other  # So that the resource dynamics can still 
    # use the current value
    for (component in names(params@other_dynamics)) {
      n_other[[component]] <-
        other_dynamics_fns[[component]](
          params,
          n = n,
          n_pp = n_pp,
          n_other = n_other_current,
          rates = r,
          t = t,
          dt = dt,
          component = component,
          ...
        )
    }
    
    # Update resource
    n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
                                 n_other = n_other_current, rates = r,
                                 t = t, dt = dt, ...)
    
    # Iterate species one time step forward:
    # a_{ij} = - g_i(w_{j-1}) / dw_j dt
    a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                      params@dw[idx], "/")
    # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
    b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + r$mort * dt
    # S_{ij} <- N_i(w_j)
    S[,idx] <- n[, idx, drop = FALSE]
    # Update n
    # for (i in 1:no_sp) # number of species assumed small, so no need to 
    #                      vectorize this loop over species
    #     for (j in (params@w_min_idx[i]+1):no_w)
    #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
    # This is implemented via Rcpp
    n <- inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                            A = a, B = b, S = S,
                            w_min_idx = params@w_min_idx)
  }
  
  # Update first and last size groups of n
  #TODO: Make this a little less hacky
  n[1,1] <- params@species_params$Prop[1]*n_pp[length(params@w_full)-length(params@w)+1]
  n[1,w_max_idx[1]] <- 0
  for (i in 2:no_sp) {
    n[i, w_max_idx[i]] <- 0
    if(params@species_params$Type[i] != "Fish"){
      n[i,params@w_min_idx[i]] <- params@species_params$Prop[i] * sum(n[-i,params@w_min_idx[i]])
    }
    if(params@species_params$Type[i] == "Fish"){
      n[i,params@w_min_idx[i]] <- 1/3 * sum(n[1:9,params@w_min_idx[i]])
    }
  }
  
  return(list(n = n, n_pp = n_pp, n_other = n_other, rates = r))
}

new_Encounter <- function(params, n, n_pp, n_other, t, ...) {
  
  # idx_sp are the index values of params@w_full such that
  # params@w_full[idx_sp] = params@w
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  
  # Note: removed the FFT code because it does not apply to this case.
  
  # If the feeding kernel does not have a fixed predator/prey mass ratio
  # then the integral is not a convolution integral and we can not use fft.
  # In this case we use the code from mizer version 0.3
  
  # n_eff_prey is the total prey abundance by size exposed to each
  # predator (prey not broken into species - here we are just working out
  # how much a predator eats - not which species are being eaten - that is
  # in the mortality calculation
  # \sum_j \theta_{ij} N_j(w_p) w_p dw_p
  assim_prey <- assim_eff * params@interaction
  n_eff_prey <- sweep(assim_prey %*% n, 2, 
                      params@w * params@dw, "*", check.margin = FALSE) 
  # pred_kernel is predator species x predator size x prey size
  # So multiply 3rd dimension of pred_kernel by the prey biomass density
  # Then sum over 3rd dimension to get consumption rate of each predator by 
  # predator size
  # This line is a bottle neck
  phi_prey_species <- rowSums(sweep(
    params@pred_kernel[, , idx_sp, drop = FALSE],
    c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
  # Eating the background
  # This line is a bottle neck
  phi_prey_background <- params@species_params$interaction_resource *
    rowSums(sweep(
      params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
      "*", check.margin = FALSE), dims = 2)
  encounter <- params@search_vol * (phi_prey_species + phi_prey_background)
  
  dimnames(encounter) <- dimnames(params@metab)
  
  # Add contributions from other components
  for (i in seq_along(params@other_encounter)) {
    encounter <- encounter + 
      do.call(params@other_encounter[[i]], 
              list(params = params,
                   n = n, n_pp = n_pp, n_other = n_other,
                   component = names(params@other_encounter)[[i]], ...))
  }
  
  return(encounter)
}
