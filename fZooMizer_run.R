phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa=params@resource_params$kappa, lambda=params@resource_params$lambda, ... ) {
  npp <- kappa*params@w_full^(1-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}

setZooMizerConstants <- function(params, Groups, sst){
  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  SearchVol <- getSearchVolume(params)
  M_sb <- getExtMort(params)
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
      M_sb[i, 10^(params@species_params$w_max[i]) < params@w] <- 0
      M_sb[i, 10^(params@species_params$w_mat[i]) > params@w] <- 0
    }
    
    ### Search volume
    SearchVol[i,] <- (params@species_params$gamma[i])*(params@w^(params@species_params$q[i]))
    SearchVol[i, 10^(params@species_params$w_max[i]) < params@w] <- 0
    SearchVol[i, 10^(params@species_params$w_min[i]) > params@w] <- 0
    
    ### Predation Kernels
    if(is.na(params@species_params$PPMRscale[i]) == FALSE){ # If group has an m-value (zooplankton)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
      D.z <- 2*(3*params@w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas <- (exp(0.02*log(D.z)^2 - params@species_params$PPMRscale[i] + 1.832))^3 # Wirtz's equation
      beta_mat <- matrix(betas, nrow = length(params@w), ncol = length(params@w_full))
      
      # Calculate feeding kernels
      pred_kernel[i, , ] <- exp(-0.5*(log((beta_mat*prey_weight_matrix)/
                                            pred_weight_matrix)/params@species_params$FeedWidth[i])^2)/
        sqrt(2*pi*params@species_params$FeedWidth[i]^2)
      # The feeding kernel of filter feeders is not expected to change  with increasing size so we fix it here
      
      # if (param$fixed_filterPPMR == TRUE){
      if(i == 3){
        pred_kernel[i, , ] <- matrix(pred_kernel[i,44,], nrow = length(params@w), ncol = length(params@w_full), byrow = TRUE)
      }
      if(i == 8){
        pred_kernel[i, , ] <- matrix(pred_kernel[i,61,], nrow = length(params@w), ncol = length(params@w_full), byrow = TRUE)
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
  
  
  params@initial_n_pp <- params@resource_params$kappa * params@w_full^(1-params@resource_params$lambda)
  params@initial_n_pp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  
  a_dynam <- (params@resource_params$kappa)*(params@w[1]^(2-params@resource_params$lambda)) # calculate coefficient for initial dynamic spectrum, so that N(w_phyto) equals N(w_dynam) at w[1]
  
  # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
  tempN <- matrix(a_dynam*(params@w)^-1, nrow = nrow(params@species_params), ncol = length(params@w), byrow = TRUE)
  props_z <- params@species_params$Prop[params@species_params$Type=="Zooplankton"] # Zooplankton proportions
  tempN[params@species_params$Type=="Zooplankton",] <- props_z * tempN[params@species_params$Type=="Zooplankton",] # Set abundances of diff zoo groups based on smallest size class proportions
  tempN[params@species_params$Type=="Fish",] <- (1/sum(mf.params@species_params$Type=="Fish")) * tempN[params@species_params$Type=="Fish",] # Set abundandances of fish groups based on smallest size class proportions
  
  # For each group, set densities at w > Winf and w < Wmin to 0
  tempN[unlist(tapply(round(log10(param$w), digits = 2), 1:length(param$w), function(wx,Winf) Winf < wx, Winf = param$Groups$Wmax))] <- 0
  tempN[unlist(tapply(round(log10(param$w), digits = 2), 1:length(param$w), function(wx,Wmin) Wmin > wx, Wmin = param$Groups$W0))] <- 0
  params@initial_n <- tempN
  
  params <- setExtMort(params, z0 = M_sb)
  params <- setSearchVolume(params, SearchVol)
  params <- setPredKernel(params, pred_kernel)
  
  return(params)
}

setassim_eff <- function(groups){
assim_eff = matrix(groups$GrossGEscale * groups$Carbon, nrow = nrow(groups), ncol = nrow(groups))
get_filterfeeders <- which(groups$FeedType == "FilterFeeder")

for (i in get_filterfeeders) {
  assim_eff[,i] <- assim_eff[,i] / groups$Carbon[i]
}
return(assim_eff)
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
  assim_prey <- params@other_params$assim_eff * params@interaction
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
  phi_prey_background <- assim_phyto * params@species_params$interaction_resource *
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

new_PredRate <- function (params, n, n_pp, n_other, t, feeding_level, ...) 
{
  no_sp <- dim(params@interaction)[1]
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  
  if (length(params@ft_pred_kernel_p) == 1) {
    n_total_in_size_bins <- sweep(n, 2, params@dw, "*", 
                                  check.margin = FALSE)
    pred_rate <- sweep(params@pred_kernel, c(1, 2), (1 - 
                                                       feeding_level) * params@search_vol * n_total_in_size_bins * params@other_params$temp_eff, 
                       "*", check.margin = FALSE)
    pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
    return(pred_rate)
  }
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
  Q[, idx_sp] <- sweep((1 - feeding_level) * params@search_vol * params@other_params$temp_eff *
                         n, 2, params@dw, "*")
  pred_rate <- Re(t(mvfft(t(params@ft_pred_kernel_p) * mvfft(t(Q)), 
                          inverse = TRUE)))/no_w_full
  return(pred_rate * params@ft_mask)
}

new_EReproAndGrowth <- function (params, n, n_pp, n_other, t, encounter, feeding_level, ...) 
{
  sweep((1 - feeding_level) * encounter * params@other_params$temp_eff, 1, params@species_params$alpha, 
        "*", check.margin = FALSE) - params@metab
}

newFeedingLevel <- function (params, n, n_pp, n_other, t, encounter, ...) 
{
  return(encounter)
}

fZooMizer_run <- function(groups, input){

  kappa = 10^(input$phyto_int)
  lambda = 1-input$phyto_slope
  chlo = input$chlo
  sst = input$sst
  dt = input$dt
  tmaxx = input$tmaxx
  
  #data
  groups$w_min <- 10^groups$w_min #convert from log10 values
  groups$w_inf <- 10^groups$w_inf
  groups$w_mat <- 10^groups$w_mat
  groups$h <- 1e50 # should be Inf, but that breaks the calculations. Massive value still works out to effectively unlimited feeding as allowed in ZooMSS
  groups$ks <- 0 #turn off standard metabolism
  #todo - ramp up constant repro for coexistence
  
  mf.params <- newMultispeciesParams(species_params=groups,
                                     interaction=NULL, #NULL sets all to 1, no strict herbivores
                                     no_w = 178, #number of zoo+fish size classes;
                                     min_w_pp = 10^(-14.4), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
                                     w_pp_cutoff = 10^(input$phyto_max), #maximum phyto size
                                     n = 0.7, #The allometric growth exponent used in ZooMSS
                                     z0pre = 1, #external mortality (senescence)
                                     z0exp = 0.3,
                                     resource_dynamics = "phyto_fixed",
                                     kappa = kappa, 
                                     lambda = lambda,
                                     RDD = constantRDD(species_params = groups) #first go at this
                                     #pred_kernel = ... #probably easiest to just import this/pre-calculate it, once dimensions are worked out
  )
  
  temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  
  mf.params@other_params$assim_eff <- setassim_eff(groups)
  
  mf.params@other_params$temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  
  mf.params <- setZooMizerConstants(params = mf.params, Groups = groups, sst = input$sst)
  
  mf.params <- setParams(mf.params)
  
  mf.params <- setRateFunction(mf.params, "PredRate", "new_PredRate")
  mf.params <- setRateFunction(mf.params, "EReproAndGrowth", "new_EReproAndGrowth")
  mf.params <- setRateFunction(mf.params, "FeedingLevel", "newFeedingLevel")
  mf.params <- setRateFunction(mf.params, "Encounter", "new_Encounter")
  
  mf.params <- setReproduction(mf.params, repro_prop = matrix(0, nrow = nrow(mf.params@psi), ncol = ncol(mf.params@psi)))
  
  sim <- project(mf.params, dt = dt, t_max = tmaxx)
  
  return(sim)
}
