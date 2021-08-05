phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa=params@resource_params$kappa, lambda=params@resource_params$lambda, ... ) {
  npp <- kappa*params@w_full^(1-lambda) / params@dw_full #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff* (1 - 1e-06)] <- 0
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

  for (i in 1:nrow(params@species_params)) {
    ## Senescence mortality
    if (params@species_params$Type[i] == "Zooplankton") {
      M_sb[i,] <- ZSpre*(params@w/(params@species_params$w_mat[i]))^ZSexp
     # M_sb[i, params@species_params$w_inf[i] < params@w * (1 + 1e-06)] <- 0 # don't care about senescence of size classes not in the model anyway
      M_sb[i, params@species_params$w_mat[i] > params@w * (1 + 1e-06)] <- 0 #turn off "senescence" for smaller size classes
    }

    if (params@species_params$Type[i] == "Fish") {
      M_sb[i,] <- 0.1*ZSpre*(params@w/(params@species_params$w_mat[i]))^ZSexp
    #  M_sb[i, params@species_params$w_inf[i] < params@w * (1 + 1e-06)] <- 0 # don't care about senescence of size classes not in the model anyway
      M_sb[i, params@species_params$w_mat[i] > params@w * (1 + 1e-06)] <- 0
    }

    ### Search volume
    SearchVol[i,] <- (params@species_params$gamma[i])*(params@w^(params@species_params$q[i]))
    # SearchVol[i, params@species_params$w_inf[i] < params@w * (1 + 1e-06)] <- 0
    # SearchVol[i, params@species_params$w_min[i] > params@w * (1 + 1e-06)] <- 0

    ### Predation Kernels
    if (is.na(params@species_params$PPMRscale[i]) == FALSE){ # If group has an m-value (zooplankton)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)
      D.z <- 2*(3*params@w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas <- (exp(0.02*log(D.z)^2 - params@species_params$PPMRscale[i] + 1.832))^3 # Wirtz's equation
      beta_mat <- matrix(betas, nrow = length(params@w), ncol = length(params@w_full))

      # Calculate feeding kernels
      pred_kernel[i, , ] <- exp(-0.5*(log((beta_mat*prey_weight_matrix) /
                                            pred_weight_matrix)/params@species_params$FeedWidth[i])^2) /
        sqrt(2*pi*params@species_params$FeedWidth[i]^2)
      # The feeding kernel of filter feeders is not expected to change  with increasing size so we fix it here

      # if (param$fixed_filterPPMR == TRUE){
      if (Groups$FeedType[i] == "FilterFeeder") {
        pred_kernel[i, , ] <- matrix(pred_kernel[i,params@w_min_idx[i],], nrow = length(params@w), ncol = length(params@w_full), byrow = TRUE)
      }
          # }

    } else { # If group does not have an m-value (fish)
      beta_mat <- matrix(params@species_params$PPMR[i], nrow = length(params@w), ncol = length(params@w_full))

      # Calculate feeding kernels
      pred_kernel[i, , ] <- exp(-0.5*(log((beta_mat*prey_weight_matrix) /
                                            pred_weight_matrix) / params@species_params$FeedWidth[i])^2) /
        sqrt(2*pi*params@species_params$FeedWidth[i]^2)
    }

  }
  # SearchVol[12,178] <- (params@species_params$gamma[12])*(params@w[178]^(params@species_params$q[12])) #adding last size class by hand

  #temperature effect
  M_sb <- params@other_params$temp_eff * M_sb * 10 # Incorporate temp effect on senscence mortality


  params@initial_n_pp <- params@resource_params$kappa * params@w_full^(1 - params@resource_params$lambda)/params@dw_full
  params@initial_n_pp[params@w_full > params@resource_params$w_pp_cutoff] <- 0


  a_dynam <- (params@resource_params$kappa)*(params@w[1]^(2 - params@resource_params$lambda))#/params@dw[1] # calculate coefficient for initial dynamic spectrum, so that N(w_phyto) equals N(w_dynam) at w[1]

  # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
  tempN <- matrix(a_dynam*(params@w)^(-1)/params@dw, nrow = nrow(params@species_params), ncol = length(params@w), byrow = TRUE)
  props_z <- params@species_params$Prop[params@species_params$Type == "Zooplankton"] # Zooplankton proportions
  tempN[params@species_params$Type == "Zooplankton",] <- props_z * tempN[params@species_params$Type == "Zooplankton",] # Set abundances of diff zoo groups based on smallest size class proportions
  tempN[params@species_params$Type == "Fish",] <- (1/sum(params@species_params$Type == "Fish")) * tempN[params@species_params$Type=="Fish",] # Set abundandances of fish groups based on smallest size class proportions

  # For each group, set densities at w > Winf and w < Wmin to 0
  params@species_params$w_min <- params@w[params@w_min_idx]
  tempN[unlist(tapply(round(log10(params@w), digits = 2), 1:length(params@w), function(wx,Winf) Winf < wx, Winf = log10(params@species_params$w_inf)))] <- 0
  tempN[unlist(tapply(params@w, 1:length(params@w), function(wx,Wmin) Wmin > wx, Wmin = params@species_params$w_min))] <- 0
  #dimnames(tempN) <- dimnames(params@initial_n)
  params@initial_n[] <- tempN
  
  #SearchVol <- readRDS("data/SearchVol.rds")
  
  params <- setExtMort(params, z0 = M_sb)
  params <- setSearchVolume(params, search_vol = SearchVol)
  params <- setPredKernel(params, pred_kernel)

  return(params)
}

setassim_eff <- function(groups){
  assim_eff = matrix(groups$GrossGEscale * groups$Carbon, nrow = nrow(groups), ncol = nrow(groups))
  #get_filterfeeders <- which(groups$FeedType == "FilterFeeder")

  #for (i in get_filterfeeders) {
  #  assim_eff[,i] <- assim_eff[,i] / groups$Carbon[i]
  #}
  return(t(assim_eff))
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
    w_max_idx[i] <- which(round(log10(params@w),2) == round(log10(params@species_params$w_inf[i]),2))
  }
  
  fish_grps <- which(params@species_params$Type == "Fish")
  
  if(sum(params@species_params$Type == "Zooplankton") > 1){ # If there's only one zoo group, then you do not need w0idx. All this stuff gives you info about all zoo groups except the smallest zoo group.
    w0idx <- which(params@w_min_idx > min(params@w_min_idx) & is.na(params@species_params$Prop) == FALSE)
    w0mins <- params@w_min_idx[w0idx]
    props_z <- params@species_params$Prop[w0idx] # Zooplankton proportions
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
    # for (component in names(params@other_dynamics)) {
    #   n_other[[component]] <-
    #     other_dynamics_fns[[component]](
    #       params,
    #       n = n,
    #       n_pp = n_pp,
    #       n_other = n_other_current,
    #       rates = r,
    #       t = t,
    #       dt = dt,
    #       component = component,
    #       ...
    #     )
    # }

    # Update resource
    # n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
    #                              n_other = n_other_current, rates = r,
    #                              t = t, dt = dt, ...)
    n_pp <- n_pp
    # Iterate species one time step forward:
    # a_{ij} = - g_i(w_{j-1}) / dw_j dt
    a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                      params@w[idx-1], "/") * 10^(-0.1) / log(10)
    # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
    b[] <- 1 + sweep(r$e_growth * dt , 2, params@w, "/") /log(10) + r$mort * dt * 0.1
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
  

  # Update first and last size groups of n
  #TODO: Make this a little less hacky
  #n[1,1] <- params@species_params$Prop[1]*n_pp[length(params@w_full)-length(params@w)+1]
  if(sum(params@species_params$Type == "Zooplankton") > 1){ # If you only have one zoo group, it will be locked to phyto spectrum so you do not need to do this
    for(i in 1:length(w0idx)){
      w_min_curr <- w0mins[i]
      exclude_mins <- w0idx[which(w0mins == w_min_curr)]
      n[w0idx[i], w_min_curr] <- props_z[i] * sum(n[-exclude_mins, w_min_curr])
    }
  }
  
  fish_mins <- unlist(params@w_min_idx[params@species_params$Type == "Fish"])
  
  
  if(sum(params@species_params$Type == "Fish") > 1 && sum(params@species_params$Type == "Zooplankton") > 1){
    n[fish_grps,fish_mins] <- (1/length(fish_grps))*(colSums(n[-fish_grps,fish_mins]))
  }else{
    n[fish_grps, fish_mins] <- (1/length(fish_grps))*sum(n[-fish_grps, fish_mins])
  }

  for (i in 1:no_sp) {
    n[i, which(1:no_w >= w_max_idx[i])] <- 0
  }
  
    
  # n[1,1] <- params@initial_n[1,1]
  # n[1,w_max_idx[1]] <- 0
  # for (i in 2:no_sp) {
  #   n[i, w_max_idx[i]] <- 0
  #   if(params@species_params$Type[i] != "Fish"){
  #     n[i,params@w_min_idx[i]] <- params@species_params$Prop[i] * sum(n[-i,params@w_min_idx[i]])
  #   }
  #   if(params@species_params$Type[i] == "Fish"){
  #     n[i,params@w_min_idx[i]] <- 1/length(which(params@species_params$Type=="Fish")) * sum(n[which(params@species_params$Type!="Fish"),params@w_min_idx[i]])
  #   }
  # }
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
  phi_prey_background <- params@other_params$assim_phyto * params@species_params$interaction_resource *
    rowSums(sweep(
      params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
      "*", check.margin = FALSE), dims = 2)
  encounter <- params@other_params$temp_eff * params@search_vol * (phi_prey_species + phi_prey_background)
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

new_PredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...)
{
  n_total_in_size_bins <- sweep(n, 2, params@dw, "*",
                                check.margin = FALSE)
  pred_rate <- sweep(params@pred_kernel, c(1, 2), (1 - feeding_level) * params@other_params$temp_eff * params@search_vol * n_total_in_size_bins,
                     "*", check.margin = FALSE)
  pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
  return(pred_rate)
}

new_EReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter, feeding_level, ...)
{
  return(encounter - params@metab)
}

newFeedingLevel <- function (params, n, n_pp, n_other, t, encounter, ...)
{
  return(encounter * 0) #zero feeding level corresponds to type 1 feeding
}

fZooMizer_run <- function(groups, input, no_w = 178){

  kappa = 10^(input$phyto_int)
  lambda = 1-input$phyto_slope
  chlo = input$chlo
  sst = input$sst
  dt = input$dt
  tmax = input$tmax
  
    #data
  groups$w_min <- 10^groups$w_min #convert from log10 values
  groups$w_inf <- 10^groups$w_inf
  groups$w_mat <- 10^groups$w_mat
  groups$h <- 1e50 # should be Inf, but that breaks the calculations. Massive value still works out to effectively unlimited feeding as allowed in ZooMSS
  groups$ks <- 0 #turn off standard metabolism
  if (is.null(groups$interaction_resource)) {
    groups$interaction_resource <- 1
    groups$interaction_resource[which(groups$FeedType == "Carnivore")] <- 0
  }
  
  #todo - ramp up constant repro for coexistence

  mf.params <- new_newMultispeciesParams(species_params=groups,
                                         interaction=NULL, #NULL sets all to 1, no strict herbivores
                                         min_w = 10^(-10.7),
                                         max_w = 10^7* (1 + 1e-06),
                                         no_w = no_w, #number of zoo+fish size classes;
                                         # w_full = 10^seq(from = -14.5, to = (log10(max(groups$w_inf)) + 0.1), by = 0.1),
                                         min_w_pp = 10^(-14.5), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
                                         w_pp_cutoff = 10^(input$phyto_max)* (1 + 1e-06), #maximum phyto size
                                         n = 0.7, #The allometric growth exponent used in ZooMSS
                                         z0pre = 1, #external mortality (senescence)
                                         z0exp = 0.3,
                                         resource_dynamics = "phyto_fixed",
                                         kappa = kappa,
                                         lambda = lambda,
                                         RDD = constantRDD(species_params = groups), #first go at this
                                         #pred_kernel = ... #probably easiest to just import this/pre-calculate it, once dimensions are worked out
  )

  mf.params@species_params$w_min <- groups$w_min  #fix Mizer setting the egg weight to be one size larger for some groups.
  #mf.params@initial_n[] <- readRDS("data/initialn.RDS")

  temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))

  mf.params@other_params$assim_eff <- setassim_eff(groups)
  cc_phyto <- 0.1 #carbon content of phytoplankton
  mf.params@other_params$assim_phyto <- mf.params@species_params$GrossGEscale * cc_phyto #assimilation efficiency when eating phytoplankton

  mf.params@other_params$temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))

  mf.params <- setZooMizerConstants(params = mf.params, Groups = groups, sst = input$sst)
  #mf.params@initial_n[] <- readRDS("data/initialn.RDS")

  #mf.params <- setParams(mf.params)

  # mf.params <- setRateFunction(mf.params, "PredRate", "new_PredRate")
  mf.params <- setRateFunction(mf.params, "EReproAndGrowth", "new_EReproAndGrowth")
  mf.params <- setRateFunction(mf.params, "FeedingLevel", "newFeedingLevel")
  mf.params <- setRateFunction(mf.params, "Encounter", "new_Encounter")
  mf.params <- setRateFunction(mf.params, "PredRate", "new_PredRate")
  mf.params <- setReproduction(mf.params, repro_prop = matrix(0, nrow = nrow(mf.params@psi), ncol = ncol(mf.params@psi)))

  #mf.params@search_vol[] <- readRDS("data/SearchVol.rds")
  #mf.params <- setSearchVolume(mf.params)
  #mf.params <- setmort_test(mf.params, sst)
  # M_sb <- getExtMort(mf.params)
  # M_sb[] <- readRDS("data/mu_b.RDS")
  # temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  # M_sb <- temp_eff * M_sb *10 # Incorporate temp effect on senscence mortality

  # mf.params <- setExtMort(mf.params, z0 = M_sb)

  sim <- project(mf.params, dt = dt, t_max = tmax, t_save = 1) #TODO: make t_save an input to fZooMizer_run

  return(sim)
}


### try to fix issues with weight vectors disagreeing with ZooMSS - allow for w_full to be specified.
# This involves slight editing to newMultispeciesParams and un-commenting lines in emptyParams functions

# new_newMultispeciesParams <- function(
#   species_params,
#   interaction = NULL,
#   no_w = 100,
#   #w_full = NA,         #added this line
#   min_w = 0.001,
#   max_w = NA,
#   min_w_pp = NA,
#   # setPredKernel()
#   pred_kernel = NULL,
#   # setSearchVolume()
#   search_vol = NULL,
#   # setMaxIntakeRate()
#   intake_max = NULL,
#   # setMetabolicRate()
#   metab = NULL,
#   p = 0.7,
#   # setExtMort
#   z0 = NULL,
#   z0pre = 0.6,
#   z0exp = n - 1,
#   # setReproduction
#   maturity = NULL,
#   repro_prop = NULL,
#   RDD = "BevertonHoltRDD",
#   # setResource
#   resource_rate = NULL,
#   resource_capacity = NULL,
#   n = 2 / 3,
#   r_pp = 10,
#   kappa = 1e11,
#   lambda = 2.05,
#   w_pp_cutoff = 10,
#   resource_dynamics = "resource_semichemostat",
#   # setFishing
#   gear_params = data.frame(),
#   selectivity = NULL,
#   catchability = NULL,
#   initial_effort = NULL) {
#   no_sp <- nrow(species_params)
# 
#   ## For backwards compatibility, allow r_max instead of R_max
#   if (!("R_max" %in% names(species_params)) &&
#       "r_max" %in% names(species_params)) {
#     names(species_params)[names(species_params) == "r_max"] <- "R_max"
#   }
# 
#   ## Create MizerParams object ----
#   params <- new_emptyParams(species_params,
#                             gear_params,
#                             # w_full = w_full,
#                             no_w = no_w,
#                             min_w = min_w,
#                             max_w = max_w,
#                             min_w_pp = min_w_pp)
# 
#   ## Fill the slots ----
#   params <- params %>%
#     set_species_param_default("n", n) %>%
#     set_species_param_default("p", p)
#   params <- set_species_param_default(params, "q",
#                                       lambda - 2 + params@species_params$n)
#   if (is.null(interaction)) {
#     interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
#   }
#   params <-
#     setParams(params,
#               # setInteraction
#               interaction = interaction,
#               # setPredKernel()
#               pred_kernel = pred_kernel,
#               # setSearchVolume()
#               search_vol = search_vol,
#               # setMaxIntakeRate()
#               intake_max = intake_max,
#               # setMetabolicRate()
#               metab = metab,
#               # setExtMort
#               z0 = z0,
#               z0pre = z0pre,
#               z0exp = z0exp,
#               # setReproduction
#               maturity = maturity,
#               repro_prop = repro_prop,
#               RDD = RDD,
#               # setResource
#               resource_rate = resource_rate,
#               resource_capacity = resource_capacity,
#               r_pp = r_pp,
#               kappa = kappa,
#               lambda = lambda,
#               n = n,
#               w_pp_cutoff = w_pp_cutoff,
#               resource_dynamics = resource_dynamics,
#               # setFishing
#               gear_params = gear_params,
#               selectivity = selectivity,
#               catchability = catchability,
#               initial_effort = initial_effort)
# 
#   params@initial_n <- get_initial_n(params)
#   params@initial_n_pp <- params@cc_pp
#   params@A <- rep(1, nrow(species_params))
# 
#   return(params)
# }
# 
# new_emptyParams <- function(species_params,
#                             gear_params = data.frame(),
#                             no_w = 100,
#                             min_w = 0.001,
#                             w_full = NA,
#                             max_w = NA,
#                             min_w_pp = 1e-12) {
#   assert_that(is.data.frame(species_params),
#               is.data.frame(gear_params),
#               no_w > 10)
# 
#   ## Set defaults ----
#   if (is.na(min_w_pp)) min_w_pp <- 1e-12
#   species_params <- set_species_param_default(species_params, "w_min", min_w)
#   min_w <- min(species_params$w_min)
# 
#   species_params <- validSpeciesParams(species_params)
#   gear_params <- validGearParams(gear_params, species_params)
# 
#   if (is.na(max_w)) {
#     max_w <- max(species_params$w_inf)
#   } else {
#     if (max(species_params$w_inf) > max_w * (1 + 1e-6)) { # The fudge factor
#       # is there to avoid false alerts due to rounding errors.
#       too_large <- species_params$species[max_w < species_params$w_inf]
#       stop("Some of your species have an maximum size larger than max_w: ",
#            toString(too_large))
#     }
#   }
# 
#   # Set up grids ----
#   if (missing(w_full)) {
#     # set up logarithmic grids
#     dx <- log10(max_w / min_w) / (no_w - 1)
#     # Community grid
#     w <- 10^(seq(from = log10(min_w), by = dx, length.out = no_w))
#     # dw[i] = w[i+1] - w[i]. Following formula works also for last entry dw[no_w]
#     dw <- (10^dx - 1) * w
#     # To avoid issues due to numerical imprecision
#     min_w <- w[1]
# 
#     # For fft methods we need a constant log bin size throughout.
#     # Therefore we use as many steps as are necessary so that the first size
#     # class includes min_w_pp.
#     x_pp <- rev(seq(from = log10(min_w),
#                     to = log10(min_w_pp),
#                     by = -dx)) - dx
#     w_full <- c(10^x_pp, w)
#     # If min_w_pp happened to lie exactly on a grid point, we now added
#     # one grid point too much which we need to remove again
#     if (w_full[2] == min_w_pp) {
#       w_full <- w_full[2:length(w_full)]
#     }
#     no_w_full <- length(w_full)
#     dw_full <- (10^dx - 1) * w_full
#   } else {
#     #     # use supplied w_full
#     no_w_full <- length(w_full) - 1
#     dw_full <- diff(w_full)
#     w_full <- w_full[seq_along(dw_full)]
#     #     # Check that sizes are increasing
#     if (any(dw_full <= 0)) {
#       stop("w_full must be increasing.")
#     }
#     w_min_idx <- match(min_w, w_full)
#     if (is.na(w_min_idx)) {
#       stop("w_min must be contained in w_full.")
#     }
#     w <- w_full[w_min_idx:no_w_full]
#     dw <- dw_full[w_min_idx:no_w_full]
#     no_w <- length(w)
#     min_w_pp <- w_full[1]
#   }
# 
#   # Basic arrays for templates ----
#   no_sp <- nrow(species_params)
#   species_names <- as.character(species_params$species)
#   gear_names <- unique(gear_params$gear)
#   mat1 <- array(0, dim = c(no_sp, no_w),
#                 dimnames = list(sp = species_names, w = signif(w,3)))
#   ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
#                           dimnames = list(sp = species_names, k = 1:no_w_full))
#   ft_mask <- plyr::aaply(species_params$w_inf, 1,
#                          function(x) w_full < x, .drop = FALSE)
# 
#   selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w),
#                        dimnames = list(gear = gear_names, sp = species_names,
#                                        w = signif(w, 3)))
#   catchability <- array(0, dim = c(length(gear_names), no_sp),
#                         dimnames = list(gear = gear_names, sp = species_names))
#   initial_effort <- rep(0, length(gear_names))
#   names(initial_effort) <- gear_names
# 
#   interaction <- array(1, dim = c(no_sp, no_sp),
#                        dimnames = list(predator = species_names,
#                                        prey = species_names))
# 
#   vec1 <- as.numeric(rep(NA, no_w_full))
#   names(vec1) <- signif(w_full, 3)
# 
#   # Round down w_min to lie on grid points and store the indices of these
#   # grid points in w_min_idx
#   w_min_idx <- as.vector(suppressWarnings(
#     tapply(species_params$w_min, 1:no_sp,
#             function(w_min, wx) max(which(wx <= w_min)), wx = w)))
#   #         function(w_min, wx) max(which(round(log10(wx), digits = 2) <= round(log10(w_min), digits = 2))), wx = w)))
#   # Due to rounding errors this might happen:
#   w_min_idx[w_min_idx == -Inf] <- 1
#   names(w_min_idx) <- species_names
#   species_params$w_min <- w[w_min_idx]
#   
#   # if (any(species_params$w_min < w[w_min_idx]) || any(species_params$w_min > w[w_min_idx + 1])) {
#   #   msg <- "The `w_min_idx` should point to the start of the size bin containing the egg size `w_min`."
#   #   errors <- c(errors, msg)
#   # }
# 
#   # Colour and linetype scales ----
#   # for use in plots
#   # Colour-blind-friendly palettes
#   # From http://dr-k-lo.blogspot.co.uk/2013/07/a-color-blind-friendly-palette-for-r.html
#   # cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00",
#   #                 "#CC79A7", "#F0E442")
#   # From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
#   # cbbPalette <- c("#E69F00", "#56B4E9", "#009E73",
#   #                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   # Random palette gemerated pm https://medialab.github.io/iwanthue/
#   colour_palette <- c("#815f00",
#                       "#6237e2",
#                       "#8da600",
#                       "#de53ff",
#                       "#0e4300",
#                       "#430079",
#                       "#6caa72",
#                       "#ee0053",
#                       "#007957",
#                       "#b42979",
#                       "#142300",
#                       "#a08dfb",
#                       "#644500",
#                       "#04004c",
#                       "#b79955",
#                       "#0060a8",
#                       "#dc8852",
#                       "#007ca9",
#                       "#ab003c",
#                       "#9796d9",
#                       "#472c00",
#                       "#b492b0",
#                       "#140000",
#                       "#dc8488",
#                       "#005c67",
#                       "#5c585a")
#   # type_palette <- c("solid", "dashed", "dotdash", "longdash",
#   #                   "twodash")
#   type_palette <- c("solid")
# 
#   if ("linecolour" %in% names(species_params)) {
#     linecolour <- species_params$linecolour
#     # If any NA's first fill them with unused colours
#     linecolour[is.na(linecolour)] <-
#       setdiff(colour_palette, linecolour)[1:sum(is.na(linecolour))]
#     # if there are still NAs, start from beginning of palette again
#     linecolour[is.na(linecolour)] <-
#       colour_palette[1:sum(is.na(linecolour))]
#   } else {
#     linecolour <- rep(colour_palette, length.out = no_sp)
#   }
#   names(linecolour) <- as.character(species_names)
#   linecolour <- c(linecolour, "Total" = "black", "Resource" = "green",
#                   "Background" = "grey", "Fishing" = "red")
# 
#   if ("linetype" %in% names(species_params)) {
#     linetype <- species_params$linetype
#     linetype[is.na(linetype)] <- "solid"
#   } else {
#     linetype <- rep(type_palette, length.out = no_sp)
#   }
#   names(linetype) <- as.character(species_names)
#   linetype <- c(linetype, "Total" = "solid", "Resource" = "solid",
#                 "Background" = "solid", "Fishing" = "solid")
# 
#   # Make object ----
#   # Should Z0, rrPP and ccPP have names (species names etc)?
#   params <- new(
#     "MizerParams",
#     w = w,
#     dw = dw,
#     w_full = w_full,
#     dw_full = dw_full,
#     w_min_idx = w_min_idx,
#     maturity = mat1,
#     psi = mat1,
#     initial_n = mat1,
#     intake_max = mat1,
#     search_vol = mat1,
#     metab = mat1,
#     mu_b = mat1,
#     ft_pred_kernel_e = ft_pred_kernel,
#     ft_pred_kernel_p = ft_pred_kernel,
#     pred_kernel = array(),
#     gear_params = gear_params,
#     selectivity = selectivity,
#     catchability = catchability,
#     initial_effort = initial_effort,
#     rr_pp = vec1,
#     cc_pp = vec1,
#     sc = w,
#     initial_n_pp = vec1,
#     species_params = species_params,
#     interaction = interaction,
#     other_dynamics = list(),
#     other_encounter = list(),
#     other_mort = list(),
#     rates_funcs = list(
#       Rates = "mizerRates",
#       Encounter = "mizerEncounter",
#       FeedingLevel = "mizerFeedingLevel",
#       EReproAndGrowth = "mizerEReproAndGrowth",
#       PredRate = "mizerPredRate",
#       PredMort = "mizerPredMort",
#       FMort = "mizerFMort",
#       Mort = "mizerMort",
#       ERepro = "mizerERepro",
#       EGrowth = "mizerEGrowth",
#       ResourceMort = "mizerResourceMort",
#       RDI = "mizerRDI",
#       RDD = "BevertonHoltRDD"),
#     resource_dynamics = "resource_semichemostat",
#     other_params = list(),
#     initial_n_other = list(),
#     A = as.numeric(rep(NA, no_sp)),
#     linecolour = linecolour,
#     linetype = linetype,
#     ft_mask = ft_mask
#   )
# 
#   return(params)
# }

# setmort_test <- function(params, sst){
#   #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
#   M_sb <- getExtMort(params)
#   M_sb[] <- readRDS("data/mu_b.RDS")
#   ZSpre <- 1 # senescence mortality prefactor
#   ZSexp <- 0.3 # senescence mortality exponent
#   #
#   # for(i in 1:nrow(params@species_params)){
#   #   ## Senescence mortality
#   #   if(params@species_params$Type[i] == "Zooplankton"){
#   #     M_sb[i,] <- ZSpre*(params@w/(params@species_params$w_mat[i]))^ZSexp
#   #     M_sb[i, params@species_params$w_inf[i] < params@w] <- 0
#   #     M_sb[i, params@species_params$w_mat[i] > params@w] <- 0
#   #   }
#   #
#   #   if(params@species_params$Type[i] == "Fish"){
#   #     M_sb[i,] <- 0.1*ZSpre*(params@w/(params@species_params$w_mat[i]))^ZSexp
#   #     M_sb[i, params@species_params$w_inf[i] < params@w] <- 0
#   #     M_sb[i, params@species_params$w_mat[i] > params@w] <- 0
#   #   }
#   #   }
# 
#   #temperature effect
# 
#   M_sb <- params@other_params$temp_eff * M_sb # Incorporate temp effect on senscence mortality
# 
#   params <- setExtMort(params, z0 = M_sb)
# 
#   return(params)
# }