#' Get Biomass of zooplankton and fish groups
#'
#' @param sim An object of class MizerSim containing the fish model
#' @param zoo_params (Optional) An object of class MizerParams containing the zooplankton model parameters
#' @param ... 
#'
#' @return An array (time x species) containing the biomass in grams.
#' @export
#'
#' @examples
getBiomass_ZooMizer <- function (sim, zoo_params = sim@params@other_params$zoo$params, ...) 
{
  assert_that(is(sim, "MizerSim"))
  assert_that(is(zoo_params, "MizerParams"))
  zoo_size_range <- get_size_range_array(zoo_params, ...)
  zoo_n <- aperm(array(unlist(sim@n_other), c(nrow(zoo_params@species_params),length(zoo_params@w),length(sim@n_other))), c(3,1,2))
  dimnames(zoo_n) <- list(getTimes(sim), zoo_params@species_params$species, zoo_params@w)
  zoo_biomass <- apply(sweep(sweep(zoo_n, c(2, 3), zoo_size_range, 
                                   "*"), 3, zoo_params@w * zoo_params@dw, "*"), 
                       c(1, 2), sum)
  fish_size_range <- get_size_range_array(sim@params, ...)
  fish_biomass <- apply(sweep(sweep(sim@n, c(2, 3), fish_size_range, 
                               "*"), 3, sim@params@w * sim@params@dw, "*"), 
                   c(1, 2), sum)
  return(cbind(zoo_biomass, fish_biomass))
}

#' Get abundance of zooplankton and fish groups
#'
#' @param sim An object of class MizerSim containing the fish model
#' @param zoo_params (Optional) An object of class MizerParams containing the zooplankton model parameters
#' @param ... 
#'
#' @return An array (time x species) containing the abundances.
#' @export
#'
#' @examples
getAbundance_ZooMizer <- function (sim, zoo_params = sim@params@other_params$zoo$params, ...) 
{
  assert_that(is(sim, "MizerSim"))
  assert_that(is(zoo_params, "MizerParams"))
  zoo_size_range <- get_size_range_array(zoo_params, ...)
  zoo_n <- aperm(array(unlist(sim@n_other), c(nrow(zoo_params@species_params),length(zoo_params@w),length(sim@n_other))), c(3,1,2))
  dimnames(zoo_n) <- list(getTimes(sim), zoo_params@species_params$species, zoo_params@w)
  zoo_abunds <- apply(sweep(sweep(zoo_n, c(2, 3), zoo_size_range, 
                                   "*"), 3, zoo_params@dw, "*"), 
                       c(1, 2), sum)
  fish_size_range <- get_size_range_array(sim@params, ...)
  fish_abunds <- apply(sweep(sweep(sim@n, c(2, 3), fish_size_range, 
                                    "*"), 3, sim@params@dw, "*"), 
                        c(1, 2), sum)
  return(cbind(zoo_abunds, fish_abunds))
}

#' Get diet per predator/predator size/prey/prey size
#' Original code by Romain Forestier
#'
#' @param params An object of class MizerParams
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list with the abundances of other components
#' @param proportion Boolean whether diet composition should be returned as a proportion.
#'     Defaults to true.
#'
#' @return A 4-dimensional array (predator species x predator weight x prey x prey weight)
#' @export
#'
#' @examples
getDietComp <- function (params, n = initialN(params), n_pp = initialNResource(params),
                         n_other = initialNOther(params), proportion = TRUE)
{
  params <- validParams(params)
  pred_kernel(params) <- pred_kernel(params) #required to make this work
  species <- params@species_params$species
  no_sp <- length(species)
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  no_other <- length(params@other_encounter)
  other_names <- names(params@other_encounter)
  assertthat::assert_that(identical(dim(n), c(no_sp, no_w)), length(n_pp) ==
                            no_w_full)
  
  n_tot <- sweep(n, 2, params@w * params@dw, "*")
  diet_comp <- array(0, dim = c(dim(n_tot),no_sp + 1, no_w_full),
                     dimnames = list("Predator" = species,
                                     "wPredator" = params@w,
                                     "Prey" = c(as.character(species), "Resource"),#, other_names),
                                     "wPrey" = params@w_full))
  
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  # if (!is.null(comment(params@pred_kernel))) {
    ## code to keep prey size | works only without other backgrounds now
    
    for(iPredator in 1:dim(diet_comp)[1])
    {
      for(wPredator in 1:dim(diet_comp)[2])
      {
        for(iPrey in 1:no_sp)
        {
          for(wPrey in 1:no_w) # assuming that w is the tail of w_full, won't work if w_full gets larger than w
          {
            diet_comp[iPredator,wPredator,iPrey,wPrey] <- params@pred_kernel[iPredator,wPredator, (idx_sp[wPrey])] *
              n_tot[iPrey,wPrey]
            
          }
        }
      }
    }
    diet_comp[,,no_sp+1,] <- sweep(params@pred_kernel,
                                   3, params@dw_full * params@w_full * n_pp, "*")
    
    ###
    
  # }
  # else {
  #   prey <- matrix(0, nrow = no_sp + 1, ncol = no_w_full)
  #   prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw,
  #                                  "*")
  #   prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
  #   ft <- array(rep(params@ft_pred_kernel_e, times = no_sp +
  #                     1) * rep(mvfft(t(prey)), each = no_sp), dim = c(no_sp,
  #                                                                     no_w_full, no_sp + 1))
  #   ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
  #   ae <- array(Re(mvfft(ft, inverse = TRUE)/no_w_full),
  #               dim = c(no_w_full, no_sp, no_sp + 1))
  #   ae <- ae[idx_sp, , , drop = FALSE]
  #   ae <- aperm(ae, c(2, 1, 3))
  #   ae[ae < 1e-18] <- 0
  #   diet_comp[, , 1:(no_sp + 1),] <- ae
  # }
  
  inter <- cbind(params@interaction, params@species_params$interaction_resource)
  
  diet_comp[, , 1:(no_sp + 1),] <- sweep(sweep(diet_comp[, , 1:(no_sp + 1),, drop = FALSE], c(1, 3), inter, "*"), c(1, 2), params@search_vol, "*")
  
  for (i in seq_along(params@other_encounter)) {
    diet_comp[, , no_sp + 1 + i,] <- do.call(params@other_encounter[[i]],
                                             list(params = params, n = n, n_pp = n_pp, n_other = n_other,
                                                  component = names(params@other_encounter)[[i]]))
  }
  f <- getFeedingLevel(params, n, n_pp)
  fish_mask <- n > 0
  diet_comp <- sweep(diet_comp, c(1, 2), (1 - f) * fish_mask, "*")
  if (proportion) {
    total <- rowSums(diet_comp, dims = 2)
    diet_comp <- sweep(diet_comp, c(1, 2), total, "/")
    diet_comp[is.nan(diet_comp)] <- 0
  }
  return(diet_comp)
}

#' Extract realised PPMR for each fish species weighted by predator/prey densities
#' Original code by Romain Forestier 
#'
#' @param params An object of class `MizerParams`
#'
#' @return A named vector of realised PPMRs for each species.
#' @export
#'
#' @examples
#' \donttest{
#' getRealisedPPMRs(NS_params)
#' }
getRealisedPPMRs <- function(params) {
  SpIdx <- 1:nrow(params@species_params)
  rPPMR <- vector("numeric", length = length(SpIdx)) # SpIdx is the species index (numeric in my model)
  
  DietComp <- getDietComp(params)
  
  for(iSpecies in SpIdx) {
    speciesDat <- DietComp[iSpecies,,,] # this is diet per predator/predator size/prey/prey size  
    # no need to know prey identity
    
    # speciesDat <- apply(speciesDat,c(1,3),sum)
    speciesDat <- colSums(aperm(speciesDat, c(2,1,3))) # slightly faster
    
    # for each size class need PPMR value
    speciesPPMR <- NULL
    for(iSize in 1:length(dimnames(speciesDat)$wPred)) {
      
      if(sum(speciesDat[iSize,])) {# if there is something at this size
        sizeDat <- speciesDat[iSize,]
        # for now the easy way
        # PreferredSizeClass <- which(sizeDat == max(sizeDat))  # in my model all diet profiles are bell curves on the prey size       
        PreferredSizeClass <- weighted.mean(params@w_full,sizeDat, na.rm = TRUE)
        sizePPMR <- as.numeric(dimnames(speciesDat)$wPred[iSize])/PreferredSizeClass
        speciesPPMR <- c(speciesPPMR,sizePPMR)
      }
    }
    rPPMR[iSpecies] <- weighted.mean(speciesPPMR,params@initial_n[iSpecies,which(params@initial_n[iSpecies,]>0)],na.rm = T)
  }
  names(rPPMR) <- params@species_params$species
  return(rPPMR)
}
