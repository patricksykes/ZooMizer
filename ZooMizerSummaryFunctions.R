#' Get Biomass of zooplankton and fish groups
#'
#' @param sim An object of class MizerSim containing the fish model
#' @param zoo_params An object of class MizerParams containing the zooplankton model parameters
#' @param ... 
#'
#' @return An array (time x species) containing the biomass in grams.
#' @export
#'
#' @examples
getBiomass_ZooMizer <- function (sim, zoo_params, ...) 
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
#' @param zoo_params An object of class MizerParams containing the zooplankton model parameters
#' @param ... 
#'
#' @return An array (time x species) containing the abundances.
#' @export
#'
#' @examples
getAbundance_ZooMizer <- function (sim, zoo_params, ...) 
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
