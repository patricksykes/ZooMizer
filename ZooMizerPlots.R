#' Title
#'
#' @param fish_object An object of class \linkS4class{MizerSim} or 
#'   \linkS4class{MizerParams}.
#' @param zoo_object (Optional) An object of class MizerParams that contains the
#'   zooplankton model parameters
#' @inheritParams valid_species_arg
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step. Ignored when called with a \linkS4class{MizerParams}
#'   object.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param power The abundance is plotted as the number density times the weight
#' raised to `power`. The default \code{power = 1} gives the biomass 
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param biomass `r lifecycle::badge("deprecated")`
#'  Only used if `power` argument is missing. Then
#'   \code{biomass = TRUE} is equivalent to \code{power=1} and 
#'   \code{biomass = FALSE} is equivalent to \code{power=0}
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param resource A boolean value that determines whether resource is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#'
#' @return
#' @export
#'
#' @examples
plotSpectra_ZooMizer <- function(fish_object, zoo_object,
                        species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE, 
                        background = TRUE,
                        highlight = NULL, ...) {
  # to deal with old-type biomass argument
  if (missing(power)) {
    power <- as.numeric(biomass)
  }
  species <- valid_species_arg(fish_object, species)
  if (is(fish_object, "MizerSim")) {
    if (missing(time_range)) {
      time_range  <- max(as.numeric(dimnames(fish_object@n)$time))
    }
    if (missing(zoo_object)) {
      zoo_object <- fish_object@params@other_params$zoo$params
    }
    time_elements <- get_time_elements(fish_object, time_range)
    n <- apply(fish_object@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    n_pp <- zoo_object@initial_n_pp
    n_zoo <- aperm(simplify2array(fish_object@n_other), c(3,1,2))
    n_zoo <- apply(n_zoo[time_elements, , , drop = FALSE], c(2, 3), mean)
    ps <- plot_spectra_zoomizer(fish_object@params, zoo_object,
                       n = n, n_pp = n_pp, n_zoo = n_zoo,
                       species = species, wlim = wlim,
                       ylim = ylim, power = power,
                       total = total, resource = resource,
                       background = background, highlight = highlight)
    return(ps)
  } else {
    if (missing(zoo_object)) {
      zoo_object <- fish_object@other_params$zoo$params
    }
    ps <- plot_spectra_zoomizer(fish_object, zoo_params = fish_object@other_params$zoo$params,
                       n = fish_object@initial_n,
                       n_zoo = zoo_object@initial_n,
                       n_pp = zoo_object@initial_n_pp,
                       species = species, wlim = wlim, ylim = ylim,
                       power = power,
                       total = total, resource = resource,
                       background = background, highlight = highlight)
    return(ps)
  }
}


plot_spectra_zoomizer <- function(fish_params, zoo_params, n, n_pp, n_zoo,
                         species, wlim, ylim, power,
                         total, resource, background, highlight) {
  fish_params <- validParams(fish_params)
  if (is.na(wlim[1])) {
    wlim[1] <- min(fish_params@w) / 100
  }
  if (is.na(wlim[2])) {
    wlim[2] <- max(fish_params@w_full)
  }
  # Need to keep species in order for legend
  species_levels <- c(dimnames(fish_params@initial_n)$sp,
                      "Background", dimnames(zoo_params@initial_n)$sp, "Total")
  fish_idx <- (length(fish_params@w_full) - length(fish_params@w) + 1):length(fish_params@w_full)
  zoo_idx <- (length(zoo_params@w_full) - length(zoo_params@w) + 1):length(zoo_params@w_full)
  if (total) {
    # Calculate total community abundance
    total_n <- n_pp
    total_n[zoo_idx] <- total_n[zoo_idx] + colSums(n_zoo)
    total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
    total_n <- total_n * zoo_params@w_full^power
  }
  species <- c(valid_species_arg(zoo_params), valid_species_arg(fish_params, species))
  # Deal with power argument
  if (power %in% c(0, 1, 2)) {
    y_label <- c("Number density [1/g]", "Biomass density",
                 "Biomass density [g]")[power + 1]
  } else {
    y_label <- paste0("Number density * w^", power)
  }
  n <- sweep(n, 2, fish_params@w^power, "*")
  n_zoo <- sweep(n_zoo, 2, zoo_params@w^power, "*")
  n_mat <- matrix(0, nrow = length(species), ncol = length(zoo_params@w_full))
  n_mat[1:nrow(n_zoo), zoo_idx] <- n_zoo
  n_mat[(nrow(n_zoo)+1):length(species), fish_idx] <- n
  rownames(n_mat) <- species
  colnames(n_mat) <- zoo_params@w_full
  # Select only the desired species
  spec_n <- n_mat[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  
  # Make data.frame for plot
  plot_dat <- data.frame(value = c(spec_n),
                         # ordering of factor is important for legend
                         Species = factor(dimnames(spec_n)[[1]],
                                          levels = species_levels),
                         w = rep(zoo_params@w_full,
                                 each = dim(spec_n)[[1]]))
  
  if (resource) {
    resource_sel <- (zoo_params@w_full >= wlim[1]) & 
      (zoo_params@w_full <= wlim[2])
    # Do we have any resource to plot?
    if (sum(resource_sel) > 0) {
      w_resource <- zoo_params@w_full[resource_sel]
      plank_n <- n_pp[resource_sel] * w_resource^power
      plot_dat <- rbind(plot_dat,
                        data.frame(value = c(plank_n),
                                   Species = "Resource",
                                   w = w_resource))
    }
  }
  if (total) {
    plot_dat <- rbind(plot_dat,
                      data.frame(value = c(total_n),
                                 Species = "Total",
                                 w = zoo_params@w_full))
  }
  # lop off 0s and apply wlim
  plot_dat <- plot_dat[(plot_dat$value > 0) & 
                         (plot_dat$w >= wlim[1]) &
                         (plot_dat$w <= wlim[2]), ]
  # Impose ylim
  if (!is.na(ylim[2])) {
    plot_dat <- plot_dat[plot_dat$value <= ylim[2], ]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- 1e-20
  }
  plot_dat <- plot_dat[plot_dat$value > ylim[1], ]
  
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, fish_params@linecolour)
  
  # Create plot
  

  p <- ggplot(plot_dat, aes(x = w, y = value)) +
    scale_x_continuous(name = "Size [g]", trans = "log10",
                       breaks = log_breaks()) +
    scale_y_continuous(name = y_label, trans = "log10",
                       breaks = log_breaks()) +
    scale_colour_manual(values = linecolours) +
    scale_linetype_manual(values = c(zoo_params@linetype[1:nrow(zoo_params@species_params)], fish_params@linetype))
  if (background) {
    back_n <- n[is.na(fish_params@A), , drop = FALSE]
    plot_back <- data.frame(value = c(back_n),
                            Species = as.factor(dimnames(back_n)[[1]]),
                            w = rep(zoo_params@w_full,
                                    each = dim(back_n)[[1]]))
    # lop off 0s and apply wlim
    plot_back <- plot_back[(plot_back$value > 0) & 
                             (plot_back$w >= wlim[1]) &
                             (plot_back$w <= wlim[2]), ]
    # Impose ylim
    if (!is.na(ylim[2])) {
      plot_back <- plot_back[plot_back$value <= ylim[2], ]
    }
    plot_back <- plot_back[plot_back$value > ylim[1], ]
    if (nrow(plot_back) > 0) {
      # Add background species
      p <- p +
        geom_line(aes(group = Species),
                  colour = fish_params@linecolour["Background"],
                  linetype = fish_params@linetype["Background"],
                  data = plot_back)
    }
  }
  linesize <- rep(0.8, length(c(zoo_params@linetype, fish_params@linetype))-4)
  names(linesize) <- c(zoo_params@species_params$species, names(fish_params@linetype))
  linesize[highlight] <- 1.6
  p <- p + scale_size_manual(values = linesize) + 
    geom_line(aes(colour = Species, linetype = Species, size = Species))
  return(p)
}

#' @rdname plotSpectra_ZooMizer
#' @export
plotlySpectra_ZooMizer <- function(fish_object, zoo_object, species = NULL,
                          time_range,
                          wlim = c(NA, NA), ylim = c(NA, NA),
                          power = 1, biomass = TRUE,
                          total = FALSE, resource = TRUE, 
                          background = TRUE,
                          highlight = NULL, ...) {
  argg <- as.list(environment())
  ggplotly(do.call("plotSpectra_ZooMizer", argg))
}


getBiomassFrame_ZooMizer <- function (sim, zoo_params, species = NULL, start_time = as.numeric(dimnames(sim@n)[[1]][1]), 
          end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]), 
          ylim = c(NA, NA), total = FALSE, ...) 
{
  species <- valid_species_arg(sim@params, species)
  b <- getBiomass_ZooMizer(sim, zoo_params, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= 
                                                           end_time), , drop = FALSE]
  b_total <- rowSums(b)
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >= 
                                      ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[2]), ]
  names(bm) <- c("Year", "Species", "Biomass")
  species_levels <- c(zoo_params@species_params$species, dimnames(sim@n)$sp, "Background", 
                      "Resource", "Total")
  bm$Species <- factor(bm$Species, levels = species_levels)
  # bm <- bm[bm$Species %in% species, ]
  return(bm)
}

#' Plot the biomass of species through time
#'
#' @param sim An object of class MizerSim contianing the fish model
#' @param zoo_params An object of class MizerParams containing the zooplankton model
#' @param species A list of fish species to include in the plot
#' @param start_time 
#' @param end_time 
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param ylim 
#' @param total 
#' @param background 
#' @param highlight 
#' @param ... 
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
plotBiomass_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
          ylim = c(NA, NA), total = FALSE, background = FALSE, highlight = NULL, 
          ...) 
{
  species <- valid_species_arg(sim, species)
  if (missing(zoo_params)) {
    zoo_params <- sim@params@other_params$zoo$params
  }
  if (missing(start_time)) 
    start_time <- as.numeric(dimnames(sim@n)[[1]][1])
  if (missing(end_time)) 
    end_time <- as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
  bm <- getBiomassFrame_ZooMizer(sim, zoo_params, species = dimnames(sim@n)$sp, 
                        start_time = start_time, end_time = end_time, ylim = ylim, 
                        total = total)#, ...)
  spec_bm <- bm[bm$Species %in% c("Total", zoo_params@species_params$species, species), ]
  x_label <- "Year"
  y_label <- "Biomass [g]"
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, sim@params@linecolour)
  p <- ggplot(spec_bm, aes(x = Year, y = Biomass)) +
    scale_y_continuous(trans = "log10", 
                       breaks = mizer:::log_breaks(n = y_ticks),
                       labels = prettyNum, 
                       name = y_label) +
    scale_x_continuous(name = x_label) + 
    scale_colour_manual(values = linecolours)
    # scale_linetype_manual(values = sim@params@linetype)
  if (background) {
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    if (nrow(back_bm) > 0) {
      p <- p + geom_line(aes(group = Species), data = back_bm, 
                         colour = sim@params@linecolour["Background"], 
                         linetype = sim@params@linetype["Background"])
    }
  }
  linesize <- rep(0.8, length(c(zoo_params@species_params$species,sim@params@linetype)))
  names(linesize) <- c(zoo_params@species_params$species, names(sim@params@linetype))
  linesize[highlight] <- 1.6
  p <- p + geom_line(aes(colour = Species ,linetype = Species , size = Species)) +
         scale_linetype_manual(values = c(zoo_params@linetype[1:nrow(zoo_params@species_params)], sim@params@linetype))+
         scale_size_manual(values = linesize)
  return(p)
}


#' @rdname plotBiomass_ZooMizer
#' @export
plotlyBiomass_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
          ylim = c(NA, NA), total = FALSE, background = TRUE, highlight = NULL, 
          ...) 
{
  argg <- c(as.list(environment()), list(...))
  ggplotly(do.call("plotBiomass_ZooMizer", argg))
}

getAbundanceFrame_ZooMizer <- function (sim, zoo_params, species = NULL, start_time = as.numeric(dimnames(sim@n)[[1]][1]), 
                                      end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]), 
                                      ylim = c(NA, NA), total = FALSE, ...) 
{
  species <- valid_species_arg(sim@params, species)
  b <- getAbundance_ZooMizer(sim, zoo_params, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= 
                                                           end_time), , drop = FALSE]
  b_total <- rowSums(b)
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >= 
                                      ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[2]), ]
  names(bm) <- c("Year", "Species", "Biomass")
  species_levels <- c(zoo_params@species_params$species, dimnames(sim@n)$sp, "Background", 
                      "Resource", "Total")
  bm$Species <- factor(bm$Species, levels = species_levels)
  # bm <- bm[bm$Species %in% species, ]
  return(bm)
}

#' Plot the abundance of species through time
#'
#' @param sim An object of class MizerSim contianing the fish model
#' @param zoo_params An object of class MizerParams containing the zooplankton model
#' @param species A list of fish species to include in the plot
#' @param start_time 
#' @param end_time 
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param ylim 
#' @param total 
#' @param background 
#' @param highlight 
#' @param ... 
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
plotAbundance_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
                                  ylim = c(NA, NA), total = FALSE, background = FALSE, highlight = NULL, 
                                  ...) 
{
  species <- valid_species_arg(sim, species)
  if (missing(start_time)) 
    start_time <- as.numeric(dimnames(sim@n)[[1]][1])
  if (missing(end_time)) 
    end_time <- as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
  bm <- getAbundanceFrame_ZooMizer(sim, zoo_params, species = dimnames(sim@n)$sp, 
                                 start_time = start_time, end_time = end_time, ylim = ylim, 
                                 total = total)#, ...)
  spec_bm <- bm[bm$Species %in% c("Total", zoo_params@species_params$species, species), ]
  x_label <- "Year"
  y_label <- "Abundance"
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, sim@params@linecolour)
  p <- ggplot(spec_bm, aes(x = Year, y = Biomass)) +
    scale_y_continuous(trans = "log10", 
                       breaks = mizer:::log_breaks(n = y_ticks),
                       labels = prettyNum, 
                       name = y_label) +
    scale_x_continuous(name = x_label) + 
    scale_colour_manual(values = linecolours)
  # scale_linetype_manual(values = sim@params@linetype)
  if (background) {
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    if (nrow(back_bm) > 0) {
      p <- p + geom_line(aes(group = Species), data = back_bm, 
                         colour = sim@params@linecolour["Background"], 
                         linetype = sim@params@linetype["Background"])
    }
  }
  linesize <- rep(0.8, length(sim@params@linetype))
  names(linesize) <- names(sim@params@linetype)
  linesize[highlight] <- 1.6
  p <- p + # scale_size_manual(values = linesize) +
    geom_line(aes(colour = Species)) #, linetype = Species)) #, size = Species))
  return(p)
}

getDiet_ZooMizer <- function(fish_params,
                             zoo_params = fish_params@other_params$zoo$params,
                             n = initialN(fish_params), 
                             n_pp = initialNResource(fish_params),
                             n_other = initialNOther(fish_params),
                             proportion = TRUE) {
  # The code is based on that for getEncounter()
  fish_params <- validParams(fish_params)
  zoo_params <- validParams(zoo_params)
  species <- fish_params@species_params$species
  zoo_species <- zoo_params@species_params$species
  no_sp <- length(species)
  no_zoo <- length(zoo_species)
  no_w <- length(fish_params@w)
  no_w_zoo <- length(zoo_params@w)
  no_w_full <- length(fish_params@w_full)
  no_other <- length(fish_params@other_encounter)
  other_names <- names(fish_params@other_encounter)
  assert_that(identical(dim(n), c(no_sp, no_w)),
              length(n_pp) == no_w_full)
  diet <- array(0, dim = c(no_sp, no_w, no_sp + no_zoo + no_other),
                dimnames = list("predator" = species,
                                "w" = dimnames(fish_params@initial_n)$w,
                                "prey" = c(as.character(species), 
                                           as.character(zoo_species),
                                           other_names)))
  # idx_sp are the index values of object@w_full such that
  # object@w_full[idx_sp] = object@w
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  idx_zoo <- (no_w_full - no_w_zoo + 1):no_w_full
    
  # If the user has set a custom kernel we can not use fft.
  if (!is.null(comment(fish_params@pred_kernel))) {
    # pred_kernel is predator species x predator size x prey size
    # We want to multiply this by the prey abundance, which is
    # prey species by prey size, sum over prey size. We use matrix
    # multiplication for this. Then we multiply 1st and 3rd 
    ae <- matrix(fish_params@pred_kernel[, , idx_sp, drop = FALSE],
                 ncol = no_w) %*%
      t(sweep(n, 2, fish_params@w * fish_params@dw, "*"))
    diet[, , 1:no_sp] <- ae
    # Eating the resource
    diet[, , no_sp + 1:no_sp+no_zoo] <- rowSums(sweep(
      fish_params@pred_kernel, 3, fish_params@dw_full * fish_params@w_full * fish_params@initial_n_other$zoo, "*"), 
      dims = 2)
  } else {
    prey <- matrix(0, nrow = no_sp + no_zoo, ncol = no_w_full)
    prey[1:no_sp, idx_sp] <- sweep(n, 2, fish_params@w * fish_params@dw, "*")
    prey[(no_sp + 1):(no_sp+no_zoo), idx_zoo] <- sweep(fish_params@initial_n_other$zoo, 2, fish_params@w_full[idx_zoo] * fish_params@dw_full[idx_zoo], "*")
    ft <- array(rep(fish_params@ft_pred_kernel_e, times = no_sp + no_zoo) *
                  rep(mvfft(t(prey)), each = no_sp),
                dim = c(no_sp, no_w_full, no_sp + no_zoo))
    # We now have an array predator x wave number x prey
    # To Fourier transform back we need a matrix of wave number x everything
    ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
    ae <- array(Re(mvfft(ft, inverse = TRUE) / no_w_full), 
                dim = c(no_w_full, no_sp, no_sp + no_zoo))
    ae <- ae[idx_sp, , , drop = FALSE]
    ae <- aperm(ae, c(2, 1, 3))
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    ae[ae < 1e-18] <- 0
    diet[, , 1:(no_sp + no_zoo)] <- ae
  }
  # Multiply by interaction matrix, including resource, and then by 
  # search volume
  inter <- cbind(fish_params@interaction, fish_params@species_params$interaction_resource)
  diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 1), drop = FALSE],
                                         c(1, 3), inter, "*"), 
                                   c(1, 2), fish_params@search_vol, "*")
  
  # Add diet from other components
  for (i in seq_along(fish_params@other_encounter)) {
    diet[, , no_sp + 1 + i] <-
      do.call(fish_params@other_encounter[[i]], 
              list(params = fish_params,
                   n = n, n_pp = n_pp, n_other = n_other,
                   component = names(fish_params@other_encounter)[[i]]))
  }
  
  # Correct for satiation and keep only entries corresponding to fish sizes
  f <- getFeedingLevel(fish_params, n, n_pp)
  fish_mask <- n > 0
  diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
  if (proportion) {
    total <- rowSums(diet, dims = 2)
    diet <- sweep(diet, c(1, 2), total, "/")
    diet[is.nan(diet)] <- 0
  }
  return(diet)
}

#' Plot diet of fish including plankton group
#'
#' @param fish_object An object of class MizerSim or MizerParams with the fish model
#' @param zoo_params (Optional) An object of class MizerPrarams with the zooplankton model parameters
#' @param species A list of  species to include in the plot
#' @param return_data A boolean value that determines whether the formatted data used for the plot 
#'  is returned instead of the plot itself. Default value is FALSE
#'
#' @return
#' @export
#'
#' @examples
plotDiet_ZooMizer <- function (fish_object, zoo_params = NULL, species = NULL, return_data = FALSE) 
{
  assert_that(is.flag(return_data))
  if (is(fish_object, "MizerSim")) {
    fish_params <- fish_object@params
    fish_params <- setInitialValues(fish_params, fish_object)
  } else if (is(fish_object, "MizerParams")) {
    fish_params <- validParams(fish_object)
  }
  if (is.null(zoo_params)){
    zoo_params <- fish_params@other_params$zoo$params
  }
  fish_species <- valid_species_arg(fish_object, species, return.logical = TRUE)
  zoo_species <- valid_species_arg(zoo_params, NULL, return.logical = TRUE)
  diet <- getDiet_ZooMizer(fish_params)[fish_species, , ]
  prey <- dimnames(diet)$prey
  prey <- factor(prey, levels = rev(prey))
  plot_dat <- data.frame(w = fish_params@w, Proportion = c(diet), 
                         Prey = rep(prey, each = length(fish_params@w)))
  plot_dat <- plot_dat[plot_dat$Proportion > 0.001, ]
  if (return_data) 
    return(plot_dat)
  no_zoo <- length(zoo_params@species_params$species)
  zoo_params@linecolour[1:no_zoo] <- zoo_params@species_params$PlotColour
  legend_levels <- intersect(names(c(fish_params@linecolour, zoo_params@linecolour)), plot_dat$Prey)
  ggplot(plot_dat) + geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
    scale_x_log10() + labs(x = "Size [g]") + scale_fill_manual(values = c(fish_params@linecolour, zoo_params@linecolour)[legend_levels])
}

#' Plot the rate of background mortality of each species against size
#'
#' @param object An object of class MizerSim or MizerParams with the fish model
#' @param species A list of  species to include in the plot
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step. Ignored when called with a \linkS4class{MizerParams}
#'   object.
#' @param all.sizes If TRUE, then predation mortality is plotted also for sizes outside a species' size range. Default FALSE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data used for the plot is returned instead of the plot itself. Default value is FALSE
#' @param ... Other arguments (currently unused)
#'
#' @return A ggplot2 object, unless return_data = TRUE, in which case a data frame with the three variables 'w', 'value', 'Species' is returned.
#' @export
#'
#' @examples 
plotBackgroundMort <- function (object, species = NULL, time_range, all.sizes = FALSE, 
                                highlight = NULL, return_data = FALSE, ...) 
{
  assert_that(is.flag(all.sizes), is.flag(return_data))
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- object@params
  }
  else {
    params <- validParams(object)
  }
  mu_b <- params@mu_b
  
  species <- valid_species_arg(params, species)
  mu_b <- mu_b[as.character(dimnames(mu_b)[[1]]) %in% 
                           species, , drop = FALSE]
  plot_dat <- data.frame(w = rep(params@w, each = length(species)), 
                         value = c(mu_b), Species = species)
  if (!all.sizes) {
    for (sp in species) {
      plot_dat$value[plot_dat$Species == sp & (plot_dat$w < 
                                                 params@species_params[sp, "w_min"] | plot_dat$w > 
                                                 params@species_params[sp, "w_inf"])] <- NA
    }
    plot_dat <- plot_dat[complete.cases(plot_dat), ]
  }
  if (return_data) 
    return(plot_dat)
  p <- plotDataFrame(plot_dat, params, xlab = "Size [g]", 
                     xtrans = "log10", highlight = highlight)
  suppressMessages(p <- p + scale_y_continuous(labels = prettyNum, 
                                               name = "Background mortality [1/year]", limits = c(0, 
                                                                                                 max(plot_dat$value))))
  p
}

#' @rdname plotBackgroundMort
#' @export
plotlyBackgroundMort <- function (object, species = NULL, time_range, highlight = NULL, 
                                  ...) 
{
  argg <- as.list(environment())
  ggplotly(do.call("plotBackgroundMort", argg), tooltip = c("Species", 
                                                      "w", "value"))
}

summary_plot <- function(object, time_range) {
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params,object)
    if (missing(time_range)) {
      time_range  <- max(as.numeric(dimnames(object@n)$time))
    }
    if (length(time_range) > 1){
    times <- max(getTimes(object))
    params@initial_n[] <- colMeans(object@n[time_range,,])
    nothers <- array(0,c(times,dim(params@initial_n_other[[1]])))
    for (i in 1:times) {
      nothers[i,,]  <- object@n_other[[i]]
    }
    params@initial_n_other$zoo[] <- colMeans(nothers[time_range,,])
    params@initial_n_pp <- colMeans(object@n_pp[time_range,])
    }
    }
  else {
    params <- validParams(object)
  }
  if (missing(time_range)) {
    time_range  <- max(as.numeric(dimnames(object@n)$time))
  }
  
  ts <- plotBiomass_ZooMizer(object)
  ss <- plotSpectra_ZooMizer(object, time_range = time_range, resource = FALSE, wlim = c(1e-14, NA))+theme(legend.position = "none")
  diet <- plotDiet_ZooMizer(params, species = tail(object@params@species_params$species, 1))+theme(legend.position = "none")+labs(y="Proportion of fish diet")
  bg <- plotBackgroundMort_ZooMizer(params)+theme(legend.position = "none")
  pred <- plotPredMort_ZooMizer(params)+theme(legend.position = "none")
  tm <- plotTotalMort_ZooMizer(params)+theme(legend.position = "none")
  gr <- plotZooGrowth(params)
  gc <- plotGrowthCurves(params@other_params$zoo$params, percentage = TRUE, max_age = 2)+
    scale_colour_manual(values =  params@other_params$zoo$params@species_params$PlotColour)+
    theme(legend.position = "none")
  surv <- plotSurvivalCurves_ZooMizer(params, max_age = 1)
  
  plot <- (ts + ss + diet) / (bg + pred + tm) / (gr + gc + surv) + plot_layout(guides = "collect")

  return(plot)
}

plotBackgroundMort_ZooMizer <- function(object, species = NULL, time_range, all.sizes = FALSE, 
                                                 highlight = NULL, return_data = FALSE, ...) {
  assert_that(is.flag(all.sizes), is.flag(return_data))
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- object@params
  }
  else {
    params <- validParams(object)
  }
  zoobg <- plotBackgroundMort(params@other_params$zoo$params, species = species, time_range = time_range, all.sizes = all.sizes, highlight = highlight,
                              return_data = TRUE)
  fishbg <- plotBackgroundMort(params, species = species, time_range = time_range, all.sizes = all.sizes, highlight = highlight,
                              return_data = TRUE)
  
  frame <- rbind(zoobg, fishbg)
  
  if (return_data)
    return(frame)
  
  # p <- plotDataFrame(plot_dat, params, xlab = "Size [g]", 
  #                    xtrans = "log10", highlight = highlight)
  
  var_names <- names(frame)
  x_var <- var_names[[1]]
  y_var <- var_names[[2]]
  group_var <- var_names[[3]]
  frame$Legend <- frame[[group_var]]
  legend_var <- "Legend"
  legend_levels <- intersect(c(names(params@other_params$zoo$params@linecolour),names(params@linecolour)), frame[[legend_var]])
  frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
  zoocols <- params@other_params$zoo$params@species_params$PlotColour
  names(zoocols) <- params@other_params$zoo$params@species_params$species
  linecolour <- c(zoocols,params@linecolour)[legend_levels]
  linetype <- c(params@other_params$zoo$params@linetype,params@linetype)[legend_levels]
  linesize <- rep_len(0.8, length(legend_levels))
  names(linesize) <- legend_levels
  linesize[highlight] <- 1.6
  xbreaks <- log_breaks()
  ybreaks <- waiver()
  
  p <- ggplot(frame, aes(group = .data[[group_var]])) +
    scale_y_continuous(trans = "identity", breaks = ybreaks, labels = prettyNum, name = waiver()) + 
    scale_x_continuous(trans = "log10", name = waiver()) +
    scale_colour_manual(values = linecolour) + 
    scale_linetype_manual(values = linetype) + scale_size_manual(values = linesize) + 
    geom_line(aes(x = .data[[x_var]], y = .data[[y_var]], 
                  colour = .data[[legend_var]], linetype = .data[[legend_var]], 
                  size = .data[[legend_var]]))
  
  suppressMessages(p <- p + scale_y_continuous(labels = prettyNum, 
                                               name = "Background mortality [1/year]",
                                               limits = c(0,max(frame$value))))
  p
}

getPredMort_ZooMizer <- function(object) {
  if (is(object, "MizerSim")) {
        params <- setInitialValues(object@params, object)
        }
      else {
        params <- validParams(object)
      }
 
  total_mort_from_fish <- getResourceMort(params)
  
  zoo_params <- params@other_params$zoo$params
  zoo_idx <- (length(zoo_params@w_full) - length(zoo_params@w) + 1):length(zoo_params@w_full)
  fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  
  no_grps_in_size_class <- colSums(params@initial_n_other$zoo > 0)
  no_grps_in_size_class[which(no_grps_in_size_class == 0)] <- 1
  no_grps_in_size_class[which(is.na(no_grps_in_size_class))] <- 1
  
  if(any(no_grps_in_size_class == 0)) {stop("NaNs produced by mort / # of groups") }
  
  n_eff <- colSums(sweep(params@initial_n_other$zoo, 1,
                         zoo_params@species_params$Carbon * zoo_params@species_params$GrossGEscale,
                         "*")) / params@species_params$alpha[1]
  
  scaling <- sapply(as.data.frame(rbind(colSums(params@initial_n_other$zoo) / n_eff, 0)), max, na.rm = TRUE) # re-scale by carbon content, set NAs to 0.
  resmort <- matrix(total_mort_from_fish[zoo_idx] / no_grps_in_size_class * scaling, byrow = TRUE, nrow = nrow(zoo_params@species_params), ncol = length(zoo_params@w))
  predmort <- resmort + getPredMort(params@other_params$zoo$params)
  
  fpm <- getPredMort(params, drop = FALSE)
  
  
  
  fish_species <- valid_species_arg(params, NULL)
  zoo_species <- valid_species_arg(zoo_params, NULL)
  species <- union(zoo_species, fish_species)
  
  for (sp in 1:length(fish_species)) {
    fpm[sp,params@w < params@species_params[sp, "w_min"] | params@w > params@species_params[sp, "w_inf"]] <- 0
  }
  
  for (sp in 1:length(zoo_species)) {
    predmort[sp,zoo_params@w < zoo_params@species_params[sp, "w_min"] | zoo_params@w > zoo_params@species_params[sp, "w_inf"]] <- 0
  }
  
  pred_mort <- matrix(0, nrow=length(species), ncol = length(zoo_params@w_full))
  dimnames(pred_mort)[[1]] <- as.character(species)
  dimnames(pred_mort)[[2]] <- as.character(params@w_full)
  pred_mort[1:length(zoo_species), zoo_idx] <- predmort
  pred_mort[(length(zoo_species)+1):length(species), fish_idx] <- fpm
  
  return(pred_mort)
}

plotPredMort_ZooMizer <- function (object, species = NULL, time_range, all.sizes = FALSE,
                                   highlight = NULL, return_data = FALSE, ...)
{
  assert_that(is.flag(all.sizes), is.flag(return_data))
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
      params <- setInitialValues(object@params, object)
    }
  } else {
    params <- validParams(object)
    }

pred_mort <- getPredMort_ZooMizer(object)

zoo_params <- params@other_params$zoo$params
zoo_idx <- (length(zoo_params@w_full) - length(zoo_params@w) + 1):length(zoo_params@w_full)
fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

  if (length(dim(pred_mort)) == 3) {
    pred_mort <- apply(pred_mort, c(2, 3), mean)
  }
  fish_species <- valid_species_arg(params, species)
  zoo_species <- valid_species_arg(zoo_params, NULL)
  species <- union(zoo_species, fish_species)
  pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in%
                           species, , drop = FALSE]

  plot_dat <- data.frame(w = rep(params@w_full, each = length(species)),
                         value = c(pred_mort), Species = species)
  if (!all.sizes) {
    specieslist <- merge(params@species_params, params@other_params$zoo$params@species_params, all = TRUE)
    for (sp in species) {
      plot_dat$value[plot_dat$Species == sp & (plot_dat$w <
                                                 specieslist[sp, "w_min"] | plot_dat$w >
                                                 specieslist[sp, "w_inf"])] <- NA
    }
    plot_dat <- plot_dat[complete.cases(plot_dat), ]
  }
  if (return_data)
    return(plot_dat)
  
  frame <- plot_dat
  
  var_names <- names(frame)
  x_var <- var_names[[1]]
  y_var <- var_names[[2]]
  group_var <- var_names[[3]]
  frame$Legend <- frame[[group_var]]
  legend_var <- "Legend"
  legend_levels <- intersect(c(names(params@other_params$zoo$params@linecolour),names(params@linecolour)), frame[[legend_var]])
  frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
  zoocols <- params@other_params$zoo$params@species_params$PlotColour
  names(zoocols) <- params@other_params$zoo$params@species_params$species
  linecolour <- c(zoocols,params@linecolour)[legend_levels]
  linetype <- c(params@other_params$zoo$params@linetype,params@linetype)[legend_levels]
  linesize <- rep_len(0.8, length(legend_levels))
  names(linesize) <- legend_levels
  linesize[highlight] <- 1.6
  xbreaks <- log_breaks()
  ybreaks <- waiver()
  
  p <- ggplot(frame, aes(group = .data[[group_var]])) +
    scale_y_continuous(trans = "identity", breaks = ybreaks, labels = prettyNum, name = waiver()) + 
    scale_x_continuous(trans = "log10", name = waiver()) +
    scale_linetype_manual(values = linetype) + scale_size_manual(values = linesize) + 
    geom_line(aes(x = .data[[x_var]], y = .data[[y_var]], 
                  colour = .data[[legend_var]], linetype = .data[[legend_var]], 
                  size = .data[[legend_var]]))+
    scale_colour_manual(values = linecolour)
  
  suppressMessages(p <- p + scale_y_continuous(labels = prettyNum, 
                                               name = "Predation mortality [1/year]",
                                               limits = c(0,max(plot_dat$value))))
  p
  }

plotTotalMort_ZooMizer <- function (object, species = NULL, time_range, all.sizes = FALSE,
                                   highlight = NULL, return_data = FALSE, ...)
{
  assert_that(is.flag(all.sizes), is.flag(return_data))
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
      params <- setInitialValues(object@params, object)
    }
  } else {
    params <- validParams(object)
  }
  
  pred <- plotPredMort_ZooMizer(object, return_data = TRUE, all.sizes = TRUE) %>% 
    reshape2::acast(Species ~ w)
  zoo_idx <- which(params@w_full >= min(params@other_params$zoo$params@w))
  pred <- pred[,zoo_idx]
  pred[is.na(pred)] <- 0
  
  bg <- plotBackgroundMort_ZooMizer(object, return_data = TRUE, all.sizes = TRUE) %>% 
    reshape2::acast(Species ~ w)
    bg[is.na(bg)] <- 0
  
  frame <- pred + bg
  frame <- melt(frame)
  
  names(frame) <- list("Species", "w", "Mort")
  
  var_names <- names(frame)
  x_var <- var_names[[2]]
  y_var <- var_names[[3]]
  group_var <- var_names[[1]]
  frame$Legend <- frame[[group_var]]
  legend_var <- "Legend"
  legend_levels <- intersect(c(names(params@other_params$zoo$params@linecolour),names(params@linecolour)), frame[[legend_var]])
  frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
  zoocols <- params@other_params$zoo$params@species_params$PlotColour
  names(zoocols) <- params@other_params$zoo$params@species_params$species
  linecolour <- c(zoocols,params@linecolour)[legend_levels]
  linetype <- c(params@other_params$zoo$params@linetype,params@linetype)[legend_levels]
  linesize <- rep_len(0.8, length(legend_levels))
  names(linesize) <- legend_levels
  linesize[highlight] <- 1.6
  xbreaks <- log_breaks()
  ybreaks <- waiver()
  
  p <- ggplot(frame, aes(group = .data[[group_var]])) +
    scale_y_continuous(trans = "identity", breaks = ybreaks, labels = prettyNum, name = waiver()) + 
    scale_x_continuous(trans = "log10", name = waiver()) +
    scale_linetype_manual(values = linetype) + scale_size_manual(values = linesize) + 
    geom_line(aes(x = .data[[x_var]], y = .data[[y_var]], 
                  colour = .data[[legend_var]], linetype = .data[[legend_var]], 
                  size = .data[[legend_var]]))+
    scale_colour_manual(values = linecolour)

  suppressMessages(p <- p + scale_y_continuous(labels = prettyNum, 
                                               name = "Total mortality [1/year]",
                                               limits = c(0,max(frame$Mort))))
  p
}


getZooMort <- function(params, n = initialN(params), n_pp = initialNResource(params), 
                       n_other = initialNOther(params), effort = getInitialEffort(params), 
                       t = 0, ...) {

zoo_params <- params@other_params$zoo$params

# get predation mortality imposed by the fish - note that this is scaled by carbon content and must be rescaled
total_mort_from_fish <- getResourceMort(params)

zoo_idx <- (length(zoo_params@w_full) - length(zoo_params@w) + 1):length(zoo_params@w_full)
fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

no_grps_in_size_class <- colSums(n_other$zoo > 0)
no_grps_in_size_class[which(no_grps_in_size_class == 0)] <- 1
no_grps_in_size_class[which(is.na(no_grps_in_size_class))] <- 1

if(any(no_grps_in_size_class == 0)) {stop("NaNs produced by mort / # of groups") }

n_eff <- colSums(sweep(n_other$zoo, 1, 
                       zoo_params@species_params$Carbon * zoo_params@species_params$GrossGEscale, 
                       "*")) / params@species_params$alpha[1]

scaling <- sapply(as.data.frame(rbind(colSums(n_other$zoo) / n_eff, 0)), max, na.rm = TRUE) # re-scale by carbon content, set NAs to 0.
mort_from_fish <- matrix(total_mort_from_fish[zoo_idx] / no_grps_in_size_class * scaling, byrow = TRUE, nrow = nrow(zoo_params@species_params), ncol = length(zoo_params@w))

# add mortality of fish eating zoo to the external mortality. Better if this was done as a fishing mortality?
return(getMort(zoo_params) + mort_from_fish)
}

getGrowthCurves_ZooMizer <- function (object, species = NULL, max_age = 10, percentage = FALSE) 
{
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getGrowthCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }
  zoo_params <- params@other_params$zoo$params
  
  species <- valid_species_arg(zoo_params, species)
  idx <- which(zoo_params@species_params$species %in% species)
  species <- zoo_params@species_params$species[idx]
  age <- seq(0, max_age, length.out = 50)
  dt <- rep(age[2]-age[1], 50)
  ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
                                                                     Age = age))
  # get vector of phytoplankton and fish resource
  total_fish_n <- params@w_full*0
  total_fish_n[fish_idx] <- colSums(n)
  
  g <- getEGrowth(zoo_params, n_pp = zoo_params@initial_n_pp + total_fish_n)
  m <- getZooMort(params)
  for (j in seq_along(species)) {
    i <- idx[j]
    g_fn <- stats::approxfun(c(log(zoo_params@w), log(zoo_params@species_params$w_inf[[i]])), 
                             c(g[i, ], 0))
    myodefun <- function(t, state, parameters) {
      return(list(g_fn(state)))
    }
    ws[j, ] <- deSolve::ode(y = log(zoo_params@w[zoo_params@w_min_idx[i]]), 
                            times = age, func = myodefun)[, 2]
    if (percentage) {
      ws[j, ] <- ws[j, ]/zoo_params@species_params$w_inf[j] * 100
    }
  }
  
  return(ws)
}

getSurvivalCurves_ZooMizer<- function(object, species = NULL, max_age = 10, percentage = FALSE) 
{
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getGrowthCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }
  zoo_params <- params@other_params$zoo$params
  
  species <- valid_species_arg(zoo_params, species)
  idx <- which(zoo_params@species_params$species %in% species)
  species <- zoo_params@species_params$species[idx]
  age <- seq(0, max_age, length.out = 50)
  # dt <- rep(age[2]-age[1], 50)
  ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
                                                                     Age = age))
  # get vector of phytoplankton and fish resource
  total_fish_n <- params@w_full*0
  fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  total_fish_n[fish_idx] <- colSums(params@initial_n)
  
  g <- getEGrowth(zoo_params, n_pp = zoo_params@initial_n_pp + total_fish_n)
  m <- getZooMort(params)
  for (j in seq_along(species)) {
    i <- idx[j]
    m_fn <- stats::approxfun(c(zoo_params@w, zoo_params@species_params$w_inf[[i]]), 
                              c(m[i, ], 0))
    
    survival <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <- - m_fn(w) * S
        list(c(w,dS))
      })
    }
    
    ws[j, ] <- deSolve::ode(y=c(w=zoo_params@species_params$w_min[j], S=1),
                            times = age, func = survival,
                            parms = 0)[,3]
  }
  # if (percentage) {
  #   ws[, ] <- ws[, ]/ws[,1] * 100
  # }
  
  return(ws)
}

# getSurvivaltosize_ZooMizer <- function(object, species = NULL, max_age = 10, percentage = FALSE) 
#   {
#     if (is(object, "MizerSim")) {
#       params <- object@params
#       params <- setInitialValues(params, object)
#     }
#     else if (is(object, "MizerParams")) {
#       params <- validParams(object)
#     }
#     else {
#       stop("The first argument to `getGrowthCurves()` must be a ", 
#            "MizerParams or a MizerSim object.")
#     }
#     zoo_params <- params@other_params$zoo$params
#     
#     species <- valid_species_arg(zoo_params, species)
#     idx <- which(zoo_params@species_params$species %in% species)
#     species <- zoo_params@species_params$species[idx]
#     age <- seq(0, max_age, length.out = 50)
#     # dt <- rep(age[2]-age[1], 50)
#     ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
#                                                                        Weight = w))
#     # get vector of phytoplankton and fish resource
#     total_fish_n <- params@w_full*0
#     total_fish_n[fish_idx] <- colSums(n)
#     
#     g <- getEGrowth(zoo_params, n_pp = zoo_params@initial_n_pp + total_fish_n)
#     m <- getZooMort(params)
#     for (j in seq_along(species)) {
#       i <- idx[j]
#       mg_fn <- stats::approxfun(c(zoo_params@w, zoo_params@species_params$w_inf[[i]]), 
#                                 c(-m[i, ]/g[i, ], 0))
#       
#       survival <- function(t, state, parameters) {
#         with(as.list(c(state, parameters)), {
#           dS <- mg_fn(w) * S
#           list(c(w,dS))
#         })
#       }
#       
#       ws[j, ] <- deSolve::ode(y=c(w=zoo_params@species_params$w_min[j], S=1),
#                               times = age, func = survival,
#                               parms = 0)[,3]
#     }
#     if (percentage) {
#     ws <- ws * 100
# 
#     return(ws)
#   }

    
plotSurvivalCurves_ZooMizer <- function(object, species = NULL, max_age = 10, percentage = FALSE, return_data = FALSE, highlight = NULL, species_panel = FALSE){
  
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getGrowthCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }
  zoo_params <- params@other_params$zoo$params
  species <- valid_species_arg(zoo_params, species)
  sp_sel <- zoo_params@species_params$species %in% species
    
  ws <- getSurvivalCurves_ZooMizer(object, species, max_age, percentage)
  
  colours <- zoo_params@species_params$PlotColour[sp_sel]
  names(colours) <- zoo_params@species_params$species

  plot_dat <- melt(ws)
  plot_dat$Species <- factor(plot_dat$Species, zoo_params@species_params$species)
  plot_dat$Legend <- "model"
  
  if (return_data) 
    return(plot_dat)
  p <- ggplot(filter(plot_dat, Legend == "model")) + 
    geom_line(aes(x = Age, y = value, colour = Species, linetype = Species, 
                  size = Species))
  y_label <- if (percentage) 
    "Percent of cohort remaining"
  else "Fraction of cohort remaining"
  
  legend_levels <- intersect(c(dimnames(zoo_params@initial_n)$sp, 
                               "Background", "Resource", "Total"), 
                             plot_dat$Species)
  plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
  linesize <- rep(0.8, length(legend_levels))
  names(linesize) <- names(zoo_params@linetype[legend_levels])
  linesize[highlight] <- 1.6
  p <- p + scale_x_continuous(name = "Age [Years]") + 
    scale_y_continuous(name = y_label) + scale_colour_manual(values = colours[legend_levels]) + 
    scale_linetype_manual(values = zoo_params@linetype[legend_levels]) + 
    scale_size_manual(values = linesize)
  if (species_panel) {
    p <- ggplot(plot_dat) + geom_line(aes(x = Age, y = value, 
                                          colour = Legend)) + scale_x_continuous(name = "Age [years]") + 
      scale_y_continuous(name = y_label) + 
      facet_wrap(~Species, scales = "free_y")
  }
  return(p)
}


getSurvivalCurves <- function(object, species = NULL, max_age = 5, percentage = FALSE){
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getSurvivalCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }

  species <- valid_species_arg(params, species)
  idx <- which(params@species_params$species %in% species)
  species <- params@species_params$species[idx]
  age <- seq(0, max_age, length.out = 50)
  ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
                                                                     Age = age))
  m <- getMort(params)
  for (j in seq_along(species)) {
    i <- idx[j]
    m_fn <- stats::approxfun(c(params@w, params@species_params$w_inf[[i]]), 
                             c(m[i, ], 0))
    
    survival <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <- - m_fn(w) * S
        list(c(w,dS))
      })
    }
    
    ws[j, ] <- deSolve::ode(y=c(w=params@species_params$w_min[j], S=1),
                            times = age, func = survival,
                            parms = 0)[,3]
  }
  if (percentage) {
     ws[, ] <- ws[, ] * 100
  }
  
  return(ws)
}

plotSurvivalCurves <- function(object, species = NULL, max_age = 10, percentage = FALSE, return_data = FALSE, highlight = NULL, species_panel = FALSE){
  
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getSurvivalCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }
  species <- valid_species_arg(params, species)
  sp_sel <- params@species_params$species %in% species
  ws <- getSurvivalCurves(object, species, max_age, percentage)

  plot_dat <- reshape2::melt(ws)
  plot_dat$Species <- factor(plot_dat$Species, params@species_params$species)
  plot_dat$Legend <- "model"
  
  if (return_data) 
    return(plot_dat)
  p <- ggplot(filter(plot_dat, Legend == "model")) + 
    geom_line(aes(x = Age, y = value, colour = Species, linetype = Species, 
                  size = Species))
  y_label <- if (percentage) 
    "Percent of cohort remaining"
  else "Fraction of cohort remaining"
  
  legend_levels <- intersect(c(dimnames(params@initial_n)$sp, 
                               "Background", "Resource", "Total"), 
                             plot_dat$Species)
  plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
  linesize <- rep(0.8, length(legend_levels))
  names(linesize) <- names(params@linetype[legend_levels])
  linesize[highlight] <- 1.6
  p <- p + scale_x_continuous(name = "Age [Years]") + 
    scale_y_continuous(name = y_label) + scale_colour_manual(values = params@linecolour[legend_levels]) + 
    scale_linetype_manual(values = params@linetype[legend_levels]) + 
    scale_size_manual(values = linesize)
  if (species_panel) {
      p <- ggplot(plot_dat) + geom_line(aes(x = Age, y = value, 
                                            colour = Legend)) + scale_x_continuous(name = "Age [years]") + 
        scale_y_continuous(name = y_label) + 
        facet_wrap(~Species, scales = "free_y")
    }
  return(p)
}

plotZooGrowth <- function(params, byweight= TRUE, carbon = FALSE){
  growths <- getEGrowth(params@other_params$zoo$params)
  ylab = "Growth [g/day]"
  if(byweight){
    ws <- matrix(data = params@other_params$zoo$params@w, nrow = 9, ncol = 178, byrow=T)
    growths <- growths/as.array(ws)
    ylab = "Growth [1/day]"
    }
  #/params@other_params$zoo$params@species_params$Carbon
  for(i in 1:9) growths[i,1:params@other_params$zoo$params@w_min_idx[i]] <- 0
  growths <- melt(growths) %>% rename(Species = sp)
  growths <- growths[which(growths$value > 0),]
  colours <- params@other_params$zoo$params@species_params$PlotColour
  
  linesize <- rep(1.6, length(params@other_params$zoo$params@species_params$species))
  names(linesize) <- names(params@other_params$zoo$params@linetype[1:9])
  
  ggplot(growths, aes(x=w, y=value/365))+
    geom_line(aes(colour = Species))+
    scale_x_log10()+
    scale_color_manual(values = colours)+
    scale_linetype_manual(values = c(params@other_params$zoo$params@linetype[1:nrow(params@other_params$zoo$params@species_params)]))+
    scale_size_manual(values = linesize)+
    labs(y = ylab, x = "Weight [g]")
}

