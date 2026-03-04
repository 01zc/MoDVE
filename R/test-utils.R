#' Get the sequential index of a coordinate in a 3D matrix
#'
#' Convert an XYZ position into the sequential position (i.e., ranging between
#' 1 and dimX x dimY x dimZ) inside the matrix
#'
#' @param x x-coordinate
#' @param y y-coordinate
#' @param z z-coordinate
#' @param dimX X-dimension of the matrix
#' @param dimY Y-dimension of the matrix
#'
index_3d <- function(x, y, z, dimX, dimY) {
  return(x + dimX * (y - 1) + dimX * dimY * (z - 1))
}

#' Generate a data frame containing species traits
#'
#' Convenience wrapper that calls [draw_species_traits()] for multiple species
#' and formats the output as a dataframe.
#'
create_rnd_species_df <- function(nb_species, species_params = draw_rnd_species_params()) {

  species_df <- purrr::imap_dfr(1:nb_species, function(x, i) {
    as.data.frame(c("SpeciesID" = i, draw_species_traits(species_params)))
  })

  return(species_df)
}

#' Shorthand function for generating species traits hyperparameters
#'
#' Generates a list containing all parameters required to run
#' [draw_species_traits()] with random but legal values.
#'
#' @param max_val sets the maximum value to sample hyperparameters that otherwise
#' don't have a theoretical maximum value (the minimum is zero).
#'
draw_rnd_species_params <- function(max_val = 100) {

  sp_params <- list(
    MaxMassLogScaleRandom = runif(1) > 0.5, # coin flip
    MaxMassRandom = runif(1, 0, max_val),
    InterceptAgeMaturity = runif(1, 0, max_val),
    AgeAtMaturityDevCorr = runif(1, 0, 1),
    ScalingAgeMaturity = runif(1, 0, 1),
    MassAtMaturityRelativeRandom = runif(1, 0, 1),

    CorrelationMassRecruitment = runif(1) > 0.5,

    DispersalKernelRandom = runif(1, 0, max_val),
    DispersalKernelAsymmetryRandom = runif(1, 0, 1),

    HeightBreadthRandom = runif(1, 0, 1),

    Imax = runif(1, 0, max_val),
    kL = runif(1, 0, 1),
    LAI = runif(1, 0, min(10, max_val))
  )

  sp_params$MaxMassRandom[2] <- draw_max_value(sp_params$MaxMassRandom[1], max_val)
  sp_params$MassAtMaturityRelativeRandom[2] <- draw_max_value(
    sp_params$MassAtMaturityRelativeRandom[1], 1
  )
  sp_params$DispersalKernelRandom[2] <- draw_max_value(sp_params$DispersalKernelRandom[1], max_val)
  sp_params$DispersalKernelAsymmetryRandom[2] <- draw_max_value(
    sp_params$DispersalKernelAsymmetryRandom[1], 1
  )
  sp_params$HeightBreadthRandom[2] <- draw_max_value(sp_params$HeightBreadthRandom[1], 1)

  # Parameters that depend on the correlation mass option
  if (sp_params$CorrelationMassRecruitment) {
    sp_params$RecruitmentInvestmentRelDevCorr <- runif(1, 0, 1)
    sp_params$RecruitmentInvestmentRelMeanCorr <- runif(1, 0, 1)
  }
  else {
    sp_params$RecruitmentIncRandom <- runif(1, 0, 1)
    sp_params$RecruitmentIncRandom[2] <- draw_max_value(sp_params$RecruitmentIncRandom[1], 1)

    sp_params$RecruitmentInvestmentRelMeanRandom <- runif(1, 0, 1)
    sp_params$RecruitmentInvestmentRelMeanRandom[2] <- draw_max_value(sp_params$RecruitmentInvestmentRelMeanRandom[1], 1)
  }
  return(sp_params)
}

#' Draw max value for uniform parameters
#'
#' Given the minimum value and an upper bound, draw a random max value such that
#' the max value is larger than the min and smaller than the upper bound
#'
#' @param min_val the minimum value of the uniform distribution
#' @param upper_bound upper bound to not exceed when drawing the value
#'
draw_max_value <- function(min_val, upper_bound) {
  return(min_val + runif(1, 0, upper_bound - min_val))
}

#' Draw random values for individual initialisation parameters
#'
#' Returns a list with random values for [draw_initial_individuals()] argument
#' `distr_params`
#'
#' @return a list with four elements: `SurfaceBiomassScaling`,
#' `IndividualsPerSpecies`, `ScalingPerHa`, and `PercentageMaturePerSpecies`
#'
draw_rnd_initial_inds_params <- function() {
  return(list(
    "SurfaceBiomassScaling" = runif(1, 1e-7, 10),
    "IndividualsPerSpecies" = sample(1:100, 1),
    "ScalingPerHa" = FALSE,
    "PercentageMaturePerSpecies" = runif(1, 0, 100)
  ))
}

#' Create a microhabitat matrix of specified dimensions and fill its
#' surface area and light layer based on species requirements
#'
#' Microhabitat is set to be suitable for the *first* species in SpeciesPool
#'
#' Surface area is set such that:
#' * a random fraction of all voxels are suitable
#' * each suitable voxels can sustain between 1 and 10 individuals of maximum mass
create_rnd_microhabitat <- function(SpeciesPool, dimensions, SurfaceBiomassScaling) {

  Microhabitat <- array(0, c(dimensions, 3))
  nb_suitable_voxels <- round(prod(dimensions) * runif(1, 0, 1))
  suitable_voxels <- sample(1:prod(dimensions), nb_suitable_voxels)
  reqd_sa_per_ind <- SpeciesPool$MaximumMass[1]^(2/3) / SurfaceBiomassScaling
  surface_area_mat <- array(0, dim = dimensions)
  surface_area_mat[suitable_voxels] <- reqd_sa_per_ind *
    sample(1:10, length(suitable_voxels), replace = TRUE)
  Microhabitat[, , , 1] <- surface_area_mat
  Microhabitat[, , , 3] <- SpeciesPool$OptimumLight[1]

  return(Microhabitat)
}
