index_3d <- function(x, y, z, dimX, dimY) {
  return(x + dimX * (y - 1) + dimX * dimY * (z - 1))
}

#' Generate a data frame containing species traits
#'
#' Convenience wrapper that calls [draw_species_traits()] for multiple species
#' and formats the output as a dataframe.
#'
create_rnd_species_df <- function(nb_species, species_params) {

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
    MaxMassRandom = runif(1, 0, max_val), #
    InterceptAgeMaturity = runif(1, 0, max_val),
    AgeAtMaturityDevCorr = runif(1, 0, 1),
    ScalingAgeMaturity = runif(1, 0, 1),
    MassAtMaturityRelativeRandom = runif(1, 0, 1), ##

    CorrelationMassRecruitment = runif(1) > 0.5,

    DispersalKernelRandom = runif(1, 0, max_val), #
    DispersalKernelAsymmetryRandom = runif(1, 0, 1), ##

    HeightBreadthRandom = runif(1, 0, max_val), #

    Imax = runif(1, 0, max_val),
    kL = runif(1, 0, max_val),
    LAI = runif(1, 0, max_val)
  )

  sp_params$MaxMassRandom[2] <- draw_max_value(sp_params$MaxMassRandom[1], max_val)
  sp_params$MassAtMaturityRelativeRandom[2] <- draw_max_value(
    sp_params$MassAtMaturityRelativeRandom[1], 1
  )
  sp_params$DispersalKernelRandom[2] <- draw_max_value(sp_params$DispersalKernelRandom[1], max_val)
  sp_params$DispersalKernelAsymmetryRandom[2] <- draw_max_value(
    sp_params$DispersalKernelAsymmetryRandom[1], 1
  )
  sp_params$HeightBreadthRandom[2] <- draw_max_value(sp_params$HeightBreadthRandom[1], max_val)

  # Parameters that depend on the correlation mass option
  if (sp_params$CorrelationMassRecruitment) {
    sp_params$RecruitmentInvestmentRelDevCorr <- runif(1, 0, 1)
    sp_params$RecruitmentInvestmentRelMeanCorr <- runif(1, 0, max_val)
  }
  else {
    sp_params$RecruitmentIncRandom <- runif(1, 0, max_val)
    sp_params$RecruitmentIncRandom[2] <- draw_max_value(sp_params$RecruitmentIncRandom[1], max_val)

    sp_params$RecruitmentInvestmentRelMeanRandom <- runif(1, 0, max_val)
    sp_params$RecruitmentInvestmentRelMeanRandom[2] <- draw_max_value(sp_params$RecruitmentInvestmentRelMeanRandom[1], max_val)
  }
  return(sp_params)
}

# Draw second value for uniform parameters
draw_max_value <- function(min_val, upper_bound) {
  return(min_val + runif(1, 0, upper_bound - min_val))
}

draw_rnd_initial_inds_params <- function() {
  return(list(
    "SurfaceBiomassScaling" = runif(1, 1e-7, 10),
    "IndividualsPerSpecies" = sample(1:100, 1),
    "ScalingPerHa" = FALSE,
    "PercentageMaturePerSpecies" = runif(1, 0, 100)
  ))
}
