check_sim_params <- function(sim_params) {

  exptd_params <- c(
    "InitialTimeStep", "timeSteps", "StopCriterionHa","hasDynamicMicrohabitat",
    "massDepCompetition", "use_mass_dep_mortality", "use_wind_dispersal",
    "SurfaceBiomassScaling", "SlopeRecruitment", "InterceptRecruitment"
  )

  if (!is.logical(sim_params$massDepCompetition)) {
    stop("massDepCompetition must be TRUE or FALSE.")
  }
  if (sim_params$use_mass_dep_mortality) {
    exptd_params <- c(exptd_params, "MortRateMass", "MortRateMassScaling")
  } else {
    exptd_params <- c(exptd_params, "MortRateRandom")
  }

  missing_params <- exptd_params[!exptd_params %in% names(sim_params)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("sim_params", missing_params))
  }

  sim_params <- sim_params[exptd_params]
  param_names <- exptd_params

  # Check for missing values
  is_missing_val <- sapply(sim_params, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- param_names[is_missing_val]
    stop(paste(c(
      "The following elements of sim_params contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (sim_params$timeSteps <= 0) {
    "timeSteps must be 1 or more."
  }

  if (sim_params$StopCriterionHa < 10000) {
    stop("StopCriterionHa cannot be smaller than 10000 (i.e, 1 epiphyte per voxel)")
  }

  if (!is.logical(sim_params$hasDynamicMicrohabitat)) {
    stop("hasDynamicMicrohabitat must be TRUE or FALSE.")
  }

  if (!is.logical(sim_params$use_wind_dispersal)) {
    stop("use_wind_dispersal must be TRUE or FALSE.")
  }

  if (!is.logical(sim_params$use_mass_dep_mortality)) {
    stop("use_mass_dep_mortality must be TRUE or FALSE.")
  }

  if (sim_params$use_mass_dep_mortality) {
    if (sim_params$MortRateMass < 0 || sim_params$MortRateMass > 1) {
      stop("MortRateMass must be between 0 and 1")
    }
    if (sim_params$MortRateMassScaling >= 0) {
      stop("MortRateMassScaling must be negative or zero.")
    }
  } else {
    if (sim_params$MortRateRandom < 0 || sim_params$MortRateRandom > 1) {
      stop("MortRateRandom must be between 0 and 1")
    }
  }

  if (sim_params$SurfaceBiomassScaling <= 0) {
    stop("SurfaceBiomassScaling must be a strictly positive number.")
  }

  if (sim_params$SlopeRecruitment < 0 || sim_params$SlopeRecruitment > 1) {
    stop("SlopeRecruitment must be between 0 and 1")
  }

  if (sim_params$InterceptRecruitment < 0) {
    stop("InterceptRecruitment must be a positive number.")
  }
}
