check_species_df <- function(species_df, microclimate_opts = no_microclimate_opts()) {

  if (nrow(species_df) < 1) {
    stop("Species trait table is empty.")
  }

  exptd_params <- species_trait_names(microclimate_opts)
  exptd_params <- c("SpeciesID", exptd_params)

  missing_params <- exptd_params[!exptd_params %in% names(species_df)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("species_df", missing_params))
  }

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(species_df, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- exptd_params[is_missing_val]
    stop(paste(c(
      "The following elements of species_df contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))

  }

  is_numeric <- sapply(species_df, function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- exptd_params[!is_numeric]
    stop(paste(c(
      "The following elements of species_df contain non-numeric values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  positive_params <- exptd_params[!exptd_params %in% c(
    "LightResponseA", "LightResponseB", "LightResponseC"
    )]
  is_positive <- sapply(species_df[positive_params], function(x) all(x >= 0))
  if (any(!is_positive)) {
    wrong_params <- positive_params[!is_positive]
    stop(paste(c(
      "The following elements of species_df contain negative values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (any(duplicated(species_df$SpeciesID))) {
    stop("species_df contains multiple entries for the same species")
  }

  if (any(species_df$MaximumMass < species_df$MassAtMaturity)) {
    stop("All values in species_df column MaximumMass should be >= MassAtMaturity")
  }

  if (any(species_df$RecruitmentInvestmentRel > 1)) {
    stop("species_df element RecruitmentInvestmentRel contains elements higher than 1.")
  }

  if (any(species_df$RecruitmentInc > 1)) {
    stop("species_df element RecruitmentInc contains elements higher than 1.")
  }

  if (any(species_df$DispersalKernelAsymmetry > 1)) {
    stop("species_df element DispersalKernelAsymmetry contains elements higher than 1.")
  }

  if (any(species_df$OptimumLight < species_df$MinLight)) {
    stop("All values in species_df column OptimumLight should be >= MinLight")
  }

  if (any(species_df$MaxLight < species_df$OptimumLight)) {
    stop("All values in species_df column MaxLight should be >= MinLight")
  }

  if (any(species_df$MaxLight < species_df$OptimumLight)) {
    stop("All values in species_df column MaxLight should be >= MinLight")
  }

  if (microclimate_opts$use_wind) {
    if (any(species_df$DispersalKernelWindEffect < 0 ||
            species_df$DispersalKernelWindEffect > 1)) {
      stop("DispersalKernelWindEffect should be between 0 and 1.")
    }
  }

}
