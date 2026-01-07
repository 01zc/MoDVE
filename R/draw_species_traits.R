#' Draw species traits from generating parameters
#'
#' @param species_params a list of trait-generating parameter, which must
#' contain the following elements:
#'
#'  * `MaxMassLogScaleRandom` logical, if `TRUE` then `MaxMass` is sampled in
#'  a uniform between the log10 of both `MaxMassRandom` values.
#'  * `MaxMassRandom` a length-2 numeric vector containing the minimum and
#'  maximum MaxMass
#'  * `MassAtMaturityRelativeRandom` a length-2 numeric vector containing the
#'  minimum and maximum bounds (between 0 and 1) in which to sample what fraction
#'  of MaxMass the mass at maturity corresponds to.
#'
#'  * `InterceptAgeMaturity` numeric, the intercept of the function generating
#'  the age at maturity as a function of mass at maturity
#'  * `ScalingAgeMaturity` numeric, the exponent of the function generating
#'  the age at maturity as a function of mass at maturity
#'  * `AgeAtMaturityDevCorr` numeric, the range of possible deviation
#'  coefficients in the function generating age at maturity from mass at maturity.
#'
#'  * `CorrelationMassRecruitment` logical, controls the option used to
#'  generate the parameters of the mass-to-reproduction allocation function.
#'   If `TRUE`:
#'   * `RecruitmentInvestmentRel` is sampled in a uniform distribution
#'   with min and max values \eqn{`RecruitmentInvestmentRelMeanCorr` *
#'   (1 +- `RecruitmentInvestmentRelDevCorr`)},
#'   * `RecruitmentInc` is set to 0.
#'   If `FALSE`:
#'   * `RecruitmentInvestmentRel` is sampled in a uniform distribution
#'   with min and max values `RecruitmentInvestmentRelMeanRandom`
#'   * `RecruitmentInc` is sampled in a uniform distribution with min and max
#'   values `RecruitmentIncRandom.`
#'
#'  * `RecruitmentInvestmentRelMeanCorr` see above
#'  * `RecruitmentInvestmentRelDevCorr` see above
#'  * `RecruitmentInvestmentRelMeanRandom` see above
#'  * `RecruitmentIncRandom` see above
#'
#'  * `DispersalKernelRandom` a length-2 numeric vector, min and max possible
#'  values of the dispersal kernel.
#'  * `DispersalKernelAsymmetryRandom` a length-2 numeric vector, min and max possible
#'  values of the asymmetry coefficient, between 0 and 1.
#'  * `HeightBreadthRandom` a length-2 numeric vector, min and max possible
#'  values for the breadth of the height niche. The minimum and maximum relative
#'  heights for the species are MeanHeight +- (HeightBreadth /2), where MeanHeight
#'  is sampled randomly between 0 and 1. The resulting min and max relative
#'  heights are then used to sample the minimum and maximum of the light niche:
#'   \deqn{I_{max} e^{-k_L LAI (1 - MinHeight)}}
#'  * `Imax` numeric, maximum light intensity used to generate the light niche
#'  * `LAI` numeric, the leaf area index of this species
#'  * `kL`, the light extinction coefficient used to generate the light niche
#'
#'
#' @details
#' Additional details... Age-mass equation, light equation,...
#'
#'
#' @returns a named list of numeric containing the following traits:
#'
#'  * `MaximumMass` positive numeric, the maximum mass an individual can reach.
#'  * `MassAtMaturity` numeric between 0 and and `MaximumMass`,
#'  fraction of `MaximumMass` above which at individual can reproduce.
#'  * `GrowthRate` positive numeric, the fraction of remaining growth an individual
#'  gains in a single generation (i.e, mass increase = growth_rate * (max_mass - mass)),
#'  under optimal light conditions.
#'  * `DispersalKernel` positive numeric, the dispersal kernel.
#'  * `DispersalKernelAsymmetry` numeric between 0 and 1, the dispersal asymmetry.
#'  D_k_A = 0.5 corresponds to symmetric dispersal; with D_k_A = 1 individuals
#'  disperse stricly below themselves; with D_k_A = 0 individuals never disperse
#'  only above themselves or at their height.
#'  }
#'  * `RecruitmentInvestmentRel` numeric between 0 and 1, a coefficient scaling
#'  the mass-dependent fecundity coefficient (see [resolve_repro_dispersal()] )
#'  * `RecruitmentInc` a positive coefficient scaling how relative mass growth (see [resolve_repro_dispersal()] )
#'  increases fecundity.
#'  * `MinLight` positive numeric, minimum light conditions under which this species can survive
#'  * `MaxLight` positive numeric, maximum light conditions under which this species can survive
#'  * `OptimumLight` positive numeric, optimum light conditions under which individuals of this species
#'  grow and reproduce at the maximum rate. It is calculated as the average of `MinLight` and `MaxLight.`
#'  * `LightBreadth` positive numeric, range between `MinLight` and `MaxLight`
#'  * `LightResponseA` first term of the parabolic light-growth response function.
#'  * `LightResponseB` second term of the parabolic light-growth response function.
#'  * `LightResponseC` third term of the parabolic light-growth response function.
#'  * `MinHeightRel` minimum relative height (between 0 and 1) at which the species can survive, used to compute `MinLight`
#'  * `MaxHeightRel` maximum relative height (between 0 and 1) at which the species can survive, used to compute `MaxLight`
#'  * `MeanHeightRel` average of `MinHeightRel` and `MaxHeightRel`
#'  * `HeightBreadth` range between `MinHeightRel` and `MaxHeightRel`
#' }
#' @export
#'
draw_species_traits <- function(species_params) {

  # Unpack parameters
  #check_species_params(species_params)
  list2env(species_params, envir = environment())

  # Draw max size
  if (MaxMassLogScaleRandom) {
    MaxMassLog <- runif(1, min = log10(MaxMassRandom[1]),
                        max = log10(MaxMassRandom[2]))
    MaxMass <- 10^MaxMassLog
  } else {
    MaxMass <- runif(1, min = MaxMassRandom[1], max = MaxMassRandom[2])
  }

  # Mass at maturity is a function of the maximum size
  MassAtMaturity <- MaxMass * runif(1, min = MassAtMaturityRelativeRandom[1],
                                    max = MassAtMaturityRelativeRandom[2])

  AgeAtMaturity <- InterceptAgeMaturity * (MaxMass^ScalingAgeMaturity) *
    runif(1, min = 1 - AgeAtMaturityDevCorr, max = 1 + AgeAtMaturityDevCorr)

  # Growth rate of the Bertalanffy growth curve
  K <- -(log(1) + log(1 - (MassAtMaturity / MaxMass))) / AgeAtMaturity

  # Recruitment
  if (CorrelationMassRecruitment) {
    RecruitmentInvestmentRel <- runif(1,
      RecruitmentInvestmentRelMeanCorr * (1 - RecruitmentInvestmentRelDevCorr),
      RecruitmentInvestmentRelMeanCorr * (1 + RecruitmentInvestmentRelDevCorr)
    )
    RecruitmentInc <- 0
  } else {
    RecruitmentInvestmentRel <- runif(1,
      RecruitmentInvestmentRelMeanRandom[1],
      RecruitmentInvestmentRelMeanRandom[2]
    )
    RecruitmentInc <- runif(1, RecruitmentIncRandom[1], RecruitmentIncRandom[2])
      # Not meaningful if no correlation
  }

   # Dispersal
  DispersalKernel <- runif(1, DispersalKernelRandom[1], DispersalKernelRandom[2])
  DispersalKernelAsymmetry <- runif(1, DispersalKernelAsymmetryRandom[1], DispersalKernelAsymmetryRandom[2])

  # Height niche
  MeanHeight <- runif(1, min = 0, max = 1)  # relative height in relation to canopy height
  HeightBreadthTheoretical <- runif(1, HeightBreadthRandom[1], HeightBreadthRandom[2])
  MinHeight <- max(c(0, MeanHeight - (HeightBreadthTheoretical / 2)))
  MaxHeight <- min(c(1, MeanHeight + (HeightBreadthTheoretical / 2)))
  HeightBreadth <- MaxHeight - MinHeight

  # Light niche
  MinLight <- Imax * exp(-kL * LAI * (1 - MinHeight))
  MaxLight <- Imax * exp(-kL * LAI * (1 - MaxHeight))
  OptimumLight <- (MaxLight + MinLight) / 2
  LightBreadth <- MaxLight - MinLight
  light_resp_params <- get_light_resp_params(MinLight, MaxLight, OptimumLight)

  # Output
  sp_traits <- list(
    MaxMass, MassAtMaturity, K, DispersalKernel, DispersalKernelAsymmetry,
    RecruitmentInvestmentRel, RecruitmentInc, MinLight, MaxLight, OptimumLight,
    LightBreadth, light_resp_params[1], light_resp_params[2], light_resp_params[3],
    MinHeight, MaxHeight, MeanHeight, HeightBreadth
  )
  names(sp_traits) <- species_trait_names()
  return(sp_traits)
}

#' Check species traits hyperparameters meet the requirements
#'
#' @param species_params a list with elements
#'
#' @export
#'
check_species_params <- function(species_params) {

  exptd_params <- c(
    "AgeAtMaturityDevCorr",
    "AgeAtMaturityRandom",
    "CorrelationMassRecruitment",
    "DispersalKernelAsymmetryRandom",
    "DispersalKernelRandom",
    "HeightBreadthRandom",
    "Imax",
    "InterceptAgeMaturity",
    "LAI",
    "MassAtMaturityRelativeRandom",
    "MaxMassLogScaleRandom",
    "MaxMassRandom",
    "MaxMassRangeCorr",
    "RecruitmentIncMaxCorr",
    "RecruitmentIncRandom",
    "RecruitmentInvestmentRelDevCorr",
    "RecruitmentInvestmentRelMeanCorr",
    "RecruitmentInvestmentRelMeanRandom",
    "ScalingAgeMaturity",
    "kL"
    )

  list2env(species_params, envir = environment())

  missing_params <- exptd_params[!exptd_params %in% names(species_params)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("species_params", missing_params))
  }

  if (!is.logical(MaxMassLogScaleRandom)) {
    stop("MaxMassLogScaleRandom must be TRUE or FALSE.")
  }

  if (!is.logical(CorrelationMassRecruitment)) {
    stop("CorrelationMassRecruitment must be TRUE or FALSE.")
  }

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(species_params, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- exptd_params[is_missing_val]
    stop(paste(c(
      "The following elements of species_params are NA, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))

  }

  numeric_params <- exptd_params[!exptd_params %in% c("MaxMassLogScaleRandom", "CorrelationMassRecruitment")]
  is_numeric <- sapply(species_params[numeric_params], function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- numeric_params[!is_numeric]
    stop(paste(c(
      "The following elements of species_params must be numeric:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  # Do all uniform distribution parameters have two elements?
  uniform_params <- c("MaxMassRandom", "MassAtMaturityRelativeRandom",
                      "DispersalKernelRandom", "DispersalKernelAsymmetryRandom",
                      "HeightBreadthRandom")
  if (!CorrelationMassRecruitment) {
    uniform_params <- c(uniform_params,
                        "RecruitmentInvestmentRelMeanRandom",
                        "RecruitmentIncRandom")
  }
  is_param_pair <- sapply(species_params[uniform_params], function(x) length(x) == 2)
  if (any(!is_param_pair)) {
    wrong_params <- uniform_params[!is_param_pair]
    stop(paste(c(
      "The following elements of species_params must have two elements (min and max):",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  # Do the elements satisfy min <= max?
  is_min_max <- sapply(species_params[uniform_params], function(x) x[1] <= x[2])
  if (any(!is_min_max)) {
    wrong_params <- uniform_params[!is_min_max]
    stop(paste(c(
      "The following elements of species_params must satisfy min <= max:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }


  # Check that the following parameters are positive:
  positive_params <- c("DispersalKernelRandom", "HeightBreadthRandom",
                       "AgeAtMaturityDevCorr", "Imax", "LAI", "kL")
  if (CorrelationMassRecruitment) {
    positive_params <- c(positive_params,
                         "RecruitmentInvestmentRelMeanCorr",
                         "RecruitmentInvestmentRelDevCorr")
  } else {
    positive_params <- c(positive_params,
                         "RecruitmentInvestmentRelMeanRandom",
                         "RecruitmentIncRandom")
  }
  is_positive <- sapply(species_params[positive_params], function(x) all(x >= 0))
  if (any(!is_positive)) {
    wrong_params <- positive_params[!is_positive]
    stop(paste(c(
      "The following elements of species_params must be positive or zero:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  # Check that the following parameters are between 0 and 1:
  prop_params <- c("MassAtMaturityRelativeRandom", "DispersalKernelAsymmetryRandom")
  is_prop <- sapply(species_params[prop_params], function(x) all(x >= 0) && all(x <= 1))
  if (any(!is_prop)) {
    wrong_params <- prop_params[!is_prop]
    stop(paste(c(
      "The following elements of species_params must be between 0 and 1:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

}

