#' Draw species traits from generating hyperparameters
#'
#' Draw random values for the demographic and niche traits that define an
#' epiphyte species in the simulation.
#'
#' @param species_params a list of trait-generating parameter, which must
#' contain the following elements:
#'
#'  * `microclimate_opts` a list specifying which climatic variables are used,
#'  with the following named logical elements: `use_temperature`, `use_wind`,
#'  and `use_humidity`.
#'  * `MaxMassLogScaleRandom` logical. If `FALSE`, `MaxMass` is sampled in a
#'  uniform distribution with parameters `MaxMassRandom`. If `TRUE`, the value
#'  is instead sampled in a uniform of the log10's of `MaxMassRandom`, and the
#'  sampled value is then transformed back. This produces the same range of
#'  values, but skews the distribution towards the minimum.
#'  * `MaxMassRandom` a length-2 numeric vector containing the minimum and
#'  maximum MaxMass
#'  * `MassAtMaturityRelativeRandom` a length-2 numeric vector containing the
#'  minimum and maximum bounds (between 0 and 1) in which to sample what
#'  fraction of `MaxMass` the mass at maturity corresponds to.
#'
#'  * `InterceptAgeMaturity` numeric, the intercept of the function generating
#'  the age at maturity as a function of mass at maturity (see details)
#'  * `ScalingAgeMaturity` numeric between 0 and 1, the exponent of the function
#' generating the age at maturity as a function of mass at maturity
#'  * `AgeAtMaturityDevCorr` numeric, the range of possible deviation
#'  coefficients in the function generating age at maturity from mass at
#'  maturity.
#'
#'  * `CorrelationMassRecruitment` logical, controls the option used to
#'  generate the parameters of the mass-to-reproduction allocation function.
#'   If `TRUE`:
#'   * `RecruitmentInvestmentRel` is sampled in a uniform distribution
#'   with min and max values `RecruitmentInvestmentRelMeanCorr` *
#'   (1 \eqn{\pm} `RecruitmentInvestmentRelDevCorr`),
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
#'   * Either `LightBreadthRandom`, a length-2 numeric vector, giving the min and
#'   max possible  values of the light niche, or equivalently,
#'   `HeightBreadthRandom`, min and max possible values for the height niche.
#'   In the latter case, height niche parameters are converted into light niche parameters by the
#'   following relation:
#'   \deqn{I_{max} * e^{-k_L \times LAI (1 - Height)}}
#'  * `HumBreadthRandom` a length-2 numeric vector, min and max possible
#'  values of the humidity niche.
#'  * `TempBreadthRandom` a length-2 numeric vector, min and max possible
#'  values of the temperature niche
#'  * `WindBreadthRandom` a length-2 numeric vector, min and max possible
#'  values of the wind niche.
#'  * `Imax` numeric, maximum light intensity used to generate the light niche
#'  * `LAI` numeric, the leaf area index of this species
#'  * `kL`, the light extinction coefficient used to generate the light niche
#'
#' @details
#' * Mass at maturity is drawn as a random fraction of the maximum mass.
#' * Age at maturity is sampled as `InterceptAgeMaturity` times
#' `MaxMass^ScalingAgeMaturity` times a random deviation sampled in a uniform
#' with parameters `1 - AgeAtMaturityDevCorr`, `1 + AgeAtMaturityDevCorr`.
#' * The growth rate (\eqn{K}) is derived after the Bertalanffy growth curve,
#' \deqn{M = M_{max} (1 - e^{-K*Age})}, which we resolve for \eqn{K} at the
#' maturity age and mass.
#'
#'
#' * The center of the height niche (relative to canopy height) is
#' sampled between 0 and 1. The breadth of this height niche is sampled between
#' the values of `HeightBreadthRandom`.
#' * The light niche is derived from the height niche, according to
#' \deqn{I_{max} \times e^{k_L \times LAI (1 - height)}} which is resolved for `MinHeight`
#' and `MaxHeight` to obtain `MinLight` and `MaxLight`, with `OptimumLight`
#' being their mean. These values are use to define the parabolic light
#' response with parameters `LightResponseA`, `LightResponseB`, `LightResponseC`
#' such that `MinLight` and `MaxLight` correspond to 0 and `OptimumLight` to 1.
#'
#' * For the humidity, temperature and wind niches, the optimum value is first
#' sampled in U(min+1, max-1), the minimum value in U(min, opt), and the maximum
#' in U(opt, max). This ensures that `MinTemp < OptimumTemp < MaxTemp`
#' (and equivalently for other climatic traits).
#'
#' * The light niche must be symmetric. For this reason, `OptimumLight` is first
#' sampled in U(min+1, max-1) (where min and max are the elements of `LightBreadthRandom`).
#' Then a light niche span is sampled in U(1, min(opt-min, max-opt+10)).
#' The sampled value is respectively subtracted and added to `OptimumLight` to
#' obtain parameters `MinLight` and `MaxLight`. Finally, the three parameters are
#' used to calibrate parameters A, B, C of a parabolic response such that
#' `MinLight` and `MaxLight` correspond to 0 and `OptimumLight` to 1.
#'
#' @returns a named list of numeric containing the following traits:
#'  * `MaximumMass` positive numeric, the maximum mass an individual can reach.
#'  * `MassAtMaturity` numeric between 0 and and `MaximumMass`,
#'  fraction of `MaximumMass` above which at individual can reproduce.
#'  * `GrowthRate` numeric between 0 an 1, the fraction of remaining growth an
#'  individual gains in a single year.
#'  (i.e, \eqn{\Delta M = K \times (M_{max} - M)} under optimal light conditions.
#'  * `DispersalKernel` positive numeric, the dispersal kernel.
#'  * `DispersalKernelAsymmetry` numeric between 0 and 1, the dispersal asymmetry.
#'  \eqn{D_{k_A} = 0.5} corresponds to symmetric dispersal; with
#'  \eqn{D_{k_A} = 1} individuals disperse stricly below themselves; with
#'  \eqn{D_{k_A} = 0} individuals never disperse only above themselves or at
#'  their height.
#'  * `DispersalKernelWindEffect` numeric between 0 and 1, a factor scaling how
#'  much wind affects dispersal, if applicable. If used, wind modifies the
#'  dispersal kernel by
#'  \eqn{e^{-\frac{Dist_V * D_K}{1 + D_W + W}}}
#'  instead of the usual
#'  \eqn{e^{-Dist_V * D_K}}
#'  * `RecruitmentInvestmentRel` numeric between 0 and 1, a coefficient scaling
#'  the mass-dependent fecundity coefficient (see [resolve_repro_dispersal()]),
#'  representing the fraction of available biomass invested in fecundity
#'  * `RecruitmentInc` numeric between 0 and 1, scaling how fecundity increases
#'  with the growth stage of the epiphyte. This factor ranges from `1` when
#'  `Mass` = `MassAtMaturity`, and `2 * RecruitmentInc` when
#'  `Mass` = `MaximumMass`.
#'  * `MinLight` positive numeric, minimum light conditions under which this
#'  species can survive
#'  * `MaxLight` positive numeric, maximum light conditions under which this
#'  species can survive
#'  * `OptimumLight` positive numeric, optimum light conditions under which
#'  individuals of this species
#'  grow and reproduce at the maximum rate.
#'  * `LightResponseA` first term of the parabolic light-growth response function.
#'  * `LightResponseB` second term of the parabolic light-growth response function.
#'  * `LightResponseC` third term of the parabolic light-growth response function.
#'  * `MinTemp` positive numeric, minimum temperature conditions under which this
#'  species can survive
#'  * `MaxTemp` positive numeric, maximum temperature conditions under which this
#'  species can survive
#'  * `OptimumTemp` positive numeric, optimum temperature conditions under which
#'  individuals of this species
#'  grow and reproduce at the maximum rate.
#'  * `MinHum` positive numeric, minimum humidity conditions under which this
#'  species can survive
#'  * `MaxHum` positive numeric, maximum humidity conditions under which this
#'  species can survive
#'  * `OptimumHum` positive numeric, optimum humidity conditions under which
#'  individuals of this species
#'  grow and reproduce at the maximum rate.
#'  * `MinWind` positive numeric, minimum wind conditions under which this
#'  species can survive
#'  * `MaxWind` positive numeric, maximum wind conditions under which this
#'  species can survive
#'  * `OptimumWind` positive numeric, optimum wind conditions under which
#'  individuals of this species
#'  grow and reproduce at the maximum rate.
#'
#' @export
#'
draw_species_traits <- function(species_params) {

  check_species_params(species_params)

  sp <- species_params # shorter alias

  # Draw max size
  if (sp$MaxMassLogScaleRandom) {
    if (sp$MaxMassRandom[1] == 0) sp$MaxMassRandom[1] <- 1e-9 # otherwise NaN
    MaxMassLog <- stats::runif(1, min = log10(sp$MaxMassRandom[1]),
                        max = log10(sp$MaxMassRandom[2]))
    MaxMass <- 10^MaxMassLog
  } else {
    MaxMass <- stats::runif(1, min = sp$MaxMassRandom[1], max = sp$MaxMassRandom[2])
  }

  # Mass at maturity is a function of the maximum size
  MassAtMaturity <- MaxMass * stats::runif(1, min = sp$MassAtMaturityRelativeRandom[1],
                                    max = sp$MassAtMaturityRelativeRandom[2])

  AgeAtMaturity <- sp$InterceptAgeMaturity * (MaxMass^sp$ScalingAgeMaturity) *
    stats::runif(1, min = 1 - sp$AgeAtMaturityDevCorr, max = 1 + sp$AgeAtMaturityDevCorr)

  # Growth rate of the Bertalanffy growth curve
  # Mass = MaxMass * (1 - exp(-K*Age)) --> resolved for K
  K <- min(1, -log(1 - (MassAtMaturity / MaxMass)) / AgeAtMaturity)
  # Age cannot be more than 1, otherwise risk of overshooting max mass

  # Recruitment
  if (sp$CorrelationMassRecruitment) {
    RecruitmentInvestmentRel <- stats::runif(1,
      sp$RecruitmentInvestmentRelMeanCorr * (1 - sp$RecruitmentInvestmentRelDevCorr),
      sp$RecruitmentInvestmentRelMeanCorr * (1 + sp$RecruitmentInvestmentRelDevCorr)
    )
    RecruitmentInc <- 0
  } else {
    RecruitmentInvestmentRel <- stats::runif(1,
      sp$RecruitmentInvestmentRelMeanRandom[1],
      sp$RecruitmentInvestmentRelMeanRandom[2]
    )
    RecruitmentInc <- stats::runif(1, sp$RecruitmentIncRandom[1],
                                   sp$RecruitmentIncRandom[2])
      # Not meaningful if no correlation
  }
  # Both parameters must be between 0 and 1
  RecruitmentInvestmentRel <- min(1, max(0, RecruitmentInvestmentRel))
  RecruitmentInc <- min(1, max(0, RecruitmentInc))

  # Dispersal
  DispersalKernel <- stats::runif(1, sp$DispersalKernelRandom[1],
                                  sp$DispersalKernelRandom[2])
  DispersalKernelAsymmetry <- stats::runif(1,
                                    sp$DispersalKernelAsymmetryRandom[1],
                                    sp$DispersalKernelAsymmetryRandom[2])
  if (sp$microclimate_opts$use_wind) {
    DispersalKernelWindEffect <- runif(1, min = 0, max = 1)
  }

  # Niche values
  if (is.null(sp$LightBreadthRandom)) {
    MinLightRandom <- sp$Imax * exp(-sp$kL * sp$LAI * (1 - sp$HeightBreadthRandom[1]))
    MaxLightRandom <- sp$Imax * exp(-sp$kL * sp$LAI * (1 - sp$HeightBreadthRandom[2]))
    sp$LightBreadthRandom <- c(MinLightRandom, MaxLightRandom)
  }

  OptimumLight <- runif(1, min = sp$LightBreadthRandom[1] + 1,
                        max = sp$LightBreadthRandom[2] - 1)

  # Light niche must be symmetric
  max_light_niche_span <- min(
    OptimumLight - sp$LightBreadthRandom[1],
    sp$LightBreadthRandom[2] - OptimumLight + 10
  ) # Ensure that light breadth does not exceed the range of light values
  light_niche_span <- runif(1, min = 1, max = max_light_niche_span)
  MinLight <- OptimumLight - light_niche_span
  MaxLight <- OptimumLight + light_niche_span

  light_resp_params <- get_light_resp_params(MinLight, MaxLight, OptimumLight)

  if (sp$microclimate_opts$use_wind) {
    WindMargin <- runif(1, 0.001, 1)
    OptimumWind <- runif(1, min = sp$WindBreadthRandom[1] + WindMargin,
                         max = sp$WindBreadthRandom[2] - WindMargin)
    MinWind <- runif(1, min = sp$WindBreadthRandom[1], max = OptimumWind)
    MaxWind <- runif(1, min = OptimumWind, max = sp$WindBreadthRandom[2])
  }

  if (sp$microclimate_opts$use_temperature) {
    OptimumTemp <- runif(1, min = sp$TempBreadthRandom[1] + 1,
                         max = sp$TempBreadthRandom[2] - 1)
    MinTemp <- runif(1, min = sp$TempBreadthRandom[1], max = OptimumTemp)
    MaxTemp <- runif(1, min = OptimumTemp, max = sp$TempBreadthRandom[2])
  }

  if (sp$microclimate_opts$use_humidity) {
    OptimumHum <- runif(1, min = sp$HumBreadthRandom[1] + 1,
                        max = sp$HumBreadthRandom[2] - 1)

    MinHum <- runif(1, min = sp$HumBreadthRandom[1], max = OptimumHum)

    MaxHum <- runif(1, min = OptimumHum, max = sp$HumBreadthRandom[2])
  }

  # Output
  sp_traits <- list(
    MaxMass, MassAtMaturity, K, DispersalKernel, DispersalKernelAsymmetry,
    RecruitmentInvestmentRel, RecruitmentInc,
    MinLight, MaxLight, OptimumLight,
    light_resp_params[1], light_resp_params[2], light_resp_params[3]
  )
  if (sp$microclimate_opts$use_humidity) {
    sp_traits <- append(sp_traits, MinHum, MaxHum, OptimumHum)
  }
  if (sp$microclimate_opts$use_temperature) {
    sp_traits <- append(sp_traits, MinTemp, MaxTemp, OptimumTemp,)

  }
  if (sp$microclimate_opts$use_wind) {
    sp_traits <- append(sp_traits, MinWind, MaxWind, OptimumWind,
                        DispersalKernelWindEffect)
  }
  names(sp_traits) <- species_trait_names(sp$microclimate_opts)
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
    "CorrelationMassRecruitment",
    "DispersalKernelAsymmetryRandom",
    "DispersalKernelRandom",
    "Imax",
    "InterceptAgeMaturity",
    "LAI",
    "MassAtMaturityRelativeRandom",
    "MaxMassLogScaleRandom",
    "MaxMassRandom",
    "ScalingAgeMaturity",
    "kL"
    )

  if (is.null(species_params$microclimate_opts)) {
    stop("species_params must contain list microclimate_opts with logical elements use_wind, use_temperature and use_humidity.")
  }

  if (species_params$microclimate_opts$use_temperature) {
    exptd_params <- c(exptd_params, "TempBreadthRandom")
  }

  if (species_params$microclimate_opts$use_humidity) {
    exptd_params <- c(exptd_params, "HumBreadthRandom")
  }

  if (species_params$microclimate_opts$use_wind) {
    exptd_params <- c(exptd_params, "WindBreadthRandom")
  }

  if (species_params$CorrelationMassRecruitment) {
    exptd_params <- c(exptd_params,
                      "RecruitmentInvestmentRelMeanCorr",
                      "RecruitmentInvestmentRelDevCorr")
  } else {
    exptd_params <- c(exptd_params,
                         "RecruitmentInvestmentRelMeanRandom",
                         "RecruitmentIncRandom")
  }


  missing_params <- exptd_params[!exptd_params %in% names(species_params)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("species_params", missing_params))
  }

  if (is.null(species_params$LightBreadthRandom) &&
      is.null(species_params$HeightBreadthRandom)) {
    stop("species_params must contain either LightBreadthRandom or HeightBreadthRandom.")
  }


  if (!is.logical(species_params$MaxMassLogScaleRandom)) {
    stop("MaxMassLogScaleRandom must be TRUE or FALSE.")
  }

  if (!is.logical(species_params$CorrelationMassRecruitment)) {
    stop("CorrelationMassRecruitment must be TRUE or FALSE.")
  }

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(species_params[exptd_params], function(x)  {
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
                      "DispersalKernelRandom", "DispersalKernelAsymmetryRandom")
  if (!species_params$CorrelationMassRecruitment) {
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
  positive_params <- c("InterceptAgeMaturity", "DispersalKernelRandom", "HeightBreadthRandom",
                       "Imax", "LAI", "kL")
  if (species_params$CorrelationMassRecruitment) {
    positive_params <- c(positive_params, "RecruitmentInvestmentRelMeanCorr")
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
  prop_params <- c(
    "MassAtMaturityRelativeRandom", "DispersalKernelAsymmetryRandom",
    "AgeAtMaturityDevCorr", "ScalingAgeMaturity")
  if (species_params$CorrelationMassRecruitment) {
    prop_params <- c(prop_params, "RecruitmentInvestmentRelDevCorr")
  }
  is_prop <- sapply(species_params[prop_params], function(x) all(x >= 0) && all(x <= 1))
  if (any(!is_prop)) {
    wrong_params <- prop_params[!is_prop]
    stop(paste(c(
      "The following elements of species_params must be between 0 and 1:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

}

