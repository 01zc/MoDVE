#' Draw species traits from generating parameters
#'
#' @param species_params a list of trait-generating parameter, which must
#' contain the following elements:
#' \itemize{
#'  \item MaxMassLogScaleRandom logical, if `TRUE` then `MaxMass` is sampled in
#'  a uniform between the log10 of both `MaxMassRandom` values.
#'  \item MaxMassRandom a length-2 numeric vector containing the minimum and
#'  maximum MaxMass
#'  \item MassAtMaturityRelativeRandom a length-2 numeric vector containing the
#'  minimum and maximum bounds (between 0 and 1) in which to sample what fraction
#'  of MaxMass the mass at maturity corresponds to.
#'
#'  \item InterceptAgeMaturity numeric, the intercept of the function generating
#'  the age at maturity as a function of mass at maturity
#'  \item ScalingAgeMaturity numeric, the exponent of the function generating
#'  the age at maturity as a function of mass at maturity
#'  \item AgeAtMaturityDevCorr numeric, the range of possible deviation
#'  coefficients in the function generating age at maturity from mass at maturity.
#'
#'  \item CorrelationMassRecruitment logical, controls the option used to
#'  generate the parameters of the mass-to-reproduction allocation function.
#'   If `TRUE`: \itemize{
#'   \item RecruitmentInvestmentRel is sampled in a uniform distribution
#'   with min and max values RecruitmentInvestmentRelMeanCorr *
#'   (1 +- RecruitmentInvestmentRelDevCorr),
#'   }
#'   \item RecruitmentInc is set to 0.
#'   If `FALSE`: \itemize{
#'   \item RecruitmentInvestmentRel is sampled in a uniform distribution
#'   with min and max values RecruitmentInvestmentRelMeanRandom
#'   \item RecruitmentInc is sampled in a uniform distribution with min and max
#'   values RecruitmentIncRandom.
#'   }
#'  \item RecruitmentInvestmentRelMeanCorr see above
#'  \item RecruitmentInvestmentRelDevCorr see above
#'  \item RecruitmentInvestmentRelMeanRandom see above
#'  \item RecruitmentIncRandom
#'
#'  \item DispersalKernelRandom a length-2 numeric vector, min and max possible
#'  values of the dispersal kernel.
#'  \item DispersalKernelAsymmetryRandom a length-2 numeric vector, min and max possible
#'  values of the asymmetry coefficient.
#'  \item HeightBreadthRandom a length-2 numeric vector, min and max possible
#'  values for the breadth of the height niche. The minimum and maximum relative
#'  heights for the species are MeanHeight +- (HeightBreadth /2), where MeanHeight
#'  is sampled randomly between 0 and 1. The resulting min and max relative
#'  heights are then used to sample the minimum and maximum of the light niche:
#'   Imax x exp(-kL x LAI x (1 - MinHeight))
#'  \item Imax numeric, maximum light intensity used to generate the light niche
#'  \item LAI numeric, the leaf area index of this species
#'  \item kL, the light extinction coefficient used to generate the light niche
#' }
#'
#' @returns a named list of numeric containing the following traits:
#' \itemize{
#'  \item MaximumMass positive numeric, the maximum mass an individual can reach.
#'  \item MassAtMaturity numeric between 0 and and MaximumMass,
#'  fraction of MaximumMass above which at individual can reproduce.
#'  \item GrowthRate positive numeric, the fraction of remaining growth an individual
#'  gains in a single generation (i.e, mass increase = growth_rate * (max_mass - mass)),
#'  under optimal light conditions.
#'  \item DispersalKernel positive numeric, the dispersal kernel.
#'  \item DispersalKernelAsymmetry numeric between 0 and 1, the dispersal asymmetry.
#'  D_k_A = 0.5 corresponds to symmetric dispersal; with D_k_A = 1 individuals
#'  disperse stricly below themselves; with D_k_A = 0 individuals never disperse
#'  only above themselves or at their height.
#'  }
#'  \item RecruitmentInvestmentRel numeric between 0 and 1, a coefficient scaling
#'  the mass-dependent fecundity coefficient (see @resolve_repro_dispersal )
#'  \item RecruitmentInc a positive coefficient scaling how relative mass growth (see @resolve_repro_dispersal )
#'  increases fecundity.
#'  \item MinLight positive numeric, minimum light conditions under which this species can survive
#'  \item MaxLight positive numeric, maximum light conditions under which this species can survive
#'  \item OptimumLight positive numeric, optimum light conditions under which individuals of this species
#'  grow and reproduce at the maximum rate. It is calculated as the average of MinLight and MaxLight.
#'  \item LightBreadth positive numeric, range between MinLight and MaxLight
#'  \item LightResponseA first term of the parabolic light-growth response function.
#'  \item LightResponseB second term of the parabolic light-growth response function.
#'  \item LightResponseC third term of the parabolic light-growth response function.
#'  \item MinHeightRel minimum relative height (between 0 and 1) at which the species can survive, used to compute MinLight
#'  \item MaxHeightRel maximum relative height (between 0 and 1) at which the species can survive, used to compute MaxLight
#'  \item MeanHeightRel average of MinHeightRel and MaxHeightRel
#'  \item HeightBreadth range between MinHeightRel and MaxHeightRel
#'  \item MaxRecruits mass-dependent fecundity coefficient when mass = MaximumMass
#'  \item MaxRecruitsAtMassAtMaturity mass-dependent fecundity coefficient when mass = MassAtMaturity
#'  \item AgeAtmaturity
#' }
#' @export
#'
draw_species_traits <- function(species_params) {

  # Unpack parameters
  check_species_params(species_params)
  list2env(config, envir = environment())

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
    RecruitmentNormalizeAtSize1 <- RecruitmentNormalizeAtSize1Corr  # Factor converting the reproductive biomass to potential recruits
    SlopeRecruitment <- 0  # Slope of the correlation between mass and recruitment
    InterceptRecruitment <- RecruitmentNormalizeAtSize1
    RecruitmentInc <- 0
  } else {
    RecruitmentInvestmentRel <- runif(1,
      RecruitmentInvestmentRelMeanRandom[1],
      RecruitmentInvestmentRelMeanRandom[2]
    )
    RecruitmentNormalizeAtSize1 <- runif(1,
      RecruitmentNormalizeAtSize1Random[1],
      RecruitmentNormalizeAtSize1Random[2]
    )
    SlopeRecruitment <- 0  # No slope if no correlation is choosen
    InterceptRecruitment <- RecruitmentNormalizeAtSize1 - SlopeRecruitment
    RecruitmentInc <- runif(1, RecruitmentIncRandom[1], RecruitmentIncRandom[2])
      # Not meaningful if no correlation
  }
  # TODO: InterceptRecruitment, SlopeRecruitment are not used (pass via config instead) -> DELETE?
  # TODO: MaxRecruits and MaxRecruitsMaturity are not used -> DELETE?
  MaxRecruits <- InterceptRecruitment * RecruitmentInvestmentRel
  MaxRecruitsMaturity <- (InterceptRecruitment + SlopeRecruitment *
                            MassAtMaturity) * RecruitmentInvestmentRel
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
    MinHeight, MaxHeight, MeanHeight, HeightBreadth,
    MaxRecruits, MaxRecruitsMaturity, AgeAtMaturity
  )
  names(sp_traits) <- species_trait_names()
  return(sp_traits)
}
