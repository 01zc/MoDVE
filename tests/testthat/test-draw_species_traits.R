ref_sp_params <- parse_config("../config_a2.toml")

#' Shorthand function for generating hyperparameters
draw_rnd_species_params <- function(max_val = 100) {

  sp_params <- list(
    MaxMassLogScaleRandom = runif(1) > 0.5, # coin flip
    MaxMassRandom = runif(1, 0, max_val), #
    InterceptAgeMaturity = runif(1, 0, max_val),
    AgeAtMaturityDevCorr = runif(1, 0, 1),
    ScalingAgeMaturity = runif(1, -max_val, max_val),
    MassAtMaturityRelativeRandom = runif(1, 0, 1), ##

    CorrelationMassRecruitment = runif(1) > 0.5,

    DispersalKernelRandom = runif(1, 0, max_val), #
    DispersalKernelAsymmetryRandom = runif(1, 0, 1), ##

    HeightBreadthRandom = runif(1, 0, max_val), #

    Imax = runif(1, 0, max_val),
    kL = runif(1, 0, max_val),
    LAI = runif(1, 0, max_val)
  )

  # Draw second value for uniform parameters
  draw_max_value <- function(min_val, upper_bound) {
    return(min_val + runif(1, 0, upper_bound - min_val))
  }
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

test_that("Input is checked correctly", {

  sp_params <- ref_sp_params
  expect_silent(check_species_params(sp_params))

  sp_params$MaxMassLogScaleRandom <- 5
  expect_error(check_species_params(sp_params),
               "MaxMassLogScaleRandom must be TRUE or FALSE.")

  sp_params <- ref_sp_params
  sp_params$DispersalKernelRandom <- "missing"
  sp_params$MaxMassRandom <- "values"
  expect_error(check_species_params(sp_params),
               "The following elements of species_params must be numeric:  DispersalKernelRandom  MaxMassRandom")

  sp_params <- ref_sp_params
  sp_params$AgeAtMaturityDevCorr[2] <- NA
  expect_error(check_species_params(sp_params),
               "The following elements of species_params are NA, NaN, or NULL:  AgeAtMaturityDevCorr")

  sp_params <- ref_sp_params
  sp_params$CorrelationMassRecruitment <- FALSE
  sp_params$RecruitmentIncRandom[1] <- sp_params$RecruitmentIncRandom[2] + 1
  expect_error(check_species_params(sp_params),
               "The following elements of species_params must satisfy min <= max:  RecruitmentIncRandom")

  sp_params <- ref_sp_params
  sp_params$LAI <- -10
  expect_error(check_species_params(sp_params),
               "The following elements of species_params must be positive or zero:  LAI")

})

test_that("Species traits meet requirements", {

  sp_params <- draw_rnd_species_params(max_val = 100)
  expect_silent(species_traits <- draw_species_traits(sp_params))

  # All expected elements are present
  expect_equal(names(species_traits), species_trait_names())

  list2env(species_traits, envir = environment())
  expect_gte(MassAtMaturity, 0)
  expect_gte(MaximumMass, MassAtMaturity)
  expect_gte(GrowthRate, 0)

  expect_gte(DispersalKernel, 0)
  expect_gte(DispersalKernelAsymmetry, 0)
  expect_lte(DispersalKernelAsymmetry, 1)

  expect_gte(RecruitmentInvestmentRel, 0)
  #expect_lte(RecruitmentInvestmentRel, 1) # i believe this can be above 1
  expect_gte(RecruitmentInc, 0)

  expect_gte(MinLight, 0)
  expect_gte(OptimumLight, MinLight)
  expect_gte(MaxLight, OptimumLight)
  expect_equal(LightBreadth, MaxLight - MinLight)

  # Light parameters satisfy the parabolic response equation
  light_reponse_params <- get_light_resp_params( # TODO: is this function correct?
    MinLight, MaxLight, OptimumLight
    )
  expect_equal(LightResponseA, light_reponse_params[1])
  expect_equal(LightResponseB, light_reponse_params[2])
  expect_equal(LightResponseC, light_reponse_params[3])

  expect_gte(MinHeightRel, 0)
  expect_gte(MaxHeightRel, MinHeightRel)
  expect_lte(MaxHeightRel, 1)
  expect_gte(MeanHeightRel, 0)
  expect_lte(MeanHeightRel, 1)
  expect_equal(HeightBreadth, MaxHeightRel - MinHeightRel)

})

test_that("Growth rate satisfies its equation", {

  sp_params <- draw_rnd_species_params(max_val = 100)
  sp_params$AgeAtMaturityDevCorr <- 0 # no variation in age at maturity

  species_traits <- draw_species_traits(sp_params)

  exptd_age_maturity <- sp_params$InterceptAgeMaturity *
    (species_traits$MaximumMass ^ sp_params$ScalingAgeMaturity)

  exptd_mass_mat <- species_traits$MaximumMass *
    (1 - exp(-species_traits$GrowthRate * exptd_age_maturity))

  expect_equal(species_traits$MassAtMaturity, exptd_mass_mat)
})

