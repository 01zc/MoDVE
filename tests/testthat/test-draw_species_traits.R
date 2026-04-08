ref_sp_params <- parse_config(test_path("configs", "config_a2.toml"))

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

  # All parameters are numeric, positive etc.
  species_df <- as.data.frame(c("SpeciesID" = 1, species_traits))
  expect_silent(check_species_df(species_df))

  list2env(species_traits, envir = environment())

  # Light parameters satisfy the parabolic response equation
  light_reponse_params <- get_light_resp_params( # TODO: is this function correct?
    MinLight, MaxLight, OptimumLight
    )
  expect_equal(LightResponseA, light_reponse_params[1])
  expect_equal(LightResponseB, light_reponse_params[2])
  expect_equal(LightResponseC, light_reponse_params[3])
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

