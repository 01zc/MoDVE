ref_sp_params <- parse_config("../config_a2.toml")

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
