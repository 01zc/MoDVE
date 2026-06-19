test_that("dispersal consistent with the previous version", {

  # Generate species
  NumberOfSpecies <- 1
  SpeciesPool <- create_rnd_species_df(NumberOfSpecies)

  # Global parameters
  InterceptRecruitment <- runif(1, 0, 100)
  SlopeRecruitment <- runif(1, 0, 1)
  SurfaceBiomassScaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimPlot <- sample(2:10, 3)

  # Initialise microhabitat
  Microhabitat <- create_rnd_microhabitat(SpeciesPool, dimPlot, SurfaceBiomassScaling)

  # Initialise individuals table
  init_params <- draw_rnd_initial_inds_params()
  init_params$SurfaceBiomassScaling <- SurfaceBiomassScaling # use same as above
  init_params$PercentageMaturePerSpecies <- rep(100, NumberOfSpecies) # only adults
  E <- draw_initial_individuals(
    init_params,
    SpeciesPool,
    Microhabitat
  ) |>
    dplyr::filter(Status == 1) # exclude unplaced dead individuals
  E$IndividualID <- 1:nrow(E)
  MaxIndividualID <- max(E$IndividualID)

  # Compute the dispersal matrix
  ProbabilityMatrixNormalized <- calc_prob_disp_matrix(dimPlot, SpeciesPool)

  # 1 - Generate expectation with old version
  # Format input to the old format
  expanded_dims <- dimPlot * 2 + 1
  centralPoint <- floor(expanded_dims / 2) + 1
  # Collate species traits to individual table
  extra_col_indices <- (ncol(E) + 1):(ncol(E) + ncol(SpeciesPool) - 1)
  E[, extra_col_indices] <- SpeciesPool[, 2:ncol(SpeciesPool)]
  rng_state <- .Random.seed # use same RNG for both versions
  disp_list <- old_dispersal(
    NumberOfSpecies,
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    dimPlot,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    ProbabilityMatrixNormalized,
    SpeciesPool,
    MaxIndividualID
  )
  # Format output
  recruitment_df <- disp_list$PotentialRecruitment
  names(recruitment_df) <- c("species_index", "exptd_nb_recruits")
  recruitment_df$species_index[recruitment_df$species_index == 0] <- 1
  recruitment_df$nb_recruits <- as.numeric(disp_list$NumberRecruitsPerSpecies)
  res_exptd <- list(
    "nbIndsBeforeDisp" = disp_list$IntialNumberIndividuals,
    "E" = disp_list$E[, 1:9], # we don't care about species-level data
    "recruitment_df" = recruitment_df,
    "max_id" = disp_list$MaxIndividualID
  )

  # 2 - Run the current algorithm
  # restore seed to ensure both versions use same RNG
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  res_obs <- resolve_repro_dispersal(
    E[,1:9],
    Microhabitat,
    SurfaceBiomassScaling,
    InterceptRecruitment,
    SlopeRecruitment,
    ProbabilityMatrixNormalized,
    SpeciesPool,
    MaxIndividualID
  )
  expect_equal(res_obs, res_exptd)
})

test_that("No dispersal in non-suitable voxels", {

  # Generate species
  NumberOfSpecies <- 1
  SpeciesPool <- create_rnd_species_df(NumberOfSpecies)

  # Set fecundity so not too many newborns
  fecundity <- 10 # seedlings per plant
  mass <- 500 # grams

  # Global parameters
  InterceptRecruitment <- 0
  SlopeRecruitment <- fecundity / mass
  SurfaceBiomassScaling <- 10
  avail_sa <- mass_to_surf_area(mass, SurfaceBiomassScaling) * 2

  SpeciesPool$RecruitmentInvestmentRel <- 1
  SpeciesPool$RecruitmentInc <- 0

  dimPlot <- rep(3, 3)

  # Initialise microhabitat with all microclimatic variables
  Microhabitat <- create_empty_microhabitat(
    dimPlot,
    microclimate_vars = c("wind", "temperature", "humidity")
  )
  # For each variable, assign two random voxels to be unsuitable
  # Doesn't matter if they overlap
  nb_voxels <- prod(dimPlot)

  Microhabitat[,,,1] <- avail_sa # surface area not limiting
  layer_light <- rep(SpeciesPool$OptimumLight, nb_voxels)
  layer_hum <- rep(SpeciesPool$OptimumHum, nb_voxels)
  layer_temp <- rep(SpeciesPool$OptimumTemp, nb_voxels)
  layer_wind <- rep(SpeciesPool$OptimumWind, nb_voxels)
  unsuitable_light <- sample(1:nb_voxels, 2)
  unsuitable_humidity <- sample(1:nb_voxels, 2)
  unsuitable_wind <- sample(1:nb_voxels, 2)
  unsuitable_temperature <- sample(1:nb_voxels, 2)
  layer_light[unsuitable_light] <- SpeciesPool$MinLight / 2
  layer_hum[unsuitable_humidity] <- SpeciesPool$MinHum / 2
  layer_temp[unsuitable_temperature] <- SpeciesPool$MinTemp / 2
  layer_wind[unsuitable_wind] <- SpeciesPool$MinWind / 2
  Microhabitat[,,,3] <- layer_light
  Microhabitat[,,,4] <- layer_wind
  Microhabitat[,,,5] <- layer_temp
  Microhabitat[,,,6] <- layer_hum

  # Initialise individuals table
  init_params <- draw_rnd_initial_inds_params()
  init_params$SurfaceBiomassScaling <- SurfaceBiomassScaling # use same as above
  init_params$PercentageMaturePerSpecies <- rep(100, NumberOfSpecies) # only adults
  E <- draw_initial_individuals(
    init_params,
    SpeciesPool,
    Microhabitat
  ) |>
    dplyr::filter(Status == 1) # exclude unplaced dead individuals
  E$IndividualID <- 1:nrow(E)
  max_id <- max(E$IndividualID)
  E$Mass <- mass

  # Compute the dispersal matrix
  # No wind-mediated dispersal here
  prob_disp_matrix <- calc_prob_disp_matrix(dimPlot, SpeciesPool)

  # restore seed to ensure both versions use same RNG
  E <- resolve_repro_dispersal(
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )$E
  E <- E |> dplyr::mutate(
    "cell_index" = index_3d(X, Y, Z, dimPlot[1], dimPlot[2]),
    "in_unsuitable_cell" = cell_index %in% c(
      unsuitable_light,
      unsuitable_humidity,
      unsuitable_temperature,
      unsuitable_wind
    )
  )
  expect_true(!any(E$in_unsuitable_cell))
})

