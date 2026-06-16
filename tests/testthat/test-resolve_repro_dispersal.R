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

  species_df <- SpeciesPool
  microhab_mat <- Microhabitat
  distr_params <- init_params

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
