source("../old_funcs.R")
source("../test-utils.R")
species_params <- parse_config("../config_a2.toml")

test_that("consistent with old version", {

  nb_species <- 1
  SpeciesPool <- create_rnd_species_df(nb_species, species_params)
  SurfaceBiomassScaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimensions <- sample(2:10, 3)
  centralPoint <- find_central_point(dimensions)

  Microhabitat <- create_rnd_microhabitat(SpeciesPool, dimensions, SurfaceBiomassScaling)

  init_params <- draw_rnd_initial_inds_params()
  init_params$PercentageMaturePerSpecies <- rep(100, nb_species)
  E <- draw_initial_individuals(
    init_params,
    SpeciesPool,
    Microhabitat
  )
  max_id <- max(E$IndividualID)

  InterceptRecruitment <- runif(1, 0, 100)
  SlopeRecruitment <- runif(1, 0, 1)

  # Format input to the old version
  dims_with_corr <- dimPlot <-  dimensions * 2 + 1
  ProbabilityMatrixNormalized <- old_compute_prob_matrix_norm(
    centralPoint,
    dims_with_corr[1],
    dims_with_corr[2],
    dims_with_corr[3],
    nb_species,
    SpeciesPool
  )
  E[, (ncol(E)+1):(ncol(E)+1+ncol(SpeciesPool))] <- SpeciesPool
  {
    NumberOfSpecies <- nb_species
    dimPlot <- dims_with_corr
    MaxIndividualID <- max_id
  }
  disp_list <- old_dispersal(
    nb_species,
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    dimensions,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    ProbabilityMatrixNormalized,
    SpeciesPool,
    max_id
  )

  recruitment_df <- disp_list$PotentialRecruitment
  names(recruitment_df) <- c("species_index", "exptd_nb_recruits")
  recruitment_df$species_index[recruitment_df$species_index == 0] <- 1
  recruitment_df$nb_recruits <- as.numeric(disp_list$NumberRecruitsPerSpecies)

  res_exptd <- list(
    "nbIndsBeforeDisp" = disp_list$IntialNumberIndividuals,
    "E" = disp_list$E,
    "recruitment_df" = recruitment_df,
    "max_id" = disp_list$MaxIndividualID
  )

  prob_disp_matrix <- calc_prob_disp_matrix(centralPoint,
                                            dimensions[1],
                                            dimensions[2],
                                            dimensions[3],
                                            SpeciesPool)
  prob_disp_matrix <- ProbabilityMatrixNormalized

  res_obs <- resolve_repro_dispersal(
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )

  expect_equal(res_obs, res_exptd)
})
