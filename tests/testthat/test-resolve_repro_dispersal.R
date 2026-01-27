source("tests/old_funcs.R")
source("tests/test-utils.R")

species_params <- parse_config("tests/config_a2.toml")

test_that("consistent with old version", {

  nb_species <- 1
  SpeciesPool <- create_rnd_species_df(nb_species, species_params)
  SurfaceBiomassScaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimensions <- sample(2:10, 3)
  centralPoint <- find_central_point(dimensions)

  Microhabitat <- array(0, c(dimensions, 3))
  nb_suitable_voxels <- round(prod(dimensions) * runif(1, 0, 1))
  suitable_voxels <- sample(1:prod(dimensions), nb_suitable_voxels)
  reqd_sa_per_ind <- SpeciesPool$MaximumMass^(2/3) / SurfaceBiomassScaling
  surface_area_mat <- array(0, dim = dimensions)
  surface_area_mat[suitable_voxels] <- reqd_sa_per_ind *
    sample(1:10, length(suitable_voxels), replace = TRUE)
  Microhabitat[, , , 1] <- surface_area_mat
  Microhabitat[, , , 3] <- SpeciesPool$OptimumLight

  init_params <-draw_rnd_initial_inds_params()
  init_params$PercentageMaturePerSpecies <- rep(100, nb_species)
  E <- draw_initial_individuals(
    init_params,
    SpeciesPool,
    Microhabitat
  )

  # nb_inds <- sample(1:dimensions[1], 1)
  # nb_inds <- 20
  # E <- tibble::tibble(
  #  SpeciesID = 1,
  #  IndividualID = seq_len(nb_inds),
  #  Status = 1,
  #  # Initialise individuals randomly
  #  X = rep(sample(dimensions[1], nb_inds, replace = TRUE)),
  #  Y = rep(sample(dimensions[2], nb_inds, replace = TRUE)),
  #  Z = rep(sample(dimensions[3], nb_inds, replace = TRUE)),
  #  Mass = runif(nb_inds, 0, 1),
  #  MassAtMaturity = Mass,
  #  MaxMass = Mass,
  #  RecruitmentInvestmentRel = runif(nb_inds, 0, 1),
  #  # Compute expected surface area
  #  expected_sa = Mass^(2 / 3) / SurfaceBiomassScaling
  #)

  # E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0

  prob_disp_matrix <- calc_prob_disp_matrix(centralPoint,
                                            dimensions[1],
                                            dimensions[2],
                                            dimensions[3],
                                            SpeciesPool)

  max_id <- max(E$IndividualID)
  InterceptRecruitment <- runif(0, 100, 1)
  SlopeRecruitment <- runif(0, 1, 1)

  # Format input to the old version
  dimX <- dimensions[1] * 2 + 1
  dimY <- dimensions[2] * 2 + 1
  dimZ <- dimensions[3] * 2 + 1
  central_point_old <- c(
    floor(dimX/2) + 1,
    floor(dimY/2) + 1,
    floor(dimZ/2) + 1
  )
  prob_disp_matrix_old <- old_compute_prob_matrix_norm(
    central_point_old, dimX, dimY, dimZ, nb_species, SpeciesPool
    )
  E_old <- E
  E_old[, (ncol(E)+1):(ncol(E)+1+ncol(SpeciesPool))] <- SpeciesPool


  disp_list <- old_dispersal(
    nb_species,
    E_old,
    Microhabitat,
    SurfaceBiomassScaling,
    dimensions,
    central_point_old,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
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
