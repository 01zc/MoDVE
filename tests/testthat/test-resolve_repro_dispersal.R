source("tests/old_funcs.R")
source("tests/test-utils.R")

species_params <- parse_config("tests/config_a2.toml")

test_that("consistent with old version", {

  nb_species <- 1
  SpeciesPool <- create_rnd_species_df(nb_species, species_params)
  surface_biomass_scaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimensions <- sample(1:10, 3)
  centralPoint <- find_central_point(dimensions)

  Microhabitat <- array(0, c(dimensions, 3))
  nb_suitable_voxels <- round(prod(dimensions) * runif(1, 0, 1))
  suitable_voxels <- sample(1:prod(dimensions), nb_suitable_voxels)
  reqd_sa_per_ind <- SpeciesPool$MaximumMass^2/3 / surface_biomass_scaling
  surface_area_mat <- array(0, dim = dimensions)
  surface_area_mat[suitable_voxels] <- reqd_sa_per_ind *
    sample(1:10, length(suitable_voxels), replace = TRUE)
  Microhabitat[, , , 1] <- surface_area_mat
  Microhabitat[, , , 3] <- SpeciesPool$OptimumLight

  E <- draw_initial_individuals(
    draw_rnd_initial_inds_params(),
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
  #  expected_sa = Mass^(2 / 3) / surface_biomass_scaling
  #)

  E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
  E$MassAtMaturity
  min(E$LightInVoxel[is_sp])

  prob_disp_matrix <- calc_prob_disp_matrix(centralPoint,
                                            dimensions[1],
                                            dimensions[2],
                                            dimensions[3],
                                            SpeciesPool)

  max_id <- max(E$IndividualID)

  disp_list <- old_dispersal(
    nb_species,
    E,
    Microhabitat,
    surface_biomass_scaling,
    dimensions,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )

  recruitment_df <- disp_list$PotentialRecruitment
  names(recruitment_df) <- c("species_index", "exptd_nb_recruits")
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
    surface_biomass_scaling,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )

  expect_equal(res_obs, res_exptd)

})
