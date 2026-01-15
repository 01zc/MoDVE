source("tests/test-utils.R")

test_that("3D coordinates are converted to sequential index correctly", {

  rnd_dim <- sample(1:10, 3, replace = TRUE)
  seq_mat <- array(1:prod(rnd_dim), dim = rnd_dim)
  coords_tbl <- tidyr::expand_grid(
    x = 1:(rnd_dim[1]),
    y = 1:(rnd_dim[2]),
    z = 1:(rnd_dim[3])
  ) |>
    dplyr::mutate(
      i = index_3d(x, y, z, rnd_dim[1], rnd_dim[2])
      )
  for (x in 1:rnd_dim[1]) {
    for (y in 1:rnd_dim[2]) {
      for (z in 1:rnd_dim[3]) {
        obs_i <- coords_tbl$i[coords_tbl$x == x & coords_tbl$y == y &  coords_tbl$z == z]
        exptd_i <- seq_mat[x, y, z]
        expect_equal(obs_i, exptd_i)
      }
    }
  }
})

test_that("Initial individuals are distributed correctly", {

  rnd_params <- list(
    "SurfaceBiomassScaling" = runif(1, 1e-7, 100),
    "IndividualsPerSpecies" = sample(1:100, 1),
    "ScalingPerHa" = FALSE,
    "PercentageMaturePerSpecies" = runif(1, 0, 100)
  )
  nb_species <- sample(1:5, 1)
  total_nb_inds <- nb_species * rnd_params$IndividualsPerSpecies
  species_df <- create_rnd_species_df(nb_species, draw_rnd_species_params())

  # Create a random grid of voxels
  dimensions <- sample(2:10, 3, replace = TRUE)
  nb_voxels <- sample(1:prod(dimensions), 1)
  unlimited_sa <- 100000
  microhab_mat <- array(unlimited_sa, dim = c(dimensions, 3))

  # light not limiting for now
  species_df$MinLight <- 0
  species_df$MaxLight <- 100
  species_df$LightBreadth <- 100
  microhab_mat[,,,3] <- 50 # optimal light conditions

  {
    distr_params <- rnd_params
    largest_inds_first = FALSE
    most_surf_area_first = FALSE
  }

  # Edge case - no mature individuals
  rnd_percent <- rnd_params$PercentageMaturePerSpecies
  rnd_params$PercentageMaturePerSpecies <- 0
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  init_ind_df <- init_ind_df |> dplyr::mutate(
    "is_mature" = Mass >= species_df$MassAtMaturity[SpeciesID],
  )
  expect_equal(sum(init_ind_df$is_mature), 0)
  # Edge case - all mature individuals
  rnd_params$PercentageMaturePerSpecies <- 100
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  init_ind_df <- init_ind_df |> dplyr::mutate(
    "is_mature" = Mass >= species_df$MassAtMaturity[SpeciesID],
  )
  expect_true(all(init_ind_df$is_mature))
  rnd_params$PercentageMaturePerSpecies <- rnd_percent

  # Edge case: 1D matrix
  flat_dim <- c(1, 1, 10)
  microhab_1d <- array(unlimited_sa, dim = c(flat_dim, 3))
  microhab_1d[,,,3] <- 50 # optimal light conditions
  expect_silent(draw_initial_individuals(rnd_params, species_df, microhab_1d))
  # 2D matrix
  square_matrix <- array(unlimited_sa, dim = c(c(2, 1, 3), 3))
  square_matrix[,,,3] <- 50 # optimal light conditions
  expect_silent(draw_initial_individuals(rnd_params, species_df, square_matrix))

  # Case 1 - empty matrix, no surface area
  microhab_mat[,,,1] <- 0
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_equal(nrow(init_ind_df), total_nb_inds)
  # default coordinates is 1,1,1
  expect_true(all(init_ind_df$X == 1 & init_ind_df$Y == 1 & init_ind_df$Z == 1))
  expect_true(all(init_ind_df$Status == 2)) # dead
  # Output format is correct
  exptd_cols <- inds_input_names()
  missing_cols <- exptd_cols[!exptd_cols %in% names(init_ind_df)]
  expect_length(missing_cols, 0)

  # Case 2 - Surface area is not limiting
  surf_area_mat <- rep(0, prod(dimensions))
  suitable_voxels <- sample(1:prod(dimensions), size = nb_voxels, replace = FALSE)
  surf_area_mat[suitable_voxels] <- unlimited_sa # surface area is not limiting
  microhab_mat[,,,1] <- surf_area_mat
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_equal(nrow(init_ind_df), total_nb_inds)
  expect_true(all(init_ind_df$Status == 1)) # alive

  # All individuals are in suitable voxels
  init_ind_df <- init_ind_df |> dplyr::mutate(
    "seq_index" = index_3d(X, Y, Z, dimensions[1], dimensions[2]),
    "is_in_suitable" = seq_index %in% suitable_voxels,
    "is_mature" = Mass >= species_df$MassAtMaturity[SpeciesID],
    "exptd_mass" =  species_df$MaximumMass[SpeciesID] *
        (1 - exp(-species_df$GrowthRate[SpeciesID] * Age)
    )
  )
  expect_true(all(init_ind_df$is_in_suitable))

  # Age and mass satisfy the growth equation
  nb_mature_inds <- nb_species * round(rnd_params$IndividualsPerSpecies *
                                   rnd_params$PercentageMaturePerSpecies / 100)
  expect_equal(sum(init_ind_df$is_mature), nb_mature_inds)
  expect_true(all(init_ind_df$exptd_mass))

  # what's enough sa?
  # set all individuals to same mass
  rnd_params$PercentageMaturePerSpecies <- 100
  fixed_mass <- runif(1, 0, 100)
  species_df$MaximumMass <- species_df$MassAtMaturity <- fixed_mass
  reqd_surf_area_per_ind <- fixed_mass ^ (2/3) / rnd_params$SurfaceBiomassScaling
  reqd_surf_area_total <- reqd_surf_area_per_ind * nb_species * rnd_params$IndividualsPerSpecies
  # Distribute it evenly among voxels
  surf_area_mat[suitable_voxels] <- reqd_surf_area_total / length(suitable_voxels)
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_true(all(init_ind_df$Status == 1)) # alive

  insufficient_surf_area <- reqd_surf_area_total - 1
  surf_area_mat[suitable_voxels] <- reqd_surf_area_total / length(suitable_voxels)
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_true(!all(init_ind_df$Status == 1))
  # expect a specific amount of dead individuals?

  # threshold: not enough SA for all
  # if not enough SA, only a certain fraction of individuals is allocated
  # if largest_first, only largest allocated
  # else random

  # light
  # if no suitable light conditions all individuals status 2 even if SA if fine
  # if variable conditions individuals only allocated where light in suitable range



  # testing microhab SA
  # (light always ok)
  # individuals are only present where SA > 0
  # individuals are only present where SA > min SA reqt

  # randomly sample which voxels contain SA

  # 1 - Enough SA to fit all individuals
  exptd_max_total_mass <- sum(species_df$MaximumMass * nb_mature_inds_per_sp +
    species_df$MassAtMaturity * (rnd_params$IndividualsPerSpecies - nb_mature_inds_per_sp))

  exptd_max_sa_occupied <- exptd_max_total_mass ^ (2/3) / rnd_params$SurfaceBiomassScaling
  total_sa_occupied <- sum(init_ind_df$SurfaceAreaOccupied)
  expect_lte(total_sa_occupied, exptd_max_sa_occupied)

  # X Y Z within boundaries of the landscape
  # X Y Z only in voxels with SA > 0 and light conditions ok
  # If enough surface area then all Status 1
  # If no SA all Status 2
  # sum(Mass <= MassAtMaturity) = nb_inds * (1- percentageMature / 100)
  # sum(Mass >= MassAtMaturity) = nb_inds * (percentageMature / 100)
  # age statisfies growth equation requirement

  # IndividualID in 1:nbinds
  # Species ID has all species id

  largest_inds_first = TRUE
  most_surf_area_first = TRUE


})
