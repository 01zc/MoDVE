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

  rnd_params <- draw_rnd_initial_inds_params()
  rnd_params$ScalingPerHa <- FALSE
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
  microhab_mat[,,,3] <- 50 # optimal light conditions

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

  # Age and mass satisfy the growth equation
  nb_mature_inds <- nb_species * round(rnd_params$IndividualsPerSpecies *
                                   rnd_params$PercentageMaturePerSpecies / 100)
  expect_equal(sum(init_ind_df$is_mature), nb_mature_inds)
  expect_equal(init_ind_df$Mass, init_ind_df$exptd_mass, tolerance = 0.1)

  # Individuals are allocated only if surface area is sufficient
  # Must fix the mass to a constant to determine necessary SA per individual
  rnd_params$PercentageMaturePerSpecies <- 100
  fixed_mass <- runif(1, 0, 100)
  species_df$MaximumMass <- species_df$MassAtMaturity <- fixed_mass
  reqd_surf_area_per_ind <- fixed_mass ^ (2/3) / rnd_params$SurfaceBiomassScaling +
    0.0001
  surf_area_mat <- rep(0, prod(dimensions))
  # for the following tests there must be less individuals than voxels
  if (total_nb_inds >= prod(dimensions)) {
    rnd_params$IndividualsPerSpecies <- floor(prod(dimensions) / nb_species) - 1
    total_nb_inds <- nb_species * rnd_params$IndividualsPerSpecies
  }
  expect_lt(total_nb_inds, prod(dimensions))

  # Enough SA for all
  suitable_voxels <- sample(1:prod(dimensions), size = total_nb_inds, replace = FALSE)
  surf_area_mat[suitable_voxels] <- reqd_surf_area_per_ind
  microhab_mat[,,,1] <- surf_area_mat
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_true(all(init_ind_df$Status == 1)) # alive

  # Enough SA for all but one
  surf_area_mat[suitable_voxels[1]] <- 0
  microhab_mat[,,,1] <- surf_area_mat
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat)
  expect_equal(sum(init_ind_df$Status == 2), 1)

  # If mass-based competition is enabled, largest individuals are allocated first
  # Set some variation, with our fixed mass as a lower threshold
  species_df$MaximumMass <- fixed_mass * 1.2
  not_enough_sa <- total_nb_inds * fixed_mass ^ (2/3) /
    rnd_params$SurfaceBiomassScaling
  surf_area_mat <- rep(0, prod(dimensions))
  rand_voxel <- sample(1:length(surf_area_mat), 1)
  surf_area_mat[rand_voxel] <- not_enough_sa
  microhab_mat[,,,1] <- surf_area_mat
  init_ind_df <- draw_initial_individuals(largest_inds_first = TRUE,
    rnd_params, species_df, microhab_mat)
  # Largest individual made it, smallest one did not
  expect_equal(init_ind_df$Status[which.max(init_ind_df$Mass)], 1)
  expect_equal(init_ind_df$Status[which.min(init_ind_df$Mass)], 2)

  # Individuals are not allocated to unsuitable light conditions
  # even if surface area is sufficient
  dimensions <- rep(6, 3)
  nb_bad_voxels <- prod(dimensions) / 3
  unlimited_sa <- 100000
  microhab_mat <- array(0, dim = c(dimensions, 3))
  microhab_mat[,,,1] <- unlimited_sa
  min_light <- runif(1, 1, 49)
  max_light <- runif(1, 51, 100)
  opt_light <- 50
  species_df$MinLight <- min_light
  species_df$MaxLight <- max_light
  species_df$OptimumLight <- opt_light
  voxels_not_enough_light <- sample(1:prod(dimensions), nb_bad_voxels)
  voxels_too_much_light <- sample((1:prod(dimensions))[-voxels_not_enough_light], nb_bad_voxels)
  expect_true(!any(voxels_too_much_light %in% voxels_not_enough_light))
  light_mat <- rep(opt_light, prod(dimensions))
  light_mat[voxels_not_enough_light] <- min_light - 0.5
  light_mat[voxels_too_much_light] <- max_light + 0.5
  microhab_mat[,,,3] <- light_mat
  init_ind_df <- draw_initial_individuals(rnd_params, species_df, microhab_mat) |>
    dplyr::mutate(
      "vox_idx" = index_3d(X, Y, Z, dimensions[1], dimensions[2])
    )
  expect_true(all(init_ind_df$Status == 1))
  expect_true(!any(init_ind_df$vox_idx %in% voxels_not_enough_light))
  expect_true(!any(init_ind_df$vox_idx %in% voxels_too_much_light))
})
