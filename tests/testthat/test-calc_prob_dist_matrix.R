test_that("consistent with previous version", {

  dims <- c(30, 30, 80)
  expanded_dims <- dims * 2 + 1
  centralPoint <- dims + 1

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species)
  )
  mat_exptd <- old_compute_prob_matrix_norm(
    centralPoint, expanded_dims[1], expanded_dims[2], expanded_dims[3], nrow(SpeciesPool), SpeciesPool
    )
  mat_obs <- calc_prob_disp_matrix(dims, SpeciesPool)
  testthat::expect_equal(mat_exptd, mat_obs)
  testthat::expect_equal(sum(mat_obs[,,,1]), 1)
})

test_that("matrix of with 1-dimension are accepted", {

  dims <- rep(1, 3)
  centralPoint <- floor(dims/2) + 1

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species)
  )
  testthat::expect_silent(disp_mat <- calc_prob_disp_matrix(
    centralPoint, dims, SpeciesPool
    ))
  testthat::expect_true(is.array(disp_mat))
  testthat::expect_equal(sum(disp_mat[,,,1]), 1)
  testthat::expect_equal(sum(disp_mat), nb_species)
})


test_that("Asymmtry", {

  dims <- sample(1:10, 3, replace = TRUE)
  expanded_dims <- dims * 2 + 1
  centre <- dims + 1

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = c(0, 0.5, 1)
  )

  disp_mat <- calc_prob_disp_matrix(dims, SpeciesPool)

  # Probabilities are always horizontally symmetric
  for (sp in 1:nb_species) {

    central_horizontal_slice <- disp_mat[ , , centre[3], sp]
    expect_equal(
      central_horizontal_slice[1, 1],
      central_horizontal_slice[expanded_dims[1], expanded_dims[2]]
    )
    expect_equal(
      central_horizontal_slice[1, 1],
      central_horizontal_slice[1, expanded_dims[2]]
    )
    expect_equal(
      central_horizontal_slice[1, 1],
      central_horizontal_slice[expanded_dims[1], 1]
    )
  }

  # Probabilities are vertically asymmetric (gravity)
  central_columns <- disp_mat[centre[1], centre[2],,]

  # Asymmetry = 0 -> strictly disperse upward
  sp <- 1
  expect_equal(central_columns[1,sp], 0)

  # Asymmetry = 0.5 -> symmetric dispersal
  sp <- 2
  expect_equal(central_columns[1, sp], central_columns[expanded_dims[3], sp])

  # Asymmetry = 1 -> strictly disperse downward
  sp <- 3
  expect_equal(central_columns[expanded_dims[3], sp], 0)

  })

test_that("Wind facilitates dispersal", {

  dims <- rep(3, 3)
  centre <- dims + 1
  expanded_dims <- dims * 2 + 1

  nb_species <- 1
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species),
    "DispersalKernelWindEffect" = runif(nb_species) # positive number
  )

  # Without wind
  disp_mat <- calc_prob_disp_matrix(dims, SpeciesPool)

  # With wind
  wind_val <- runif(1, min = 0, max = 100)
  wind_layer <- array(wind_val, dim = dims)
  disp_mat_wind <- calc_prob_disp_matrix(dims, SpeciesPool, wind_layer)

  testthat::expect_equal(sum(disp_mat_wind), sum(disp_mat))

  disp_mat[centre[1], centre[2], centre[3],1] >
    disp_mat_wind[centre[1], centre[2], centre[3], 1]

  disp_mat[expanded_dims[1], expanded_dims[2], expanded_dims[3], 1] <
    disp_mat_wind[expanded_dims[1], expanded_dims[2], expanded_dims[3], 1]

  disp_mat[1:7,4,1:7,1]
  disp_mat_wind[1:7,4,1:7,1]

  })

