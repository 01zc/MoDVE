test_that("consistent with previous version", {

  dims <- c(30, 30, 80)
  centralPoint <- floor(dims/2) + 1

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species)
  )
  mat_exptd <- old_compute_prob_matrix_norm(
    centralPoint, dims[1], dims[2], dims[3], nrow(SpeciesPool), SpeciesPool
    )
  mat_obs <- calc_prob_disp_matrix(centralPoint, dims, SpeciesPool)
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
