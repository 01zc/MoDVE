source("../old_funcs.R")

test_that("consistent with previous version", {

  dimX <- dimY <- 30
  dimZ <- 80
  centralPoint <- find_central_point(c(dimX, dimY, dimZ))

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species)
  )
  mat_exptd <- old_compute_prob_matrix_norm(centralPoint, dimX, dimY, dimZ, nrow(SpeciesPool), SpeciesPool)
  mat_obs <- calc_prob_disp_matrix(centralPoint, dimX, dimY, dimZ,  SpeciesPool)
  testthat::expect_equal(mat_exptd, mat_obs)
})

test_that("matrix of with 1-dimension are accepted", {

  dimX <- dimY <- 1
  dimZ <- 1
  centralPoint <- find_central_point(c(dimX, dimY, dimZ))

  nb_species <- 3
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = seq_len(nb_species),
    "DispersalKernel" = runif(nb_species),
    "DispersalKernelAsymmetry" = runif(nb_species)
  )
  testthat::expect_silent(disp_mat <- calc_prob_disp_matrix(centralPoint, dimX, dimY, dimZ,  SpeciesPool))
  testthat::expect_true(is.array(disp_mat))
  testthat::expect_equal(sum(disp_mat[,,,1]), 1)
  testthat::expect_equal(sum(disp_mat), nb_species)
})
