test_that("Returns the correct amount of surface area", {

  surface_biomass_scaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimensions <- sample(1:20, 3)
  surface_area_mat <- array(10, dim = dimensions)
  nb_inds <- sample(1:dimensions[1], 1)

  epiphyte_tbl <- tibble::tibble(
    # Initialise individuals randomly
    X = rep(sample(dimensions[1], nb_inds, replace = FALSE)),
    # ^all individuals are in a different voxel
    Y = rep(sample(dimensions[2], nb_inds, replace = TRUE)),
    Z = rep(sample(dimensions[3], nb_inds, replace = TRUE)),
    Mass = runif(nb_inds, 0, 1),
    # Compute expected surface area
    expected_sa = Mass^(2/3) / surface_biomass_scaling
  )

  avail_sa_mat <- get_surf_area_mat(
    epiphyte_tbl,
    surface_area_mat,
    surface_biomass_scaling
    )

  # Correct for a single individual
  rand_ind <- dplyr::slice_sample(epiphyte_tbl, n = 1)
  sa_exptd <- 10 - rand_ind$expected_sa
  sa_obs <- avail_sa_mat[rand_ind$X, rand_ind$Y, rand_ind$Z]
  expect_equal(sa_obs, sa_exptd) # could be smaller

  # Correct overall
  expect_true(all(epiphyte_tbl$expected_sa < 10)) # otherwise SA could be zero
  expected_occ_sa <- sum(surface_area_mat) - sum(epiphyte_tbl$expected_sa)
  expect_equal(sum(avail_sa_mat), expected_occ_sa)

  # Surface area cannot be negative
  epiphyte_tbl$Mass <- 100000
  avail_sa_mat <- get_surf_area_mat(
    epiphyte_tbl,
    surface_area_mat,
    surface_biomass_scaling
  )
  expect_gte(min(avail_sa_mat), 0.0)

})
