test_that("Light suitability only", {

  # Initialise Microhabitat matrix
  dimensions <- sample(1:10, 3, replace = TRUE)
  nb_voxels <- prod(dimensions[1:3])
  Microhabitat <- create_empty_microhabitat(dimensions)

  # Initialise SpeciesPool
  nb_species <- 3
  SpeciesPool <- create_rnd_species_df(nb_species)

  # Light below niche mininum
  layer_light <- rep(min(SpeciesPool$MinLight) / 2, nb_voxels)
  Microhabitat[,,,3] <- layer_light
  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )
  # Nb of species, all between 0 and 1, etc...
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_equal(sum(suit_mat), 0)

  # Light above niche maximum
  layer_light <- rep(max(SpeciesPool$MaxLight) * 1.5, nb_voxels)
  Microhabitat[,,,3] <- layer_light
  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_equal(sum(suit_mat), 0)

  # Optimal
  layer_light <- rep(max(SpeciesPool$OptimumLight), nb_voxels)
  Microhabitat[,,,3] <- layer_light
  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_equal(sum(suit_mat), nb_voxels) # that is, 1 everywhere

  # Within light niche range
  layer_light <- runif(
    nb_voxels,
    min = min(SpeciesPool$MinLight),
    max  = max(SpeciesPool$MaxLight)
  )
  Microhabitat[,,,3] <- layer_light
  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_true(all(suit_mat >= 0) & all(suit_mat <= 1))

})

test_that("Microclimate suitability", {

  # Initialise Microhabitat matrix
  dimensions <- sample(1:10, 3, replace = TRUE)
  nb_voxels <- prod(dimensions[1:3])
  Microhabitat <- create_empty_microhabitat(dimensions)

  # Initialise SpeciesPool
  nb_species <- 3
  SpeciesPool <- create_rnd_species_df(nb_species)

  layer_light <- rep(min(SpeciesPool$MinLight) / 2, nb_voxels)
  Microhabitat[,,,3] <- layer_light

  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )

  # Nb of species, all between 0 and 1, etc...
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_equal(sum(suit_mat), 0)
})
