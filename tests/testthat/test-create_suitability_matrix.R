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
  which_sp <- 1
  layer_light <- rep(SpeciesPool$OptimumLight[which_sp], nb_voxels)
  Microhabitat[,,,3] <- layer_light
  suit_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE
  )
  testthat::expect_silent(check_suitability(suit_mat, dimensions, nb_species))
  expect_equal(sum(suit_mat[,,,which_sp]), nb_voxels) # that is, 1 everywhere

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
  expect_true(sum(suit_mat[,,,1]) != sum(suit_mat[,,,2]))
})

test_that("Microclimate suitability", {

  microclimate_vars <- c("temperature", "humidity", "wind")
  microclimate_vars <- microclimate_vars[sample(1:3)] # shuffle
  idx_hum <- 3 + which(microclimate_vars == "humidity")
  idx_wind <- 3 + which(microclimate_vars == "wind")
  idx_temp <- 3 + which(microclimate_vars == "temperature")

  # Initialise Microhabitat matrix
  dimensions <- sample(1:10, 3, replace = TRUE)
  nb_voxels <- prod(dimensions[1:3])
  Microhabitat <- create_empty_microhabitat(dimensions, microclimate_vars)

  # Initialise SpeciesPool
  nb_species <- 1
  SpeciesPool <- create_rnd_species_df(nb_species)

  # General conditions - intermediate suitability
  layer_light <- runif(nb_voxels, SpeciesPool$MinLight, SpeciesPool$MaxLight)
  layer_hum <- runif(nb_voxels, SpeciesPool$MinHum, SpeciesPool$MaxHum)
  layer_temp <- runif(nb_voxels, SpeciesPool$MinTemp, SpeciesPool$MaxTemp)
  layer_wind <- runif(nb_voxels, SpeciesPool$MinWind, SpeciesPool$MaxWind)

  # Special use cases in specific voxels
  # 1 - one layer is unsuitable
  rnd_unsuitable <- sample(1:nb_voxels, 1)
  layer_temp[rnd_unsuitable] <- SpeciesPool$MinTemp / 2
  # 2 - all layers are optimal
  rnd_perfect <- sample((1:nb_voxels)[-rnd_unsuitable], 1)
  layer_temp[rnd_perfect] <- SpeciesPool$OptimumTemp
  layer_hum[rnd_perfect] <- SpeciesPool$OptimumHum
  layer_wind[rnd_perfect] <- SpeciesPool$OptimumWind
  layer_light[rnd_perfect] <- SpeciesPool$OptimumLight

  Microhabitat[,,,3] <- layer_light
  Microhabitat[,,,idx_hum] <- layer_hum
  Microhabitat[,,,idx_temp] <- layer_temp
  Microhabitat[,,,idx_wind] <- layer_wind

  suitability_mat <- create_suitability_matrix(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = FALSE
  )

  testthat::expect_silent(check_suitability(suitability_mat, dimensions, nb_species))
  expect_equal(suitability_mat[rnd_unsuitable], 0)
  expect_equal(suitability_mat[rnd_perfect], 1)
  expect_true(all(suitability_mat[-c(rnd_perfect, rnd_unsuitable)] > 0))
  expect_true(all(suitability_mat[-c(rnd_perfect, rnd_unsuitable)] < 1))
})
