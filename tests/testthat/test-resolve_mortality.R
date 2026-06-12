test_that("Mortality works as expected", {

  nb_inds <- 5000
  rnd_mass <- runif(1, 1, 100)

  # Initialise Microhabitat matrix
  dimensions <- sample(1:5, 3, replace = TRUE)
  sa_loss_layer <- rep(0, prod(dimensions))
  Microhabitat <- array(0, dim = c(dimensions, 3))
  # Mortality doesn't check surface area (only loss) so we don't need to set it

  # Initialise light niche
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = 1,
    "MinLight" = sample(0:50, 1),
    "MaxLight" = MinLight * runif(1, 1.05, 2)
  )
  OptimumLight <- (SpeciesPool$MinLight + SpeciesPool$MaxLight) / 2
  Microhabitat[,,,3] <- OptimumLight

  # If no SA loss / light issue, death tally corresponds to mortality rate
  # (with 5% tolerance for stochasticity vs sample size)
  MortRateRandom <- runif(1, 0, 1)
  nb_inds <- 15000
  E_initial <- tibble::tibble(
    "X" = sample(1:dimensions[1], nb_inds, replace = TRUE),
    "Y" = sample(1:dimensions[2], nb_inds, replace = TRUE),
    "Z" = sample(1:dimensions[3], nb_inds, replace = TRUE),
    "Status" = rep(1, nb_inds),
    "SpeciesID" = rep(1, nb_inds),
    "Mass" = rnd_mass
  )
  E <- resolve_mortality(
    E_initial,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = MortRateRandom,
    MortRateMass = 0,
    MortRateMassScaling = 0
  )
  prop_dead <- sum(E$Status == 5) / nb_inds
  expect_equal(prop_dead, MortRateRandom, tolerance = 0.05)

  # When mass-dependent mortality is enabled,
  # death tally is conform to mass-dependence equation
  # (with 5% tolerance for stochasticity vs sample size)
  MortRateMass <- runif(1, 0, 1)
  MortRateMassScaling <- -rgamma(1, 5, 5)
  exptd_mortality_rate <- rnd_mass^MortRateMassScaling * MortRateMass
  E <- resolve_mortality(
    E_initial,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = TRUE,
    MortRateRandom = 0,
    MortRateMass = MortRateMass,
    MortRateMassScaling = MortRateMassScaling
  )
  prop_dead <- sum(E$Status == 5) / nb_inds
  expect_equal(prop_dead, exptd_mortality_rate, tolerance = 0.05)

  # Individuals who are already dead are not affected by mortality
  E <- E_initial
  E$Status <- sample(2:5, nrow(E), replace = TRUE)
  status_before <- E$Status
  E <- resolve_mortality(
    E,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = MortRateRandom,
    MortRateMass = 0,
    MortRateMassScaling = 0
  )
  status_after <- E$Status
  expect_equal(status_before, status_after)

  # Mortality from surface area loss and light conditions

  # Set surface area loss > 0 for two random cells
  rnd_sa_loss_vox <- sample(1:prod(dimensions), 2)
  rand_sa_loss <- sort(runif(2, 0, 1))
  sa_loss_layer[rnd_sa_loss_vox] <- rand_sa_loss
  Microhabitat[,,,2] <- sa_loss_layer

  # One of the two voxels with SA loss also has insufficient light
  light_layer <- rep(OptimumLight, prod(dimensions))
  rand_poor_light_vox <- rnd_sa_loss_vox[sample(1:2, 1)]
  light_layer[rand_poor_light_vox] <- SpeciesPool$MinLight * 0.8
  Microhabitat[,,,3] <- light_layer

  voxels_sa_loss <- arrayInd(which(Microhabitat[,,, 2] > 0), c(dimensions))
  voxel_insufficient_light <- arrayInd(which(Microhabitat[,,, 3] < SpeciesPool$MinLight), c(dimensions))

  E <- tibble::tibble(
    # Put 1000 individuals in each voxel with some SA loss
    "X" = c(rep(voxels_sa_loss[1, 1], nb_inds), rep(voxels_sa_loss[2, 1], nb_inds)),
    "Y" = c(rep(voxels_sa_loss[1, 2], nb_inds), rep(voxels_sa_loss[2, 2], nb_inds)),
    "Z" = c(rep(voxels_sa_loss[1, 3], nb_inds), rep(voxels_sa_loss[2, 3], nb_inds)),
    "Status" = rep(1, nb_inds * 2),
    "SpeciesID" = rep(1, nb_inds * 2),
    "Mass" = rnd_mass
  )

  E <- resolve_mortality(
    E,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = 0,
    MortRateMass = 0,
    MortRateMassScaling = 0
  )

  prop_dead_branch_fall <- E |> dplyr::group_by(X, Y, Z) |>
    dplyr::filter(Status == 3) |>
    dplyr::summarise(
      "prop_dead_branch_fall" = dplyr::n() / nb_inds,
      .groups = "keep"
      ) |>
    dplyr::pull(prop_dead_branch_fall) |> sort()

  # Deaths to branch fall match proportion of surface area loss
  expect_equal(prop_dead_branch_fall, rand_sa_loss, tolerance = 0.03)

  remaining_individuals <- E |> dplyr::filter(Status != 3) |>
    dplyr::mutate(
      "not_enough_light" =
        X == voxel_insufficient_light[,1] &
        Y == voxel_insufficient_light[,2] &
        Z == voxel_insufficient_light[,3]
    )

  # All other individuals in voxel with not enough light die from it
  # but surface area loss takes precedence as a source of mortality
  all_dead_from_light <- remaining_individuals |>
    dplyr::filter(not_enough_light) |>
    dplyr::pull(Status) |>
    magrittr::equals(4) |>
    all()
  expect_true(all_dead_from_light)

  # All remaining individuals in voxel with enough light are still alive
  all_alive <- remaining_individuals |>
    dplyr::filter(!not_enough_light) |>
    dplyr::pull(Status) |>
    magrittr::equals(1) |>
    all()
  expect_true(all_alive)

})

create_empty_microhab <- function(dim) {
  Microhabitat <- array(0, dim = dim)
  attr(Microhabitat, "layer_mapping") <- c("surface_area", "surface_area_loss", "light")
  return(Microhabitat)
}

test_that("Niche-related mortality", {

  climatic_vars <- c("humidity", "temperature", "wind")
  trait_names <- c("Hum", "Temp", "Wind")
  names(trait_names) <- climatic_vars
  status_dead <- 6:8
  names(status_dead) <- climatic_vars

  # Species die based on their own niche
  for (var in climatic_vars) {

    trait_name <- trait_names[var]

    # Non-overlapping climate niches
    niche_diff <- 10
    niche_min_sp1 <- runif(1, 0, 50)
    niche_max_sp1 <- niche_min_sp1 * runif(1, 1.05, 2)
    niche_min_sp2 <- niche_max_sp1 + niche_diff
    niche_max_sp2 <- niche_min_sp2 * runif(1, 1.05, 2)

    # Initialise niche
    SpeciesPool <- tibble::tibble(
      "SpeciesID" = c(1, 2),
      "MinLight" = 0,
      "MaxLight" = 0,
      "Min{trait_name}" := c(niche_min_sp1, niche_min_sp2),
      "Max{trait_name}" := c(niche_max_sp1, niche_max_sp2)
    )

    # Environment is suitable for sp 1, unsuitable for sp2
    microhabitat_value <- niche_min_sp1 + niche_diff / 2
    # Initialise Microhabitat matrix
    dimensions <- c(sample(1:5, 3, replace = TRUE), 4)
    Microhabitat <- create_empty_microhab(dimensions)
    attr(Microhabitat, "layer_mapping") <- c(attr(Microhabitat, "layer_mapping") , var)
    Microhabitat[,,,4] <- microhabitat_value

    # Distribute individuals randomly across space and species
    nb_inds <- 15000
    E <- tibble::tibble(
      "X" = sample(1:dimensions[1], nb_inds, replace = TRUE),
      "Y" = sample(1:dimensions[2], nb_inds, replace = TRUE),
      "Z" = sample(1:dimensions[3], nb_inds, replace = TRUE),
      "Status" = 1,
      # Attribute randomly to species 1 or 2
      "SpeciesID" = sample(1:2, nb_inds, replace = TRUE),
      "Mass" = 0 # needs to exist but we don't use it here
    )
    E <- dplyr::slice_sample(E, prop = 1) # shuffle

    E <- resolve_mortality(
      E,
      SpeciesPool,
      Microhabitat,
      use_mass_dep_mortality = FALSE,
      MortRateRandom = 0, MortRateMass = 0, MortRateMassScaling = 0
    )

    testthat::expect_true(all(E[E$SpeciesID == 1, "Status"] == 1))
    testthat::expect_true(all(E[E$SpeciesID == 2, "Status"] == status_dead[var]))
  }

  # Individuals die only where conditions are unsuitable
  for (var in climatic_vars) {

    trait_name <- trait_names[var]

    # Non-overlapping climate niches
    niche_min <- runif(1, 0, 50)
    niche_max <- niche_min * runif(1, 1.05, 2)

    # Initialise niche
    SpeciesPool <- tibble::tibble(
      "SpeciesID" = 1,
      "MinLight" = 0,
      "MaxLight" = 0,
      "Min{trait_name}" := niche_min,
      "Max{trait_name}" := niche_max
    )

    # Initialise Microhabitat matrix
    dimensions <- c(sample(1:5, 3, replace = TRUE), 4)
    Microhabitat <- create_empty_microhab(dimensions)
    attr(Microhabitat, "layer_mapping") <- c(attr(Microhabitat, "layer_mapping") , var)
    nb_cells <- prod(dimensions[1:3])

    # Only a random subset of cells are unsuitable
    val_suitable <- (niche_min + niche_max) / 2
    val_unsuitable <- niche_min * 0.9
    nb_unsuitable <- round(nb_cells * runif(1))
    layer_values <- rep(val_suitable, nb_cells)
    index_unsuitable <- sample(1:nb_cells, size = nb_unsuitable)
    layer_values[index_unsuitable] <- val_unsuitable
    Microhabitat[,,,4] <- layer_values

    # Distribute individuals randomly across space and species
    nb_inds <- 15000
    E <- tibble::tibble(
      "X" = sample(1:dimensions[1], nb_inds, replace = TRUE),
      "Y" = sample(1:dimensions[2], nb_inds, replace = TRUE),
      "Z" = sample(1:dimensions[3], nb_inds, replace = TRUE),
      "Status" = 1,
      # Attribute randomly to species 1 or 2
      "SpeciesID" = 1,
      "Mass" = 0 # needs to exist but we don't use it here
    )
    E <- dplyr::slice_sample(E, prop = 1) # shuffle

    E <- resolve_mortality(
      E,
      SpeciesPool,
      Microhabitat,
      use_mass_dep_mortality = FALSE,
      MortRateRandom = 0, MortRateMass = 0, MortRateMassScaling = 0
    )

    E <- E |>
      dplyr::mutate(
        "cell_index" = index_3d(X, Y, Z, dimensions[1], dimensions[2]),
        "should_die" = cell_index %in% index_unsuitable
      )

    testthat::expect_true(all(E[!E$should_die, "Status"] == 1))
    testthat::expect_true(all(E[E$should_die, "Status"] == status_dead[var]))
  }

  ##  Multiple climatic layers

  # Partially-overlapping climate niches
  niche_hum_min <- runif(1, 0, 50)
  niche_hum_max <- niche_hum_min * runif(1, 1.05, 2)
  niche_wind_min <- niche_hum_min
  niche_wind_max <- (niche_hum_min + niche_hum_max) / 2
  expect_lt(niche_wind_max, niche_hum_max)
  val_suitable_hum <- (niche_hum_min + niche_hum_max) / 2
  # Suitable for both
  val_suitable <- (niche_wind_min + niche_wind_max) / 2
  # Unsuitable for both
  val_unsuitable <- niche_hum_max * 1.5

  # Initialise niche
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = 1,
    "MinLight" = 0,
    "MaxLight" = 0,
    "MinWind" = niche_wind_min, "MaxWind" = niche_wind_max,
    "MinHum" = niche_hum_min, "MaxHum" = niche_hum_max
  )

  dimensions <- c(sample(1:5, 3, replace = TRUE), 5)

  # Distribute individuals randomly
  nb_inds <- 100
  E_init <- tibble::tibble(
    "X" = sample(1:dimensions[1], nb_inds, replace = TRUE),
    "Y" = sample(1:dimensions[2], nb_inds, replace = TRUE),
    "Z" = sample(1:dimensions[3], nb_inds, replace = TRUE),
    "Status" = 1,
    # Attribute randomly to species 1 or 2
    "SpeciesID" = 1,
    "Mass" = 0 # needs to exist but we don't use it here
  )
  E_init <- dplyr::slice_sample(E_init, prop = 1) # shuffle

  Microhabitat <- create_empty_microhab(dimensions)
  attr(Microhabitat, "layer_mapping") <- c(
    attr(Microhabitat, "layer_mapping") , "wind", "humidity"
  )

  # 1/2 - Suitable conditions for humidity but not wind
  # All individuals should die from wind mortality
  Microhabitat[,,,5] <- val_suitable
  Microhabitat[,,,4] <- val_unsuitable
  E <- resolve_mortality(
    E_init,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = 0, MortRateMass = 0, MortRateMassScaling = 0
  )
  testthat::expect_true(all(E[, "Status"] == status_dead["wind"]))


  # 2/2 - Unsuitable conditions for both
  # Humidity resolved before wind so all die from wind
  Microhabitat[,,,5] <- val_unsuitable
  Microhabitat[,,,4] <- val_unsuitable
  E <- resolve_mortality(
    E_init,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = 0, MortRateMass = 0, MortRateMassScaling = 0
  )
  testthat::expect_true(all(E[, "Status"] == status_dead["humidity"]))
})


