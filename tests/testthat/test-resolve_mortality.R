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
