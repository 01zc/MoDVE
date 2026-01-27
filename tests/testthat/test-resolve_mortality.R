source("../tests/test-utils.R")

test_that("multiplication works", {

  # Set Microhabitat matrix
  dimensions <- c(5, 4, 2)
  surf_area_loss <- rep(0, prod(dimensions))
  Microhabitat <- array(0, dim = c(dimensions, 3))
  # Mortality doesn't check actual surface area so we don't need to set it

  # Set surface area loss > 0 for two random cells
  rnd_sa_loss_vox <- sample(1:prod(dimensions), 2)
  rand_sa_loss <- sort(runif(2, 0, 1))
  surf_area_loss[rnd_sa_loss_vox] <- rand_sa_loss
  Microhabitat[,,,2] <- surf_area_loss

  # Set light niche
  SpeciesPool <- tibble::tibble(
    "SpeciesID" = 1,
    "MinLight" = sample(0:50, 1),
    "MaxLight" = MinLight * runif(1, 1.05, 2)
  )
  OptimumLight <- (SpeciesPool$MinLight + SpeciesPool$MaxLight) / 2
  light_layer <- rep(OptimumLight, prod(dimensions))

  # Insufficient light conditions in some voxels
  rand_poor_light_vox <- sample(1:prod(dimensions)[-rnd_sa_loss_vox], 2)
  light_layer[rand_poor_light_vox] <- SpeciesPool$MinLight * 0.8
  Microhabitat[,,,3] <- light_layer

  voxels_sa_loss <- arrayInd(which(Microhabitat[,,, 2] > 0), c(dimensions))

  nb_inds <- 5000
  E_initial <- tibble::tibble(
    # Put 1000 individuals in each voxel with some SA loss
    "X" = c(rep(voxels_sa_loss[1, 1], nb_inds), rep(voxels_sa_loss[2, 1], nb_inds)),
    "Y" = c(rep(voxels_sa_loss[1, 2], nb_inds), rep(voxels_sa_loss[2, 2], nb_inds)),
    "Z" = c(rep(voxels_sa_loss[1, 3], nb_inds), rep(voxels_sa_loss[2, 3], nb_inds)),
    "Status" = rep(1, nb_inds * 2),
    "SpeciesID" = rep(1, nb_inds * 2),
    "Mass" = runif(nb_inds * 2, 1, 100)
  )

  E <- resolve_mortality(
    E_initial,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = 0,
    MortRateMass = 0,
    MortRateMassScaling = 0
  )

  prop_dead_branch_fall <- E |> dplyr::group_by(X, Y, Z) |>
    dplyr::filter(Status == 3) |>
    dplyr::summarise("prop_dead_branch_fall" = dplyr::n() / nb_inds) |>
    dplyr::pull(prop_dead_branch_fall) |> sort()

  # Deaths to branch fall match proportion of surface area loss
  expect_equal(prop_dead_branch_fall, rand_sa_loss, tolerance = 0.03)

  E |> dplyr::filter()

  E <- resolve_mortality(
    E_initial,
    SpeciesPool,
    Microhabitat,
    use_mass_dep_mortality = FALSE,
    MortRateRandom = 0,
    MortRateMass = 0,
    MortRateMassScaling = 0
  )

  # Add dead individuals to make sure they don't die twice
  nb_already_dead <- 100
  E <- E |> dplyr::bind_rows(tibble::tibble(
    "X" = rep(voxels_sa_loss[1, 1], nb_already_dead),
    "Y" = rep(voxels_sa_loss[1, 2], nb_already_dead),
    "Z" = rep(voxels_sa_loss[1, 3], nb_already_dead),
    # Status 2 to 5 correspond to dead individuals, and 3 to death to branch loss
    "Status" = sample(c(2, 4, 5), nb_already_dead, replace = TRUE),
    "SpeciesID" = rep(1, nb_already_dead),
    "Mass" = runif(nb_already_dead, 1, 100)
  ))


  # for voxels where light < MinLight or > MaxLight all individuals status 4

  # sa_loss mortality takes precedence over other forms of mortality
  # if light out of conditions and sa_loss > 0
  # then fro 1000 individuals nb status 4 = (1-sa_loss)

  # assuming sa_loss 0 and !use_mass_dep_mortality
  # then nb status 5 = MortRateRandom

  # if use_mass_dep_mortality
  # Mass^MortRateMassScaling * MortRateMass


  # individuals with status > 1 are not affected by death

})
