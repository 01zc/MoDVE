test_that("Competition works as expected", {

  # Initialise habitat matrix
  dimensions <- sample(2:10, 3, replace = TRUE)
  Microhabitat <- array(0, dim = c(dimensions, 3))

  # Habitat contains two suitable voxels, one oversaturated with epiphytes...
  oversatd_vox <- c(
    sample(1:dimensions[1], 1),
    sample(1:dimensions[2], 1),
    sample(1:dimensions[3], 1)
    )
  # ...The other could host more
  undersatd_vox <- oversatd_vox
  while(identical(undersatd_vox, oversatd_vox)) { # different voxels
    undersatd_vox <- c(
      sample(1:dimensions[1], 1),
      sample(1:dimensions[2], 1),
      sample(1:dimensions[3], 1)
    )
  }

  # Initialise individuals
  nb_inds_per_vox <- sample(10:100, 1)
  E <- E_init <- tibble::tibble(
    "is_oversaturated" = c(rep(TRUE, nb_inds_per_vox), rep(FALSE, nb_inds_per_vox)), # shortcut for counting below
    "X" = c(rep(oversatd_vox[1], nb_inds_per_vox), rep(undersatd_vox[1], nb_inds_per_vox)),
    "Y" = c(rep(oversatd_vox[2], nb_inds_per_vox), rep(undersatd_vox[2], nb_inds_per_vox)),
    "Z" = c(rep(oversatd_vox[3], nb_inds_per_vox), rep(undersatd_vox[3], nb_inds_per_vox)),
    "Status" = rep(1, nb_inds_per_vox * 2),
    "IndividualID" = 1:(nb_inds_per_vox * 2),
    "SurfaceAreaOccupied" = runif(nb_inds_per_vox * 2, 0, 100)
  ) |>
    dplyr::slice_sample(prop = 1) # shuffle rows

  # Set surface area so there is too little, or enough surface area respectively
  E |> dplyr::filter(is_oversaturated) |> dplyr::pull(SurfaceAreaOccupied) |> sum()
  total_sa_oversatd <- 0.9 * sum(E$SurfaceAreaOccupied[which(E$is_oversaturated)])
  total_sa_undersatd <- 1.1 * sum(E$SurfaceAreaOccupied[which(!E$is_oversaturated)])
  surf_area_layer <- rep(0, prod(dimensions))
  index_oversatd <- index_3d(oversatd_vox[1], oversatd_vox[2], oversatd_vox[3], dimensions[1], dimensions[2])
  surf_area_layer[index_oversatd] <- total_sa_oversatd
  index_undersatd <- index_3d(undersatd_vox[1], undersatd_vox[2], undersatd_vox[3], dimensions[1], dimensions[2])
  surf_area_layer[index_undersatd] <- total_sa_undersatd
  Microhabitat[,,,1] <- surf_area_layer

  E <- resolve_competition(E_init, Microhabitat, larger_first = FALSE)

  # After competition is resolved, epiphytes in oversaturated voxel have died
  # down such that voxel is no longer oversaturated.
  nb_dead_oversatd <- E |> dplyr::filter(is_oversaturated) |>
    dplyr::pull(Status) |> magrittr::equals(2) |> sum()
  occupied_sa <- E |> dplyr::filter(is_oversaturated, Status == 1) |>
    dplyr::pull(SurfaceAreaOccupied) |> sum()
  expect_gte(nb_dead_oversatd, 1)
  expect_lte(occupied_sa, total_sa_oversatd)

  # No death have occurred in the undersaturated voxel
  nb_dead_undersatd <- E |> dplyr::filter(!is_oversaturated) |>
    dplyr::pull(Status) |> magrittr::equals(2) |> sum()
  expect_equal(nb_dead_undersatd, 0)

  # larger_first --> lightest X are dead
  # otherwise a random succession (!= sort individuals)
  # take mass off some individuals and add it to others
  E <- resolve_competition(E_init, Microhabitat, larger_first = TRUE)
  which_dead <- E$IndividualID[which(E$Status == 2)]
  smallest_inds <- E |> dplyr::filter(is_oversaturated) |>
    dplyr::arrange(SurfaceAreaOccupied) |>
    dplyr::slice_head(n = length(which_dead)) |>
    dplyr::pull(IndividualID)
  # Dead individuals are also the lightest ones
  expect_setequal(which_dead, smallest_inds)

  # Individuals that are already dead are not affected by competition, and
  # don't contribute to it
  # Dump a bunch of dead epiphytes in the undersaturated voxel
  nb_dead <- sample(1:nb_inds_per_vox, 1)
  E_dead <- E_init |> dplyr::filter(!is_oversaturated) |>
    dplyr::slice_sample(n = nb_dead)
  E_dead$Status <- sample(2:5, nb_dead, replace = TRUE)
  E <- rbind(E_init, E_dead)
  E <- resolve_competition(E, Microhabitat, larger_first = FALSE)
  nb_dead_undersatd <- E |> dplyr::filter(!is_oversaturated) |>
    dplyr::pull(Status) |> dplyr::between(2, 5) |> sum()
  expect_equal(nb_dead_undersatd, nb_dead)

  # If the available surface area matches the occupied SA, no epiphyte die
  E <- E_init |> dplyr::filter(is_oversaturated)
  surf_area_layer <- rep(0, prod(dimensions))
  surf_area_layer[index_oversatd] <- sum(E$SurfaceAreaOccupied)
  Microhabitat[,,,1] <- surf_area_layer
  E <- resolve_competition(E, Microhabitat, larger_first = FALSE)
  nb_dead <- sum(E$Status == 2)
  expect_equal(nb_dead, 0)

  # If surface area is zero, all epiphytes die
  E <- E_init
  Microhabitat[,,,1] <- 0
  E <- resolve_competition(E, Microhabitat, larger_first = FALSE)
  all_dead <- all(E$Status == 2)
  expect_true(all_dead)

  # TODO: Does it work with a 1- or 2-dimension matrix?

})
