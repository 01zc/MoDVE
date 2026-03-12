test_that("abuse cases", {
  # TODO: here try to enter incorrect input esp. for shoot and trunk tables
})

# Diameter is constant through tests, doesn't influence
# which voxels the branch goes through etc.
shoot_diameter <- trunk_diameter <- 0.01

# Shorthand functions to help building input tables for unit tests
create_empty_shoot_tbl <- function() {
  return(tibble::tibble(
    "xbegin" = numeric(),
    "ybegin" = numeric(),
    "zbegin" = numeric(),
    "xend" = numeric(),
    "yend" = numeric(),
    "zend" = numeric(),
    "length" = numeric(),
    "diameter" = numeric(),
    "shootID" = numeric()
  ))
}

add_shoot_row <- function(shoot_tbl, begin_coords, end_coords) {
  return(shoot_tbl |>
    tibble::add_row(
      "xbegin" = begin_coords[1], "xend" = end_coords[1],
      "ybegin" = begin_coords[2], "yend" = end_coords[2],
      "zbegin" = begin_coords[3], "zend" = end_coords[3],
      "length" = sqrt((xend - xbegin)^2 + (yend - ybegin)^2 + (zend - zbegin)^2),
      "diameter" = shoot_diameter,
      "shootID" = nrow(shoot_tbl) + 1
    ))
}

create_empty_trunk_tbl <- function() {
  return(tibble::tibble(
    "x" = numeric(),
    "y" = numeric(),
    "height" = numeric(),
    "diameter" = numeric(),
    "treeID" = numeric()
  ))
}

test_that("Branch surface area is calculated correctly", {

  # A 5*5 landscape with a 1-cell corridor around it
  corridor <- 1
  dim <- 5
  grid_dim <- dim + 2 * corridor
  microhab_extent <- c(corridor, dim + corridor)

  config <- list(
    calcSurfaceArea = TRUE,
    calcSurfaceAreaLoss = FALSE,
    calcLightConditions = FALSE,
    calcWeightedAngles = FALSE,
    MaxX = dim,
    MaxY = dim,
    MaxZ = 1, # 2D
    corridor = corridor
  )

  # Create a test shoots table
  # Each case is a branch segment corresponding to an edge case we want to test
  {
    shoots_dt <- create_empty_shoot_tbl()
    z <- 0.5
    intersctd_voxels <- list()

    # Along x
    # Case 1: along x: segment starts and ends in corridor, beyond microhabitat area
    begin_coords <- c("x" = 0.5, "y" = microhab_extent[1] + 0.5, "z" = z)
    end_coords <- c("x" = grid_dim, "y" = microhab_extent[1] + 0.5, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[1]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 2: along y: segment goes in descending direction
    # intersects the first segment in [1, 1]
    begin_coords <- c("x" = microhab_extent[1] + 0.5, "y" = microhab_extent[2] - 0.5, "z" = z)
    end_coords <- c("x" = microhab_extent[1] + 0.5, "y" = microhab_extent[1] + 0.5, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[2]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 3: segment is smaller than a cell
    # intersects the first segment
    begin_coords <- c("x" = microhab_extent[2] - 0.8, "y" = microhab_extent[2] - 0.8, "z" = z)
    end_coords <- c("x" = microhab_extent[2] - 0.2, "y" = microhab_extent[2] - 0.2, "z" = z)
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[3]] <- find_intersecting_voxels(begin_coords, end_coords)

    # Case 4 and 5: segments on cell limits belong to the cell below
    begin_coords <- c("x" = 1.5, "y" = 1, "z" = z)
    end_coords <- c("x" = 2.5, "y" = 1, "z" = z)
    # this one doesn't contribute to microhabitat
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[4]] <- find_intersecting_voxels(begin_coords, end_coords)
    begin_coords <- c("x" = 2.5, "y" = 6, "z" = z)
    end_coords <- c("x" = 3.5, "y" = 6, "z" = z)
    # but this one does
    shoots_dt <- shoots_dt |> add_shoot_row(begin_coords, end_coords)
    intersctd_voxels[[5]] <- find_intersecting_voxels(begin_coords, end_coords)
  }

  # Landscape viz to help working with these unit tests
  #shoots_dt |>
  #  ggplot2::ggplot() +
  #  ggplot2::geom_rect(
  #    xmin = microhab_extent[1], xmax = microhab_extent[2],
  #    ymin = microhab_extent[1], ymax = microhab_extent[2],
  #    alpha = 0.1
  #  ) +
  #  ggplot2::geom_segment(
  #    ggplot2::aes(x = xbegin, y = ybegin, xend = xend, yend = yend, colour = as.factor(shootID))
  #  ) +
  #  ggplot2::coord_cartesian(xlim = c(0, grid_dim), ylim = c(0, grid_dim)) +
  #  ggplot2::theme_linedraw()

  # Expectations
  nb_voxels <- sapply(intersctd_voxels, length)
  surf_area_exptd <- shoots_dt$length * shoots_dt$diameter * pi / 2 / nb_voxels
  expected_mat <- matrix(data = 0, nrow = dim, ncol = dim)
  expected_mat[1:5, 1] <- surf_area_exptd[1]
  expected_mat[1, 1:5] <- expected_mat[1, 1:5] + surf_area_exptd[2]
  expected_mat[5, 5] <- surf_area_exptd[3]
  expected_mat[2:3, 5] <- surf_area_exptd[5]

  # Carry out test
  microhab_mat <- create_microhabitat_mat(
    config = config,
    shoot_dt = shoots_dt,
    trunk_dt = create_empty_trunk_tbl()
  )
  expect_equal(microhab_mat[,,1,1], expected_mat)
})

test_that("Trunk surface area is calculated correctly", {

  # A 2*2*5 landscape
  dim_xy <- 2
  dim_z <- 5
  config <- list(
    calcSurfaceArea = TRUE,
    calcSurfaceAreaLoss = FALSE,
    calcLightConditions = FALSE,
    calcWeightedAngles = FALSE,
    MaxX = dim_xy,
    MaxY = dim_xy,
    MaxZ = dim_z, # 2D
    corridor = 0
  )

  trunk_dt <- create_empty_trunk_tbl() |>
    tibble::add_row(
      # One sapling in one corner
      "x" = 0.5, "y" = 0.5, "height" = 0.5,
      "diameter" = trunk_diameter, "treeID" = 1
    ) |>
    tibble::add_row(
      # A small tree in the opposite corner
      "x" = 1.5, "y" = 1.5, "height" = 2.5,
      "diameter" = trunk_diameter, "treeID" = 2
    ) |> dplyr::mutate(
      "sa" = pi * diameter / 2 * sqrt((diameter/2)^2 + height^2)
    )

  # Create a single branch to assert both branch and trunk contribute to area
  shoots_dt <- create_empty_shoot_tbl() |>
    add_shoot_row(begin_coords = c(0.1, 0.5, 0.5), end_coords = c(0.5, 0.5, 0.5)) |>
    dplyr::mutate("sa" = length * diameter * pi / 2)

  # Carry out test
  microhab_mat <- create_microhabitat_mat(
    config = config,
    shoot_dt = shoots_dt,
    trunk_dt = trunk_dt
  )[,,,1] # only retain surface area

  # Total cone volume is calculated correctly
  expect_equal(sum(microhab_mat[2,2,1:3]), trunk_dt$sa[2])
  # Cone volume is distributed among crossed voxels
  expect_true(all(microhab_mat[2,2,1:3] > 0.0))
  expect_true(microhab_mat[2,2,3] < microhab_mat[2,2,2])
  expect_true(microhab_mat[2,2,2] < microhab_mat[2,2,1])

  # Trunk and branches contribute additively to surface area
  exptd_sa_corner <- shoots_dt$sa[1] + trunk_dt$sa[1]
  expect_equal(microhab_mat[1,1,1],  exptd_sa_corner)

  # Voxels above tree height remain empty
  expect_equal(sum(microhab_mat[1,1,2:5]), 0.0)
  expect_equal(sum(microhab_mat[2,2,4:5]), 0.0)

  # Other areas remain empty
  expect_equal(sum(microhab_mat[1, 2, 1:5]), 0.0)
  expect_equal(sum(microhab_mat[2, 1, 1:5]), 0.0)
})

test_that("Available light is calculated correctly", {

  # follows a Beer-Lambert extinction law
  voxel_area <- 10000

  # A 3*3*3 grid with a corridor
  dim_xy <- 3
  dim_z <- 3
  corridor <- 1
  dim_corr <- dim_xy + 2 * corridor

  config <- list(
    calcSurfaceArea = FALSE,
    calcSurfaceAreaLoss = FALSE,
    calcLightConditions = TRUE,
    calcWeightedAngles = FALSE,
    # At first only above voxels affect light availability
    DistVoxToConsider = 0,
    kL = exp(runif(1, -4, 1)), # reasonable values
    MaxX = dim_xy,
    MaxY = dim_xy,
    MaxZ = dim_z, # 2D
    corridor = corridor
  )

  vox_dt <- tidyr::expand_grid(
    "x" = seq_len(dim_corr),
    "y" = seq_len(dim_corr),
    "z" = seq_len(dim_z)
  )
  vox_dt$leafarea <- runif(n = nrow(vox_dt), min = 0, max = 5000)

  config$DistVoxToConsider <- 1
  central_coord <- ceiling(dim_corr / 2)
  vox_dt <- vox_dt |> dplyr::mutate(
    "ring_nb" = pmax(abs(x - central_coord), abs(y - central_coord)),
    "rel_contrib" = 1 / (config$DistVoxToConsider + 1) / pmax(1, ring_nb * 8),
  )

  # Compute contributions to the central cell
  exptd_light_mat <- array(dim = c(dim_corr, dim_corr))
  for (x in 1:dim_corr) {
    for (y in 1:dim_corr) {
      voxel_column <- vox_dt$x == x & vox_dt$y == y
      total_leaf_area <- sum(vox_dt$leafarea[voxel_column])
      exptd_light_mat[x,y] <- vox_dt$rel_contrib[voxel_column & vox_dt$z == 1] *
        exp(-config$kL * total_leaf_area / voxel_area)
    }
  }
  exptd_light_avail <- sum(exptd_light_mat)

  sum_contribs <- vox_dt |>
    dplyr::filter(ring_nb <= config$DistVoxToConsider, z == 1) |>
    dplyr::pull(rel_contrib) |> sum()
  if (sum_contribs != 1.0) {
    stop("Relative contributions don't sum to 1.")
  }

  # LR = 0 only focal voxel + those above it
  microhab_mat <- create_microhabitat_mat(
    config = config,
    shoot_dt = create_empty_shoot_tbl(),
    trunk_dt = create_empty_trunk_tbl(),
    vox_dt = vox_dt
  )[,,,3] # only retain light

  # Compute cumulated leaf area
  vox_dt$cumul_leaf_area <- vox_dt$leafarea
  for (z in (dim_z - 1):1) {
    for (x in seq_len(dim_corr)) {
      for (y in seq_len(dim_corr)) {
        row_above <- vox_dt$x == x & vox_dt$y == y & vox_dt$z == (z + 1)
        row <- vox_dt$x == x & vox_dt$y == y & vox_dt$z == z
        vox_dt$cumul_leaf_area[row] <- vox_dt$cumul_leaf_area[row] + vox_dt$cumul_leaf_area[row_above]
      }
    }
  }
  vox_dt$light <- exp(-config$kL * vox_dt$cumul_leaf_area / voxel_area)
  vox_dt$contrib <- vox_dt$rel_contrib * vox_dt$light

  exptd_light <- vox_dt |>
    dplyr::filter(
      ring_nb <= config$DistVoxToConsider,
      z == 1
      ) |>
    dplyr::pull(contrib) |> sum()

  expect_equal(microhab_mat[2, 2, 1], exptd_light)
})
