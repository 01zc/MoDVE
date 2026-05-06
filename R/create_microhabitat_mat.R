#' Assemble the microhabitat matrix from forest simulation outputs
#'
#' Read in forest stand simulation data from `MoF3D`, and compute the surface area,
#' surface area loss, and/or light conditions available for epiphytes.
#'
#' @param config a list of parameters with at least the following elements:
#' * `corridor` size (in number of voxels) of the corridor, that is a band
#' of voxels on the edge of forest plots without trees.
#' * `MaxX` maximum coordinate of the forest plot along the x direction
#' * `MaxY` maximum coordinate of the forest plot along the y direction
#' * `MaxZ` maximum coordinate of the forest plot along the z direction
#' *  `Imax` maximum light intensity above the canopy
#' * `kL` light extinction coefficient
#' * `DistVoxToConsider` how far (in voxels and in every x and y direction)
#' does light diffuse horizontally?
#' * `calcSurfaceArea` TRUE/FALSE, should available surface area be
#' calculated?
#' * `calcSurfaceAreaLoss` TRUE/FALSE, should loss of surface are
#' (between timesteps) be calculated?
#' * `calcLightConditions` TRUE/FALSE, should light intensity in each voxel
#' be calculated?
#' @param shoot_dt a `data.frame` with branch information, with one row per
#' branch segment and the following columns:
#' * `xbegin` x-coordinate of the start of the segment
#' * `ybegin` y-coordinate of the start of the segment
#' * `zbegin` z-coordinate of the start of the segment
#' * `xend` x-coordinate of the end of the segment
#' * `yend` y-coordinate of the end of the segment
#' * `zend` z-coordinate of the end of the segment
#' * `length` length of the branch segment
#' * `diameter` diameter of the branch segment
#' * `shootID` unique identifier for this branch segment
#'
#' @param trunk_dt `data.frame` containing trunk information, with one row per
#' tree and the following columns:
#' * `x` the x coordinate of the tree trunk
#' * `y` the y coordinate of the tree trunk
#' * `height` height of the tree trunk
#' * `diameter` diameter of the tree trunk
#' * `treeID` unique identifier for this tree.
#'
#' @param vox_dt only required if `calcLightConditionsOpt = TRUE`,
#' a `data.frame` specifying the total leaf area in each voxel, with the
#' following columns:
#' * `x` x-coordinate of the voxel
#' * `y` y-coordinate of the voxel
#' * `z` z-coordinate of the voxel
#' * `leafarea` leaf area in this voxel
#'
#' @param path_to_output string, where to save output? Must be an `.rds` file
#' or `NULL`, in which case the result matrix is returned.
#' @param dead_branches_id integer vector containing the IDs of all branches
#' that will die this timestep
#' @param dead_trees_id integer vector containing the IDs of all trees that
#' will die this timestep
#'
#' @export
#'
create_microhabitat_mat <- function(config, shoot_dt, trunk_dt, vox_dt = NULL,
                                    path_to_output = NULL, dead_branches_id = NULL,
                                    dead_trees_id = NULL) {
  # Inputs are correct
  # check_config(config)
  # DistVoxToConsider <= corridor
  if (is.null(config$Imax))
    stop("Element Imax is missing from config list.")

  if (is.character(shoot_dt))
    utils::read.table(shoot_dt, sep = "\t",  header = TRUE, skip = 1)
  check_shoot_dt(shoot_dt)

  if (is.character(trunk_dt))
    utils::read.table(trunk_dt, sep = "\t",  header = TRUE, skip = 8)
  check_trunk_dt(trunk_dt)

  if (config$calcLightConditions) {
    if (is.character(vox_dt))
      utils::read.table(vox_dt, sep = "\t",  header = TRUE, skip = 1)
   check_vox_dt(vox_dt)
  }

  if (!is.null(path_to_output)) {
    dir_output <- dirname(path_to_output)
    if (!dir.exists(dirname(dir_output))) {
      stop(paste0("Output directory ", dir_output, " does not exist."))
    }
    if (!grepl("*.rds$", path_to_output)) {
      stop("path_to_output must be a rds file")
    }
  }

  # Set dimensions
  MaxX <- config$MaxX
  MaxY <- config$MaxY
  MaxZ <- config$MaxZ
  corridor <- config$corridor
  dimPlot <- c(MaxX, MaxY, MaxZ)
  forest_max_x <- MaxX + 2 * corridor
  forest_max_y <- MaxY + 2 * corridor
  voxel_area <- 100^2

  microhab_mat <- array(
    rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3] * 3),
    dim = c(dimPlot[1], dimPlot[2], dimPlot[3], 3)
  )

  # Element indices of the matrix
  sa_elt <- 1
  sa_loss_elt <- 2
  light_elt <- 3

  pb <- progress::progress_bar$new(
    format = "  Surface area for branches [:bar] :percent in :elapsed",
    total = nrow(shoot_dt)
    )
  pb$tick(0)

  for (s in seq_len(nrow(shoot_dt))) {

    seg_len <- shoot_dt$length[s]
    seg_diam <- shoot_dt$diameter[s]
    seg_start <- c(shoot_dt$xbegin[s] - corridor,
                   shoot_dt$ybegin[s] - corridor,
                   shoot_dt$zbegin[s])
    seg_end <- c(shoot_dt$xend[s] - corridor,
                 shoot_dt$yend[s] - corridor,
                 shoot_dt$zend[s])
    intersectd_voxels <- find_intersecting_voxels(seg_start, seg_end)

    # Calculate total surface area and split it evenly across intersected voxels
    seg_surface_area <- seg_len * seg_diam * pi / 2 / length(intersectd_voxels)

    # Drop voxels in the corridor
    is_in_microhab_area <- sapply(
      intersectd_voxels,
      function(v) v[1] <= MaxX && v[2] <= MaxY &&
        v[1] > 0 && v[2] > 0
      )
    intersectd_voxels <- intersectd_voxels[is_in_microhab_area]

    for (v in intersectd_voxels) {
      x <- v[1]
      y <- v[2]
      z <- v[3]
      voxel <- microhab_mat[x, y, z, ]

      if (config$calcSurfaceArea) {
        microhab_mat[x, y, z, sa_elt] <- voxel[sa_elt] + seg_surface_area
      }

      if (config$calcSurfaceAreaLoss && shoot_dt$shootID[s] %in% dead_branches_id) {
        microhab_mat[x, y, z, sa_loss_elt] <- voxel[sa_loss_elt] +
        seg_surface_area
      }
    }
    pb$tick()
  }

  pb <- progress::progress_bar$new(
    format = "  Surface area for trunks [:bar] :percent in :elapsed",
    total = nrow(trunk_dt)
  )
  pb$tick(0)

  for (t in seq_len(nrow(trunk_dt))) {

    x <- ceiling(trunk_dt$x[t]) - corridor
    y <- ceiling(trunk_dt$y[t]) - corridor
    trunk_height <- trunk_dt$height[t]
    trunk_diameter <- trunk_dt$diameter[t]

    SurfaceAreaTotal <- 0

    z_seq <- rev(seq_len(ceiling(trunk_height))) # top to bottom
    for (z in z_seq) {

      cone_height <- trunk_height - z + 1  # height of cylinder from top to bottom of voxel
      cone_radius <- trunk_diameter / 2  # radius of cylinder at bottom of voxel

      # Calculate total surface area in voxel for trunks
      SurfaceAreaInVoxel <- pi * cone_radius *
        sqrt(cone_radius^2 + cone_height^2) - SurfaceAreaTotal

      # Update total surface area of cylinder so far (to use in next step)
      SurfaceAreaTotal <- SurfaceAreaTotal + SurfaceAreaInVoxel

      if (config$calcSurfaceArea) {
        microhab_mat[x, y, z, sa_elt] <- microhab_mat[x, y, z, sa_elt] +
          SurfaceAreaInVoxel
      }

      # If trunk is lost during this time step, add it to lost surface
      if (trunk_dt$treeID[s] %in% dead_trees_id) {
        microhab_mat[x, y, z, sa_loss_elt] <- microhab_mat[x, y, z, sa_loss_elt] +
          SurfaceAreaInVoxel
      }
    } # z in z seq

    pb$tick()
  } # t in trunk set

  # Calculate light conditions in voxels (relative light conditions)
  if (config$calcLightConditions) {

    # TODO: parallelise and optimise this section

    light_range <- config$DistVoxToConsider

    # Total leaf area in each column
    # Must process voxels in the corridor too as they affect neighbouring voxels
    leaf_area_mat <- light_mat <- array(
      rep(0, forest_max_x * forest_max_y * MaxZ),
      dim = c(forest_max_x, forest_max_y, MaxZ)
      )

    # Store information on leaf area in matrix
    for (vx in seq_len(nrow(vox_dt))) {
      x <- vox_dt$x[vx]
      y <- vox_dt$y[vx]
      z <- vox_dt$z[vx]
      leaf_area_mat[x, y, z] <- vox_dt$leafarea[vx]
    }

    # Calculate single column light conditions based on leaf area distribution
    for (x in seq_len(forest_max_x)) {
      for (y in seq_len(forest_max_y)) {
        for (z in seq_len(MaxZ)) {
          total_leaf_area <- sum(leaf_area_mat[x, y, z:MaxZ])
          light_mat[x, y, z] <- exp(-config$kL * total_leaf_area / voxel_area)
        }
      }
    }

    # Calculate final light conditions by accounting for the light
    # conditions in adjacent voxels
    # x and y are indices in the full matrix including corridors
    x_seq <- seq(from = corridor + 1, to = forest_max_x - corridor)
    y_seq <- seq(from = corridor + 1, to = forest_max_y - corridor)
    pb <- progress::progress_bar$new(
      format = "  Calculating light conditions [:bar] :percent in :elapsed",
      total = length(x_seq) * length(y_seq)
    )
    pb$tick(0)
    for (x in x_seq) {
      for (y in y_seq) {
        for (z in seq_len(MaxZ)) {
          total_contribtn <- 0

          # loop over ring surrounding the focal voxel
          xx_seq <- seq(from = x - light_range, to = x + light_range)
          yy_seq <- seq(from = y - light_range, to = y + light_range)
          for (xx in xx_seq) {
            for (yy in yy_seq) {
              ring_index <- max(abs(xx - x), abs(yy - y)) # TODO: optimise this
              rel_contribtn <- 1 / (light_range + 1) / max(1, (ring_index * 8)) *
                light_mat[xx, yy, z] # TODO: optimise this
              total_contribtn <- total_contribtn + rel_contribtn
            }
          }
          microhab_mat[x - corridor, y - corridor, z, light_elt] <- total_contribtn
        } # z
        pb$tick()
      } # y
    } # x

    microhab_mat[,,,3] <- microhab_mat[,,,3] * config$Imax

  } # lightConditions

  if (!is.null(path_to_output)) {
    saveRDS(microhab_mat, path_to_output)
  } else {
    return(microhab_mat)
  }
}

check_shoot_dt <- function(shoot_dt) {

  exptd_cols <- c("xbegin", "ybegin", "zbegin", "xend", "yend", "zend", "length",
                  "diameter", "shootID")

  missing_params <- exptd_cols[!exptd_cols %in% names(shoot_dt)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("shoot_dt", missing_params))
  }

  col_names <- names(shoot_dt)

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(shoot_dt, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- col_names[is_missing_val]
    stop(paste(c(
      "The following elements of shoot_dt contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_numeric <- sapply(shoot_dt, function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- col_names[!is_numeric]
    stop(paste(c(
      "The following elements of shoot_dt contain non-numeric values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_positive <- sapply(shoot_dt[col_names], function(x) all(x >= 0))
  if (any(!is_positive)) {
    wrong_params <- col_names[!is_positive]
    stop(paste(c(
      "The following elements of shoot_dt contain negative values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (any(duplicated(shoot_dt$shootID))) {
    stop("shoot_dt contains multiple entries for the same branch (same shootID)")
  }
}

check_trunk_dt <- function(trunk_dt) {

  exptd_cols <- c("x", "y", "height", "diameter", "treeID")

  missing_params <- exptd_cols[!exptd_cols %in% names(trunk_dt)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("trunk_dt", missing_params))
  }

  col_names <- names(trunk_dt)

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(trunk_dt, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- col_names[is_missing_val]
    stop(paste(c(
      "The following elements of trunk_dt contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_numeric <- sapply(trunk_dt, function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- col_names[!is_numeric]
    stop(paste(c(
      "The following elements of trunk_dt contain non-numeric values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_positive <- sapply(trunk_dt[col_names], function(x) all(x >= 0))
  if (any(!is_positive)) {
    wrong_params <- col_names[!is_positive]
    stop(paste(c(
      "The following elements of trunk_dt contain negative values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (any(duplicated(trunk_dt$treeID))) {
    stop("trunk_dt contains multiple entries for the same tree (same treeID)")
  }
}

check_vox_dt <- function(vox_dt) {

  exptd_cols <- c("x", "y", "z", "leafarea")

  missing_params <- exptd_cols[!exptd_cols %in% names(vox_dt)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("vox_dt", missing_params))
  }

  col_names <- names(vox_dt)

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(vox_dt, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- col_names[is_missing_val]
    stop(paste(c(
      "The following elements of vox_dt contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_numeric <- sapply(vox_dt, function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- col_names[!is_numeric]
    stop(paste(c(
      "The following elements of vox_dt contain non-numeric values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  is_positive <- sapply(vox_dt[col_names], function(x) all(x >= 0))
  if (any(!is_positive)) {
    wrong_params <- col_names[!is_positive]
    stop(paste(c(
      "The following elements of vox_dt contain negative values:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }
}


