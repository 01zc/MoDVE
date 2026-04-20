#' Resolve the competition step of the epiphyte simulation
#'
#' Each year, the capacity (surface area) of each voxel is compared to its
#' occupancy (sum of surface area requirements of all epiphytes in the voxel).
#' If occupancy exceeds capacity, either random (if `massDepCompetiton = FALSE`)
#' or the smallest epiphytes (if `massDepCompetiton = TRUE`) are killed (i.e,
#' marked with status 2) until the voxel is no longer oversaturated.
#'
#' @param E a `data.frame` containing the individual epiphytes present in the
#' landscape.
#' @param Microhabitat the microhabitat matrix, containing surface area, loss
#' and light conditions.
#' @param massDepCompetition, `TRUE` = larger individuals get priority in
#' voxel attribution, otherwise (`FALSE`) individuals are distributed randomly.
#'
#' @returns the modified epiphyte data frame
#' @export
#'
resolve_competition <- function(E, Microhabitat, massDepCompetition) {

  # Calculate total surface area occupied by epiphytes per voxel
  dims <- dim(Microhabitat)
  occupiedSA <- array(
    rep(0, dims[1] * dims[2] * dims[3]),
    dim = c(dims[1], dims[2], dims[3])
  )
  for (w in seq_len(nrow(E))) {
    if (E$Status[w] == 1) {
      occupiedSA[E$X[w], E$Y[w], E$Z[w]] <-
        occupiedSA[E$X[w], E$Y[w], E$Z[w]] +
        E$SurfaceAreaOccupied[w]
      # for each cell, sum sa over all inds in this cell
    }
  }

  # Oversaturated voxels: occupied S.A. exceeds available S.A.
  oversat_voxels <- arrayInd(
    which(occupiedSA > Microhabitat[, , , 1]),
    dim(occupiedSA)
  )
  vox_xs <- oversat_voxels[, 1]
  vox_ys <- oversat_voxels[, 2]
  vox_zs <- oversat_voxels[, 3]

  for (i in seq_len(length(vox_xs))) {

    # Get all individuals in this voxel
    isInVoxel <- E$X == vox_xs[i] & E$Y == vox_ys[i] & E$Z == vox_zs[i]
    indsInVoxel <- E[isInVoxel & E$Status == 1, ]

    # Sort individuals
    if (massDepCompetition) { # priority to larger individuals
      ind_seq <- order(indsInVoxel$SurfaceAreaOccupied, decreasing = TRUE)
    } else { # random
      ind_seq <- sample(seq_len(nrow(indsInVoxel)))
    }
    indsInVoxel <- indsInVoxel[ind_seq, ]

    totalSurfAreaOcc <- cumsum(indsInVoxel$SurfaceAreaOccupied)
    availSurfArea <- Microhabitat[vox_xs[i], vox_ys[i], vox_zs[i], 1]

    # The voxel can support this many individuals
    capacity <- length(which(totalSurfAreaOcc <= availSurfArea))

    if (capacity < nrow(indsInVoxel)) {
      # largest n individuals live, the rest die
      seq_beyond_capacity <- int_seq(capacity + 1, nrow(indsInVoxel))
      dead_ids <- indsInVoxel[seq_beyond_capacity, "IndividualID"]
      these_die <- E$IndividualID %in% dead_ids$IndividualID
      E[these_die, "Status"] <- 2
    }
  }
  return(E)
}
