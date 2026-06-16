#' Calculate the probabilistic distribution of dispersal in a 3D matrix
#'
#' Given a focal point in a 3D matrix, calculate the probability to disperse in
#' each voxel around the focal point.
#'
#' @details The 3D matrix is an expanded version of the microhabitat matrix used
#' in the main simulation, such that the dispersal matrix covers dispersal from
#' any given point in the microhabitat matrix.
#' During the dispersal step, this matrix is centered on the focal
#' individual and cropped to the edges of the microhabitat, giving the
#' probabilities of the offspring dispersing in each voxel.
#'
#' @param dims 3-element vector with the X, Y, Z dimensions of a microhabitat matrix
#' @param SpeciesPool a data frame containing the species traits of all species
#' @param WindSpeed the wind layer of the Microhabitat matrix, either a 3-D
#' matrix with dimensions X, Y, Z or a vector of equivalent length.
#'
#' @returns a matrix with the base probability of dispersing from the central
#' point to each cell within reach. The matrix sums to 1.
#' @export
#'

calc_prob_disp_matrix <- function(dims, SpeciesPool, WindSpeed = NULL) {
  use_wind_dispersal <- !is.null(WindSpeed)

  # Generate probabilities to disperse in any direction from any voxel
  # in the matrix
  # --> we need a matrix twice as large as microhabitat
  expanded_dims <- dims[1:3] * 2 + 1
  expanded_mat_central_point <- dims + 1

  # Calculate distance to central point
  DistanceMatrix <- array(
    rep(0, prod(expanded_dims)),
    dim = expanded_dims
  )
  for (i in seq_len(expanded_dims[1])) {
    for (j in seq_len(expanded_dims[2])) {
      for (k in seq_len(expanded_dims[3])) {
        x1 <- c(i, j, k)
        x2 <- expanded_mat_central_point
        DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))
      }
    }
  }

  expanded_wind_mat <- 0
  if (use_wind_dispersal) {
    # Modify the windspeed matrix to match the larger distance matrix
    expanded_wind_mat <- array(0, dim = expanded_dims)

    # Coordinates to insert original wind field in center
    original_dims <- (expanded_dims - 1) / 2
    start <- floor(expanded_dims / 2) - floor(original_dims / 2) + 1
    end <- start + original_dims - 1

    # Embed the original wind field into the center
    expanded_wind_mat[start[1]:end[1], start[2]:end[2], start[3]:end[3]] <- WindSpeed
  }

  # Get probabilities to disperse in each voxel
  NumberOfSpecies <- nrow(SpeciesPool)
  ProbabilityMatrix <- prob_disp_matrix <- array(
    rep(0, prod(expanded_dims) * NumberOfSpecies),
    dim = c(expanded_dims, NumberOfSpecies)
  )

  for (i in seq_len(NumberOfSpecies)) {

    # Negative exponential
    exponent <- DistanceMatrix * SpeciesPool$DispersalKernel[i]
    # Scale the dispersal kernel by wind speed (depending on species specific wind dispersal)
    if (use_wind_dispersal) {
      exponent <- exponent /
        (1 + SpeciesPool$DispersalKernelWindEffect[i] * expanded_wind_mat)
    }
    ProbabilityMatrix[, , , i] <- exp(-exponent)

    # Dispersal asymmetry (probability to disperse downwards > upwards)
    dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]
    z_seq_up <- int_seq(
      from = expanded_mat_central_point[3],
      to = expanded_dims[3],
      by = 1
      )
    z_seq_down <- int_seq(
      from = 1,
      to = expanded_mat_central_point[3] - 1,
      by = 1
      )
    ProbabilityMatrix[,,z_seq_up, i] <- ProbabilityMatrix[, , z_seq_up, i] *
      2 * (1 - dispersalAsymmetry)
    ProbabilityMatrix[,,z_seq_down, i] <- ProbabilityMatrix[, , z_seq_down, i] *
      2 * dispersalAsymmetry

    # Normalize
    prob_disp_matrix[, , , i] <- ProbabilityMatrix[, , , i] /
      sum(ProbabilityMatrix[, , , i])
  }
  return(prob_disp_matrix)
}
