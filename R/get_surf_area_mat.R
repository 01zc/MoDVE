#' Get available surface area
#'
#' Compute a matrix containing the remaining available surface area for each
#' voxel after deducing for the area occupied by individuals.
#'
#' @param E a data frame containing the attributes of the individuals in the
#' MoDVE simulation
#' @param surface_area_mat a 3-D array containing the base surface area made
#' available thanks to trunks and branches
#' @param surface_biomass_scaling numeric parameter
#'
get_surf_area_mat <- function(E, surface_area_mat, surface_biomass_scaling) {
  avail_sa_matrix <- surface_area_mat
  for (i in seq_len(nrow(E))) {
    # Deduce the surface area occupied by each individual
    sa_needed <- E$Mass[i]^(2/3) / surface_biomass_scaling
    avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] <- max(
      0, avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] - sa_needed
    )
  }
  return(avail_sa_matrix)
}
