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
    sa_needed <- mass_to_surf_area(E$Mass[i], surface_biomass_scaling)
    avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] <- max(
      0, avail_sa_matrix[E$X[i], E$Y[i], E$Z[i]] - sa_needed
    )
  }
  return(avail_sa_matrix)
}

#' Get the surface area requirements of an epiphyte from its body mass
#'
#' \deqn{SA = \frac{M^{2/3}}{g_S}}
#'
#' @param mass mass of the individual (\eqn{M})
#' @param SurfaceBiomassScaling a scaling factor \eqn{g_S}
#' @export
mass_to_surf_area <- function(mass, SurfaceBiomassScaling) {
  return(mass^(2/3) / SurfaceBiomassScaling)
}
