#' Resolve the growth step of the simulation
#'
#' Each year, all living epiphytes grow by an amount corresponding to
#' \eqn{k * (M_{max} - M) * I}, where *k* is the growth rate, *M* the epiphyte's
#' mass, \eqn{M_{max}} its species maximum mass, and I is the parabolic light
#' response curve:
#' \deqn{I = I_A * (I^{XYZ})^2 + I_B * I^{XYZ} + I_C} where \eqn{I_A, I_B, I_C}
#' are species parameters derived from the species' light niche, and
#' \deqn{I^{XYZ}} is the current light intensity in the voxel.
#'
#' @param E a `data.frame` containing the individual epiphytes present in the
#' landscape
#' @param SpeciesPool a `data.frame` containing species-level traits
#' @param Microhabitat the microhabitat matrix, containing surface area, loss
#' and light conditions.
#' @param SurfaceBiomassScaling numeric parameter determining the surface area
#' occupied by an epiphyte as a function of its mass:
#' \deqn{M^{2/3} / g_S}
#'
#' @returns the modified epiphyte data frame
#' @export
#'
resolve_growth <- function(E, SpeciesPool, Microhabitat, SuitabilityMat, SurfaceBiomassScaling) {

  layer_map <- attr(Microhabitat, "layer_mapping")
  microclimate_opts <- get_microclimate_opts(Microhabitat)
  idx_sa <- which(layer_map == "surface_area")
  idx_sa_loss <- which(layer_map == "surface_area_loss")
  idx_light <- which(layer_map == "light")
  idx_temperature <- which(layer_map == "temperature")
  idx_wind <- which(layer_map == "wind")
  idx_humidity <- which(layer_map == "humidity")

  for (i in seq_len(nrow(E))) {

    vox <- Microhabitat[E$X[i], E$Y[i], E$Z[i], ]

    # TODO: could be faster without if statement => speed testing
    if (E$Status[i] == 1) {

      sp_row <- which(SpeciesPool$SpeciesID == E$SpeciesID[i])

      # Von Bertalanffy growth function
      growth_term <- SpeciesPool$GrowthRate[sp_row] *
        (SpeciesPool$MaximumMass[sp_row] - E$Mass[i])

      env_term <- SuitabilityMat[E$X[i], E$Y[i], E$Z[i], E$SpeciesID[i]]
      if (is.na(env_term) | is.nan(env_term)) {env_term <- 0}
      all_suits_prec <- c(all_suits_prec, env_term)
      E$Mass[i] <- E$Mass[i] + max(0, growth_term * env_term)
    }

    # Add info about the voxel to the epiphyte matrix
    E$SurfaceAreaOccupied[i] <- mass_to_surf_area(E$Mass[i], SurfaceBiomassScaling)
    # TODO: do we really need individual-level copies of these habitat values?
    # This is only for output, not used during simulation
    E$TotalSurfaceInVoxel[i] <- vox[idx_sa]  # Total surface in voxel
    E$SurfaceLossInVoxel[i] <- vox[idx_sa_loss]  # Percentage surface loss in this year
    E$LightInVoxel[i] <- vox[idx_light]  # Light conditions in voxel
    if (microclimate_opts$use_humidity) E$HumInVoxel[i] <- vox[idx_humidity]
    if (microclimate_opts$use_temperature) E$TempInVoxel[i] <- vox[idx_temperature]
    if (microclimate_opts$use_wind) E$WindInVoxel[i] <- vox[idx_wind]
  }

  return(E)
}
