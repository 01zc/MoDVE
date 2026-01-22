#' Title
#'
#' @param E epiphyte data frame
#' @param SpeciesPool species data frame
#' @param Microhabitat microhabitat matrix
#' @param SurfaceBiomassScaling numeric parameter
#'
#' @returns the modified epiphyte data frame
#' @export
#'
resolve_growth <- function(E, SpeciesPool, Microhabitat, SurfaceBiomassScaling) {

  for (i in seq_len(nrow(E))) {

    vox <- Microhabitat[E$X[i], E$Y[i], E$Z[i], ]

    # TODO: could be faster without if statement => speed testing
    if (E$Status[i] == 1) {

      sp_row <- which(SpeciesPool$SpeciesID == E$SpeciesID[i])

      # Von Bertalanffy growth function
      growth_term <- SpeciesPool$GrowthRate[sp_row] *
        (SpeciesPool$MaximumMass[sp_row] - E$Mass[i])

      # Parabolic light response
      light_vox <- vox[3]
      light_term <- max(0, SpeciesPool$LightResponseA[sp_row] * light_vox^2 +
        SpeciesPool$LightResponseB[sp_row] * light_vox +
        SpeciesPool$LightResponseC[sp_row])

      E$Mass[i] <- E$Mass[i] + growth_term * light_term
    }

    # Add info about the voxel to the epiphyte matrix
    E$SurfaceAreaOccupied[i] <- (E$Mass[i]^(2/3)) / SurfaceBiomassScaling
    # TODO: do we really need individual-level copies of these habitat values?
    E$TotalSurfaceInVoxel[i] <- vox[1]  # Total surface in voxel
    E$SurfaceLossInVoxel[i] <- vox[2]  # Percentage surface loss in this year
    E$LightInVoxel[i] <- vox[3]  # Light conditions in voxel
  }

  return(E)
}
