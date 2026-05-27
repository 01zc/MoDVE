#' Resolve the mortality step of the epiphyte simulation
#'
#' Each year, epiphytes may die from either of the following sources:
#' * Branch fall: each epiphytes dies (status 3) with a probability equal to the
#'  proportion of surface area loss (branch fall) in its voxel this generation
#' * Light conditions: if light in the voxels falls outside of the epiphyte's
#' light niche, it dies (status 4).
#' * Base mortality (status 5): each epiphyte dies with probability
#' `MortRateRandom`, or if `use_mass_dep_mortality == TRUE`, with probability
#' `MortRateMass` * mass ^`MortRateMassScaling.`
#'
#' @param E a `data.frame` containing the individual epiphytes present in the
#' landscape
#' @param SpeciesPool a `data.frame` containing species-level traits
#' @param Microhabitat the microhabitat matrix, containing surface area, loss
#' and light conditions.
#' @param use_mass_dep_mortality boolean, if `FALSE` individuals die randomly
#' according to `MortRateRandom`, if `TRUE` mortality is mass-dependent, using
#' `MortRateMass * (mass^MortRateMassScaling)`.
#' @param MortRateRandom numeric between 0 and 1, the probability of an
#' individual dying if `use_mass_dep_mortality == FALSE`
#' @param MortRateMass numeric, if `use_mass_dep_mortality == TRUE` the
#' coefficient for the effect of mass on the probability of death
#' @param MortRateMassScaling numeric between `-Inf` and `0`, if
#' `use_mass_dep_mortality == TRUE` the exponent for the effect of mass on the
#' probability of death. Must be negative or zero.
#'
#' @returns the modified epiphyte data frame
#' @export
#'
resolve_mortality <- function(E, SpeciesPool, Microhabitat, use_mass_dep_mortality,
                              MortRateRandom, MortRateMass, MortRateMassScaling) {

  layer_map <- attr(Microhabitat, "layer_mapping")
  idx_sa_loss <- which(layer_map == "surface_area_loss")
  idx_light <- which(layer_map == "light")
  idx_temperature <- which(layer_map == "temperature")
  idx_humidity <- which(layer_map == "humidity")
  idx_wind <- which(layer_map == "wind")

  for (i in seq_len(nrow(E))) {
    if (E$Status[i] == 1) {

      vox <- Microhabitat[E$X[i], E$Y[i], E$Z[i],]
      this_species <- SpeciesPool$SpeciesID == E$SpeciesID[i]
      min_light <- SpeciesPool$MinLight[this_species]
      max_light <- SpeciesPool$MaxLight[this_species]
      min_hum <- SpeciesPool$MinHumidity[this_species]
      max_hum <- SpeciesPool$MaxHumidity[this_species]
      min_wind <- SpeciesPool$MinWind[this_species]
      max_wind <- SpeciesPool$MaxWind[this_species]
      min_temp <- SpeciesPool$MinTemperature[this_species]
      max_temp <- SpeciesPool$MaxTemperature[this_species]

      # The following comparison would fail without the is.nan check,
      # because Microhabitat contains NaNs in some entries and
      # in R a comparison with a NaN returns NA, not a boolean.
      # Note: We call runif repeatedly intentionally. See Issue #16 on Github
      # TODO: discuss priority among the different sources of mortality
      # e.g. the rate of random mortality won't match the parameter
      # because some fraction has already died from branch fall/light
      # Once priority is clarified drawing mortality from sa loss and
      # random/mass mortality could be vectorised

      # Branch fall mortality
      if (!is.nan(vox[idx_sa_loss]) &&
          stats::runif(1, min = 0, max = 1) < vox[idx_sa_loss]) {
        E$Status[i] <- 3
      } else if ("light" %in% layer_map &&
                 vox[idx_light] < min_light | vox[idx_light] > max_light) {
        # Unsuitable light conditions
        E$Status[i] <- 4
      } else if (!use_mass_dep_mortality &&
                 stats::runif(1, min = 0, max = 1) < MortRateRandom) {  # Natural mortality rate
        # Baseline random mortality
        E$Status[i] <- 5
      } else if (use_mass_dep_mortality &&
                 stats::runif(1, min = 0, max = 1) <
                 (MortRateMass * (E$Mass[i]^MortRateMassScaling))) {
        # Mass-dependent mortality
        E$Status[i] <- 5
      } else if ("humidity" %in% layer_map && !is.na(vox[idx_humidity]) &&
                 (vox[idx_humidity] < min_hum | vox[idx_humidity] > max_hum)) {
        E$Status[i] <- 6
      } else if ("temperature" %in% layer_map && !is.na(vox[idx_temperature]) &&
                 (vox[idx_temperature] < min_temp | vox[idx_temperature] > max_temp)) {
        E$Status[i] <- 7
      } else if ("wind" %in% layer_map && !is.na(vox[idx_wind]) &&
                 min_wind > max_wind) {
        E$Status[i] <- 8
      }
    }
  }
  return(E)
}
