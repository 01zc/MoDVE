#' Find which voxels are suitable for a given species
#'
#' Given a microhabitat matrix and a set of environmental variables,
#' return the set of voxels which microclimatic conditions are suitable for this
#' species.
#'
#' @param Microhabitat a 3-D matrix with one layer per environmental variable`.
#' @param species_row a single row of a species data frame containing this
#' species attributes, including climatic niche.
#'
#' @export
#'
get_suitable_voxels <- function(Microhabitat, species_row) {

  layer_map <- attr(Microhabitat, "layer_mapping")
  microclimate_opts <- list(
    "use_light" = "light" %in% layer_map,
    "use_temperature" = "temperature" %in% layer_map,
    "use_humidity" = "humidity" %in% layer_map,
    "use_wind" = "wind" %in% layer_map,
  )

  # Base condition: Total surface area option must be present and > 0
  SuitableMask <- Microhabitat[,,,layer_map["surface_area"]] > 0

  # ---- Light ----
  if (microclimate_opts$use_light) {
    LightIdx <- layer_map["light"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , LightIdx] >= species_row$MinLight &
      Microhabitat[, , , LightIdx] <= species_row$MaxLight
  }

  # ---- Humidity ----
  if (microclimate_opts$use_humidity) {
    HumIdx <- layer_map["humidity"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , HumIdx] >= species_row$MinHum &
      Microhabitat[, , , HumIdx] <= species_row$MaxHum
  }

  # ---- Temperature ----
  if (microclimate_opts$use_temperature) {
    TempIdx <- layer_map["temperature"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , TempIdx] >= species_row$MinTemp &
      Microhabitat[, , , TempIdx] <= species_row$MaxTemp
  }

  # ---- Wind ----
  if (microclimate_opts$use_wind) {
    WindIdx <- layer_map["wind"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , WindIdx] >= species_row$MinWind &
      Microhabitat[, , , WindIdx] <= species_row$MaxWind
  }

  return(which(SuitableMask))
}
