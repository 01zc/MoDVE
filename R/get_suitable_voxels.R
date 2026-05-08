get_suitable_voxels <- function(Microhabitat,
                                  microclimate_opts,
                                  SpeciesPool,
                                  which_species) {

  layer_map <- attr(microhabitat, "layer_mapping")

  # Base condition: Total surface area option must be present and > 0
  SuitableMask <- Microhabitat[,,,layer_map["surface_area"]] > 0

  # ---- Light ----
  if (microclimate_opts$use_light) {
    LightIdx <- layer_map["light"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , LightIdx] >= SpeciesPool$MinLight[which_species] &
      Microhabitat[, , , LightIdx] <= SpeciesPool$MaxLight[which_species]
  }

  # ---- Humidity ----
  if (microclimate_opts$use_humidity) {
    HumIdx <- layer_map["humidity"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , HumIdx] >= SpeciesPool$MinHum[which_species] &
      Microhabitat[, , , HumIdx] <= SpeciesPool$MaxHum[which_species]
  }

  # ---- Temperature ----
  if (microclimate_opts$use_light$use_temperature) {
    TempIdx <- layer_map["temperature"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , TempIdx] >= SpeciesPool$MinTemp[which_species] &
      Microhabitat[, , , TempIdx] <= SpeciesPool$MaxTemp[which_species]
  }

  # ---- Wind ----
  if (microclimate_opts$use_wind) {
    WindIdx <- layer_map["wind"]
    SuitableMask <- SuitableMask &
      Microhabitat[, , , WindIdx] >= SpeciesPool$MinWind[which_species] &
      Microhabitat[, , , WindIdx] <= SpeciesPool$MaxWind[which_species]
  }

  return(which(SuitableMask))
}
