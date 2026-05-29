
create_suitability_matrix <- function(Microhabitat, SpeciesPool, path_to_output = NULL, use_parabolic_light = TRUE) {

  layer_map <- attr(Microhabitat, "layer_mapping")
  keys <- c("light" = "Light", "humidity" = "Hum", "temperature" = "Temp", "wind" = "Wind")
  active_vars <- keys[names(keys) %in% layer_map]

  nb_sp <- nrow(SpeciesPool)
  suitability_mat <- array(1, dim = c(dim(Microhabitat[1:3], nb_sp)))

  for (sp in seq_len(nb_sp)) {

    for (i in seq_len(dim(Microhabitat)[4])) {

      layer_name <- layer_map[i]
      var_name <- keys[names(keys) == layer_name]
      if (length(var_name) == 0) next # not a microclimate variable

      var <- Microhabitat[,,,i]

      if (var_name == "Light" & use_parabolic_light == "Parabolic") {
        suitability_var <- get_parabolic_resp(
          var,
          SpeciesPool$LightResponseA[sp],
          SpeciesPool$LightResponseB[sp],
          SpeciesPool$LightResponseC[sp]
        )
      } else {
        suitability_var <- calc_suitability(
          var,
          SpeciesPool[sp, paste0("Min", var_name)],
          SpeciesPool[sp, paste0("Max", var_name)],
          SpeciesPool[sp, paste0("Optimum", var_name)]
        )
      }
    }

    # Combine scores (product) across variables
    suitability_mat[,,,sp] <- suitability_mat[,,,sp] * suitability_var

  } # for each species

  if (!is.null(path_to_output)) {
    if (file.exists(path_to_output)) file.remove(path_to_output)
    rhdf5::h5createFile(path_to_output)
    rhdf5::h5write(suitability_mat, path_to_output, "EnvironmentalSuitabilityScores")
  } else {
    return(suitability_mat)
  }
}
