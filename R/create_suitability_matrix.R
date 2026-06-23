#' Create suitability matrix
#'
#' Pre-compute global environmental suitability scores between 0 and 1 for each species.
#'
#' @param Microhabitat a matrix containing the surface area, light and
#' optionally microclimatic conditions in each voxel, as created by
#' [create_microhabitat_mat()]
#' @param SpeciesPool a data frame containing the species traits, as created by
#' [draw_species_traits()].
#' @param use_parabolic_light boolean, if `TRUE` (default), the light suitability score is
#' calculated using the parabolic response curve instead of the cardinal score.
#' @param path_to_output name and path of a `.rds` or `.h5` file to save the suitability
#' matrix in. If `NULL` the matrix is returned instead.
#'
#' @details
#' Cardinal suitability scores (see [calc_suitability()]) are first calculated
#' for each environmental (light, temperature, humidity, wind) layer present in
#' `Microhabitat`. Global suitability scores are then computed as the product of
#' these scores for each species.
#'
#' @return if `path_to_output == NULL`, returns a 4 dimension array with dimensions
#' x, y, z, number of species, where each layer contains the aggregated suitability scores
#' for each species.
#'
#' @export
create_suitability_matrix <- function(
    Microhabitat,
    SpeciesPool,
    use_parabolic_light = TRUE,
    path_to_output = NULL
) {

  check_microhabitat(Microhabitat)
  check_species_df(SpeciesPool, get_microclimate_opts(Microhabitat))

  # Prepare output
  if (!is.null(path_to_output)) {
    dir_output <- dirname(path_to_output)
    if (!dir.exists(dir_output)) {
      stop(paste0("Output directory ", dir_output, " does not exist."))
    }
    file_ext <- strsplit(path_to_output, "\\.", fixed = FALSE)[[1]][2]
    if (!file_ext %in% c(".rds", ".h5")) {
      stop("path_to_output must be a .rds or .h5 file")
    }
  }

  layer_map <- attr(Microhabitat, "layer_mapping")
  keys <- c("light" = "Light", "humidity" = "Hum", "temperature" = "Temp", "wind" = "Wind")
  active_vars <- keys[names(keys) %in% layer_map]

  nb_sp <- nrow(SpeciesPool)
  suitability_mat <- array(1, dim = c(dim(Microhabitat)[1:3], nb_sp))

  for (sp in seq_len(nb_sp)) {

    for (i in seq_len(dim(Microhabitat)[4])) {

      layer_name <- layer_map[i]
      var_name <- keys[names(keys) == layer_name]
      if (length(var_name) == 0) next # not a microclimate variable

      var <- Microhabitat[,,, i]

      if (var_name == "Light" & use_parabolic_light) {
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

      # Combine scores (product) across variables
      suitability_mat[,,,sp] <- suitability_mat[,,,sp] * suitability_var
    }
  } # for each species

  if (!is.null(path_to_output)) {
    if (file.exists(path_to_output)) file.remove(path_to_output)
    if (file_ext == "h5") {
    rhdf5::h5createFile(path_to_output)
    rhdf5::h5write(suitability_mat, path_to_output, "EnvironmentalSuitabilityScores")
    } else { # rds
      saveRDS(suitability_mat, path_to_output)
    }
  } else {
    return(suitability_mat)
  }
}
