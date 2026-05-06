
#Microhabitat <- readRDS("tests/testthat/data/microhabitat/MicrohabitatMatrix1.rds")
#microclimate_mat <- readRDS("tests/testthat/data/microclimate/MicroclimateMatrix1.rds")[,,,c(1, 7, 11)]
#path_to_output

#' Expand microhabitat matrix with microclimate layers
#'
#' Add temperature, humidity and wind layers to the microhabitat matrix, from
#' a pre-assembled matrix containing this microclimate information
#'
#' @param Microhabitat a 4D matrix containing information about surface area,
#' surface area loss and light conditions in a 3D environment, assembled with
#' [create_microhabitat_mat()] or an equivalent alternative.
#' @param microclimate_mat a pre-assembled 4D array with the same X and Y
#' dimensions as `Microhabitat`, containing three layers of data corresponding
#' to the temperature, humidity and wind conditions in the voxels, respectively.
#' If the Z-dimension differs from `Microhabitat`, it is either clipped
#' (if larger), or the top layer is repeated to fill the gap (if smaller).
#' @param path_to_output either a path ending in `.rds` indicating where to save
#' the combined matrix, or `NULL`, in which case the output is returned.
#'
#' @export
#'
add_microclimate_layers <- function(Microhabitat, microclimate_mat, path_to_output = NULL) {

  dir_output <- dirname(path_to_output)
  if (!dir.exists(dir_output)) {
    stop(paste0("Output directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.rds$", path_to_output)) {
    stop("path_to_output must be a rds file")
  }

  # Load Microhabitat if not needed
  if (!is.array(Microhabitat)) {
    if (!grepl("*.rds$", Microhabitat)) {
      stop("Microhabitat should be an array or a valid path to a .rds file.")
    } else if (!file.exists(Microhabitat)) {
      stop(paste0(Microhabitat[i], " doesn't exist.\n"))
    } else {
      Microhabitat <- readRDS(Microhabitat)
    }
  }

  # Load Microclimate if needed
  if (!is.array(microclimate_mat)) {
    if (!grepl("*.rds$", microclimate_mat)) {
      stop("Microhabitat should be an array or a valid path to a .rds file.")
    } else if (!file.exists(microclimate_mat)) {
      stop(paste0(microclimate_mat[i], " doesn't exist.\n"))
    } else {
      microclimate_mat <- readRDS(microclimate_mat)
    }
  }

  dims_mcc <- dim(microclimate_mat)
  dims_mhb <- dim(Microhabitat)

  if (dims_mcc[4] != 3) {
    stop("microclimate_mat must contain 3 layers of data (temperature, humidity, and wind).")
  }

  if (!all(dim(microclimate_mat)[1:2] == dims_mhb[1:2])) {
    stop("Microhabitat and Microclimate matrices must have the same dimensions.")
  }
  microclimate_mat <- microclimate_mat[1:dims_mhb[1], 1:dims_mhb[2],,]

  # Initialise updated matrix
  dims_new <- dims_mhb
  dims_new[4] <- dims_new[4] + 3
  new_microhab_mat <- array(rep(as.numeric(NA), prod(dims_new)), dim = dims_new)
  new_microhab_mat[,,,1:dims_mhb[4]] <- Microhabitat
  rm(Microhabitat)

  mcc_indices <- (dims_mhb[4] + 1):(dims_mhb[4] + dims_mcc[4])

  # If Z does not match, issue warning and populate with top layer
  z_diff <- dims_mhb[3] - dims_mcc[3]
  if (z_diff < 0) {
    warning("Microclimate matrix Z-dimension is larger than Microhabitat. Clipping for assembly.")
    new_microhab_mat[,,1:dims_mcc[3], mcc_indices] <- microclimate_mat[,,1:dims_mcc[3],]
  } else if (z_diff > 0) {

    warning(paste0(
      "Microclimate matrix Z-dimension is smaller than Microhabitat.
Top layer is duplicated ", z_diff, " times to fill the gap."
    ))
    new_microhab_mat[,,1:dims_mcc[3], mcc_indices] <- microclimate_mat[,,1:dims_mcc[3],]

    top_layer <- microclimate_mat[,,dims_mcc[3],]
    z_range <- (dims_mcc[3] + 1):(dims_mcc[3] + z_diff)
    new_microhab_mat[,,z_range, mcc_indices] <- top_layer
  }

  if (!is.null(path_to_output)) {
    saveRDS(new_microhab_mat, path_to_output)
  } else {
    return(new_microhab_mat)
  }
}
