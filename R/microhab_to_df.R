
#' Convert the microhabitat matrix into a data frame
#'
#' @param Microhabitat the matrix to convert
#' @param keep_empty_voxels if `FALSE`, only include voxels with positive
#' surface area
#'
#' @export
microhabitat_to_df <- function(Microhabitat, keep_empty_voxels = FALSE) {
  if (keep_empty_voxels) {
    which_voxels <- 1:length(Microhabitat)
  } else {
    which_voxels <- which(Microhabitat[,,,1] > 0)
  }
  microhab_tbl <- map_dfr(which_voxels, index_to_3d_coords, dimX = 80, dimY = 80)
  microhab_tbl$sa <- Microhabitat[,,,1][which_voxels]
  microhab_tbl$sa_loss <- Microhabitat[,,,2][which_voxels]
  microhab_tbl$light <- Microhabitat[,,,3][which_voxels]
}
