check_microhabitat <- function(microhab_mat) {
   dimensions <- dim(microhab_mat)
   if (length(dimensions) != 4) {
     stop("Microhabitat matrix must be a 4D array (3D matrix + (surface area, surface area loss, light and microclimate layers))")
   }

   if (!all(dimensions[1:3] >= 1)) {
     stop("Microhabitat matrix has a dimension of size 0.")
   }

   if (any(is.na(microhab_mat))) {
     stop("Microhabitat matrix contains NAs.")
   }

   if (any(is.null(microhab_mat))) {
     stop("Microhabitat matrix contains NULL elements.")
   }

   if(any(is.nan(microhab_mat))) {
     stop("Microhabitat matrix contains NaNs.")
   }

   if (any(microhab_mat < 0)) {
     stop("Microhabitat matrix contains negative elements.")
   }

   layer_map <- attr(microhab_mat, "layer_mapping")
   if (is.null(layer_map)) {
     stop("Microhabitat matrix must have an attribute named layer_mapping indicating the names (and position) of environmental layers.")
   }
   if (length(layer_map) != dimensions[4]) {
     stop("Layer map should have one element for each layer of the Microhabitat matrix.")
   }

   reqd_layers <- c(
     "surface_area", "surface_area_loss", "light"
   )
   missing_layers <- reqd_layers[!reqd_layers %in% layer_map]
   if (length(missing_layers > 0)) {

     stop(paste(c("Microhabitat matrix must have the following layers: ",
                  missing_layers), rep(" ", length(missing_layers) + 1)))
   }


}
