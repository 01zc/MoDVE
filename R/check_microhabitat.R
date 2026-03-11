check_microhabitat <- function(microhab_mat) {
   dimensions <- dim(microhab_mat)
   if (length(dimensions) != 4) {
     stop("Microhabitat matrix must be a 4D array (3D matrix + (surface area, surface area loss, light))")
   }

   if (dimensions[4] != 3) {
     stop("Microhabitat matrix must be a 4D array (3D matrix + (surface area, surface area loss, light))")
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

}
