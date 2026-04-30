
climate_var_indices <- c(
  "temperature" = 1,
  "humidity" = 7,
  "wind" = 11
)

# Load Microclimate
microclimate_mat <- readRDS(
  "tests/testthat/data/microclimate/MicroclimateMatrix1.rds"
  )[,,, climate_var_indices] # Subset to variables of interest
dims_mcc <- dim(microclimate_mat)

Microhabitat <- readRDS("tests/testthat/data/microhabitat/MicrohabitatMatrix1.rds")
dims_mhb <- dim(Microhabitat)

if (!all(dim(microclimate_mat)[1:2] == dims_mhb[1:2])) {
  stop("Microhabitat and Microclimate matrices must have the same dimensions.")
}
microclimate_mat <- microclimate_mat[1:dims_mhb[1], 1:dims_mhb[2],,]

# Initialise updated matrix
dims_new <- dims_mhb
dims_new[4] <- dims_new[4] + length(climate_var_indices)
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


path_to_output <- ".rds"
saveRDS(new_microhab_mat, path_to_output)

# Merge
# Save
