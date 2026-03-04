#include <Rcpp.h>
using namespace Rcpp;

int index_3d(const int& x, const int& y, const int& z, const int& rows, const int& cols) {
  return x * rows * cols + y * cols + z;
}

/*
// [[Rcpp::export]]
NumericVector calc_light_matrix(const DataFrame& vox_dt, const List& config) {

  const int light_range = config["DistVoxToConsider"];
  const double kL = config["kL"];
  const int max_x = config["MaxX"];
  const int max_y = config["MaxY"];
  const int max_z = config["MaxZ"];
  const int corridor = config["corridor"];
  const int max_x_corr = max_x + 2 * corridor;
  const int max_y_corr = max_y + 2 * corridor;

  const IntegerVector col_x = vox_dt["x"];
  const IntegerVector col_y = vox_dt["y"];
  const IntegerVector col_z = vox_dt["z"];
  const NumericVector col_leaf_area = vox_dt["leafarea"];

  NumericVector leaf_area_mat(Dimension(max_x_corr, max_y_corr, max_z));

  int x, y, z;
  double leaf_area;
  for (int i = 0; i < col_x.size(); i++) {
    // faster implementation possible with STL?
    x = col_x[i];
    y = col_y[i];
    z = col_z[i];
    leaf_area = col_leaf_area[i];
    leaf_area_mat[index_3d(x, y, z, max_x_corr, max_y_corr)] = leaf_area;
  }
  return leaf_area_mat;

# Total leaf area in each column
# Must process voxels in the corridor too as they affect neighbouring voxels
   leaf_area_mat <- light_mat <- array(
   rep(0, forest_max_x * forest_max_y * MaxZ),
   dim = c(forest_max_x, forest_max_y, MaxZ)
   )

# Store information on leaf area in matrix
   for (vx in seq_len(nrow(vox_dt))) {
   x <- vox_dt$x[vx]
   y <- vox_dt$y[vx]
   z <- vox_dt$z[vx]
   leaf_area_mat[x, y, z] <- vox_dt$leafarea[vx]
   }

# Calculate single column light conditions based on leaf area distribution
   for (x in seq_len(forest_max_x)) {
   for (y in seq_len(forest_max_y)) {
   for (z in seq_len(MaxZ)) {
   total_leaf_area <- sum(leaf_area_mat[x, y, z:MaxZ])
   light_mat[x, y, z] <- exp(-config$kL * total_leaf_area / voxel_area)
   }
   }
   }

# Calculate final light conditions by accounting for the light
# conditions in adjacent voxels
# x and y are indices in the full matrix including corridors
   for (x in seq(from = corridor + 1, to = forest_max_x - corridor)) {
   for (y in seq(from = corridor + 1, to = forest_max_y - corridor)) {
   for (z in seq_len(MaxZ)) {
   total_contribtn <- 0

# loop over ring surrounding the focal voxel
   xx_seq <- seq(from = x - light_range, to = x + light_range)
   yy_seq <- seq(from = y - light_range, to = y + light_range)
   for (xx in xx_seq) {
   for (yy in yy_seq) {
   ring_index <- max(abs(xx - x), abs(yy - y))
   rel_contribtn <- 1 / (light_range + 1) / max(1, (ring_index * 8)) *
   light_mat[xx, yy, z]
   total_contribtn <- total_contribtn + rel_contribtn
   }
   }
   microhab_mat[x - corridor, y - corridor, z, light_elt] <- total_contribtn
   } # z
   } # y
   } # x

}
*/


