#' Extract individual coordinates in vector form
#'
#' Reads a count matrix and extracts the x, y and z coordinates of all
#' individuals
#'
#' @param count_matrix a 3d integer matrix containing tallies of individuals
#' present in each cell
#'
#' @returns a list of 3 integer vectors `x`, `y` and `z` containing the
#' corresponding coordinates of all individuals
#' @export
#'
extract_ind_coords <- function(count_matrix) {

  nb_inds <- sum(count_matrix)

  ids <- arrayInd(
    which(count_matrix > 0),
    dim(count_matrix)
  )
  x_inds <- ids[, 1]
  y_inds <- ids[, 2]
  z_inds <- ids[, 3]

  while (nb_inds > length(x_inds)) {

    # recursively distribute coordinates until counts of remaining
    # unprocessed individuals reaches zero
    tmp_ids <- arrayInd(
      which(count_matrix > 0),
      dim(count_matrix)
    )

    # decrement count
    count_matrix[tmp_ids] = count_matrix[tmp_ids] - 1

    tmp_ids <- arrayInd(
      which(count_matrix > 0),
      dim(count_matrix)
    )

    x_inds <- append(x_inds, tmp_ids[, 1])
    y_inds <- append(y_inds, tmp_ids[, 2])
    z_inds <- append(z_inds, tmp_ids[, 3])
  }

  return(list(x  = x_inds, y = y_inds, z = z_inds))
}
