#' Combine branch segments to return one row per branch
#'
#' @param shoots_tbl a table containing branch segments
#'
#'@export
combine_branch_segments <- function(shoots_tbl) {
  shoots_tbl |> dplyr::group_by(treeID, branchID, order) |>
    dplyr::summarise(
      "length" = sum(length),
      "diameter" = max(diameter),
      "xbegin" = ifelse(xbegin[1] < xend[1], min(xbegin), max(xbegin)),
      "ybegin" = ifelse(ybegin[1] < yend[1], min(ybegin), max(ybegin)),
      "zbegin" = ifelse(zbegin[1] < zend[1], min(zbegin), max(zbegin)),
      "xend" = ifelse(xbegin[1] < xend[1], max(xend), min(xend)),
      "yend" = ifelse(ybegin[1] < yend[1], max(yend), min(yend)),
      "zend" = ifelse(zbegin[1] < zend[1], max(zend), min(zend)),
    )
}
