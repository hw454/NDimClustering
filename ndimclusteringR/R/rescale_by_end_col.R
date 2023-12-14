#' Multiply all by the final terms in a row by the final term.
#' Repeat for all rows.
#'
#' @export
#' @family utility
rescale_by_end_col <- function(mat){
  rescale_row <- function(r, mat){
    nc <- ncol(mat)
    nfin <- nc - 1
    row_new <- mat[r, 1:nfin] * mat[r, nc]
    return(row_new)
  }
  nr <- nrow(mat)
  new_row_list <- lapply(1:nr, rescale_row,
                         mat = mat)
  mat <- Reduce(rbind, new_row_list)
  return(mat)
}