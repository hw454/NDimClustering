#' Multiply all by the final terms in a row by the final term.
#' Repeat for all rows.
#'
#' @param matin matrix whose rows are to be rescaled by the final column
#' @export
#' @family preconditioning_functions
rescale_by_end_col <- function(matin) {
  rescale_row <- function(r, matin) {
    nc <- ncol(matin)
    nfin <- nc - 1
    row_new <- matin[r, 1:nfin] * matin[r, nc]
    return(row_new)
  }
  nr <- nrow(matin)
  new_row_list <- lapply(1:nr, rescale_row,
                         mat = matin)
  mat <- Reduce(rbind, new_row_list)
  rownames(mat) <- rownames(matin)
  return(mat)
}