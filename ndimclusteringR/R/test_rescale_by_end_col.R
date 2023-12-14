#' Test the `rescale_by_end_col` function
#'
#' @export
#'
#' @family tests
test_rescale_by_end_col <- function() {
  # Test whether the function for checking the NaNs in a column works.
  dummy_traits <- c("T1", "T2", "T3", "T4", "T5", "T6")
  dummy_snps <- c("rs35662", "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  b_mat <- matrix(stats::runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  se_mat <- matrix(stats::runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  colnames(se_mat) <- dummy_traits
  rownames(se_mat) <- dummy_snps
  mat_list <- list("beta" = b_mat, "se" = se_mat)
  test1 <- rescale_by_end_col(mat_list$beta)
  expec_ncols <- ncol(b_mat) - 1
  expec_nrows <- nsnps
  testit::assert("Wrong number of columns after crop", {
    ncol(test1) == expec_ncols
  })
  testit::assert("Wrong number of rows after crop", {
    nrow(test1) == expec_nrows
  })
}