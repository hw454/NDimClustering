#' Test for `find_low_percent_col`
#'
#' @description Test the function which is supposed to find the
#' row names for the points in the lowest `col_percent` of each column.
#'
#' * Check the output is list
#' * Check the number of terms is 10 percent of less.
#' * Check the min value is in the list
#' * check all outputs are in the original row names
#'
#' @family tests
#'
#' @export
test_find_low_percent_col <- function() {
    r <- 20
    c <- 3
    p <- 0.1
    col <- 3
    na_rm <- TRUE
    data_mat <- matrix(runif(r * c, 0, 20), nrow = r)
    rownames(data_mat) <- 1:r
    colnames(data_mat) <- 1:c
    low_rows <- find_low_percent_col(data_mat,
                                col = col,
                                col_percent = p,
                                na_rm = na_rm)
    testit::assert("...Should find character", is.character(low_rows))
    testit::assert("...More than 10% of terms found",
        length(low_rows) <= 0.1 * r)
    testit::assert("...All rows should be rownames in the data_mat",
        low_rows %in% rownames(data_mat))
    min_row <- rownames(data_mat)[which.min(data_mat[, col])]
    testit::assert("...Min should be in lowest percent",
        min_row %in% low_rows)
}