#' Test the `crop_matrices` function
#'
#' @export
#'
#' @family tests
test_crop_matrices <- function() {
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
  trait_df <- list("phenotype" = dummy_traits)
  out_pheno <- "T1"
  exp_pheno <- "T2"
  n_rows <- 2
  n_col0 <- 1
  n_col1 <- 4
  print("... Testing the exp and outcome within range")
  test1 <- crop_matrices(mat_list = mat_list,
    trait_df = trait_df,
    out_pheno = out_pheno,
    exp_pheno = exp_pheno,
    n_rows = n_rows,
    n_col0 = n_col0,
    n_col1 = n_col1
  )
  expec_ncols <- 1 + n_col1 - n_col0
  expec_nrows <- n_rows
  testit::assert("Wrong number of columns after crop", {
    ncol(test1$beta) == expec_ncols
  })
  testit::assert("Wrong number of rows after crop", {
    nrow(test1$beta) == expec_nrows
  })
  testit::assert("Exposure and outcome not in columns after crop", {
    c(exp_pheno, out_pheno) %in% colnames(test1$beta)
  })
  print("... Testing the exp within range and outcome outside range")
  out_pheno <- "T6"
  test1 <- crop_matrices(mat_list = mat_list,
    trait_df = trait_df,
    out_pheno = out_pheno,
    exp_pheno = exp_pheno,
    n_rows = n_rows,
    n_col0 = n_col0,
    n_col1 = n_col1
  )
  expec_ncols <- 2 + n_col1 - n_col0
  expec_nrows <- n_rows
  testit::assert("Wrong number of columns after crop", {
    ncol(test1$beta) == expec_ncols
  })
  testit::assert("Wrong number of rows after crop", {
    nrow(test1$beta) == expec_nrows
  })
  testit::assert("Exposure and outcome not in columns after crop", {
    c(exp_pheno, out_pheno) %in% colnames(test1$beta)
  })
  print("... Testing the exp within range and outcome outside range")
  exp_pheno <- "T5"
  out_pheno <- "T6"
  test1 <- crop_matrices(mat_list = mat_list,
    trait_df = trait_df,
    out_pheno = out_pheno,
    exp_pheno = exp_pheno,
    n_rows = n_rows,
    n_col0 = n_col0,
    n_col1 = n_col1
  )
  expec_ncols <- 3 + n_col1 - n_col0
  expec_nrows <- n_rows
  testit::assert("Wrong number of columns after crop", {
    ncol(test1$beta) == expec_ncols
  })
  testit::assert("Wrong number of rows after crop", {
    nrow(test1$beta) == expec_nrows
  })
  testit::assert("Exposure and outcome not in columns after crop", {
    c(exp_pheno, out_pheno) %in% colnames(test1$beta)
  })
  print("... Testing the nrows outside range")
  n_rows <- length(dummy_snps) + 2
  test1 <- crop_matrices(mat_list = mat_list,
    trait_df = trait_df,
    out_pheno = out_pheno,
    exp_pheno = exp_pheno,
    n_rows = n_rows,
    n_col0 = n_col0,
    n_col1 = n_col1
  )
  expec_ncols <- 3 + n_col1 - n_col0
  expec_nrows <- min(n_rows, length(dummy_snps))
  testit::assert("Wrong number of columns after crop", {
    ncol(test1$beta) == expec_ncols
  })
  testit::assert("Wrong number of rows after crop", {
    nrow(test1$beta) == expec_nrows
  })
  testit::assert("Exposure and outcome not in columns after crop", {
    c(exp_pheno, out_pheno) %in% colnames(test1$beta)
  })
}