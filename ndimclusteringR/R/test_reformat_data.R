#' Test the `reformat_data` function
#'
#' @export
#'
#' @family tests
test_reformat_data <- function() {
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
  data_matrices <- append(mat_list,
                          "trait_info" = trait_df)

  print("... Testing with angles off and pca off")
  fmt_mat_list <- reformat_data(data_matrices = data_matrices,
                                bin_angles = 0,
                                pca_type = "none")
  expec_nterms <- 8
  testit::assert("Wrong number of terms in list", {
    length(fmt_mat_list) == expec_nterms
  })
  expec_terms <- c("beta", "se", "pval", "trait_info",
    "beta_pc", "se_pc", "pval_pc", "transform"
  )
  testit::assert("Wrong terms in list", {
    all(expec_terms %in% names(fmt_mat_list))
  })

}