#' Test the function which finds the cluster centroids.
#'
#' @description
#' Tests:
#' * Each cluster number should have a corresponding centroid
#' * The dimensions of the centroid should be the same as the data matrix.
#'
#' @family tests
#' @export
test_calc_clust_cent <- function(){
  nclust <- 3
  snp_id <- "rs35662"
  dummy_traits <- c("T1", "T2", "T3")
  dummy_snps <- c(snp_id, "rs301884", "rs69696",
                  "rs4096", "rs646464", "rs1234",
                  "rs203052", "rs646416", "rs98765",
                  "rs124532", "rs161616", "rs001122",
                  "rs111144", "rs222224", "rs01234",
                  "rs409664", "rs323223", "rs998877",
                  "rs309664", "rs453223", "rs668877",
                  "rs209664", "rs673223", "rs558877",
                  "rs109664", "rs893223", "rs448877")
  expec_cols <- c("clust_dist")
  num_axis <- length(dummy_traits)
  nsnps <- length(dummy_snps)
  # Create the test data matrix
  b_mat <- matrix(runif(nsnps * num_axis, 0, 10), nrow = nsnps)
  colnames(b_mat) <- dummy_traits
  rownames(b_mat) <- dummy_snps
  # Create the test p-val matrix
  p_mat <- matrix(runif(nsnps * num_axis, 0, 1), nrow = nsnps)
  colnames(p_mat) <- dummy_traits
  rownames(p_mat) <- dummy_snps

  # Assign points to clusters
  c_num_rnd <- sample.int(nclust,
    size = nsnps,
    replace = TRUE
  )
  cluster_df <- data.frame(
    row.names = dummy_snps,
    clust_num = c_num_rnd
  )
  print("... Test without p-values")
  bin_p_clust <- 0
  cent_df <- calc_clust_cent(clustnum_df = cluster_df,
    b_mat = b_mat,
    p_mat = p_mat,
    bin_p_clust = bin_p_clust
  )
  expec_rows <- length(unique(cluster_df$clust_num))
  expec_cols <- ncol(b_mat)
  testit::assert("... Output is not a dataframe",
                 is.data.frame(cent_df))
  testit::assert("... Wrong number of rows",
                 nrow(cent_df) == expec_rows)
  testit::assert("... Wrong number of columns",
                 ncol(cent_df) == expec_cols)
  testit::assert("... Wrong column label",
                 colnames(cent_df) == colnames(b_mat))
  print("... Test with p-values")
  bin_p_clust <- 1
  cent_df <- calc_clust_cent(clustnum_df = cluster_df,
    b_mat = b_mat,
    p_mat = p_mat,
    bin_p_clust = bin_p_clust
  )
  testit::assert("... Output is not a dataframe",
                 is.data.frame(cent_df))
  testit::assert("... Wrong number of rows",
                 nrow(cent_df) == expec_rows)
  testit::assert("... Wrong number of columns",
                 ncol(cent_df) == expec_cols)
  testit::assert("... Wrong column label",
                 colnames(cent_df) == colnames(b_mat))
}