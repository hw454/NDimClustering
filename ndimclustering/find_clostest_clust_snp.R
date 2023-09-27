#' Find the cluster whose centroid is closest to snp_id
#'
#' @description
#'   Find the distance between the snp and each cluster centroid.
#'   If the distance is less then the currently assigned cluster then
#'   reassign the cluster number and cluster distance.
#'
#' @param snp_id the label for the snp
#' @param b_mat the matrix of the data. The rows correpond to the snps.
#' @param cluster_df dataframe of cluster membership. Each row is a snp,
#'   it has columns for\:
#'    * "clust_num" the number of the allocated cluster.
#'    * "clust_dist" the distane to the centre of the cluster.
#'    * "clust_prob" probability the snp is in the assigned cluster.
#' @param centroids_df dataframe of the centroid co-ordinates. Each row
#'   corresponds to a cluster number and the columns are the trait axes.
#' @param norm_typ The type of norm to be used to calculate distance. The
#'    default is the Froebenius norm "F".
#' @param max_dist The maximum distance between any two points. Used
#'   for initialisation if snp is not already in the cluster assignment.
#'
#' @details
#'   The output dataframe "snp_clust_df"  has a row labelled by the "snp_id".
#'   The columns are\:
#'   * "clust_num" The number of the closest cluster
#'   * "clust_dist" The distance from the snp to the centroid of that cluster.
#'   * "clust_prob" probability the snp is in the cluster it's been assigned.
#'
#' @return snp_clust_df
#'
#' @export
find_closest_clust_snp <- function(snp_id, b_mat, cluster_df,
                              centroids_df,
                              norm_typ = "F",
                              max_dist = 10.0) {
  snp_score <- b_mat[snp_id, ]
  if (!(snp_id %in% row.names(cluster_df))){
    snp_dist <- max_dist
    snp_clust_num <- 0
  } else {
    snp_dist <- cluster_df[snp_id, "clust_dist"]
    snp_clust_num <- cluster_df[snp_id, "clust_num"]
  }
  for (c_num in rownames(centroids_df)) {
    dist <- norm(data.matrix(
              na.omit(centroids_df[c_num, ] - snp_score)), norm_typ)
    if (dist < snp_dist) {
      snp_dist <- dist
      snp_clust_num <- c_num
    }
   }
  snp_cluster_df <- data.frame(
    row.names = snp_id,
    clust_num = snp_clust_num,
    clust_dist = snp_dist,
    clust_prob = 1.0
   )
  return(snp_cluster_df)
}