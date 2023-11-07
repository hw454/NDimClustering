devtools::install(("ndimclusteringR"))
library("ndimclusteringR")

number_row_col_names <- function(mat) {
    colnames(mat) <- seq_len(ncol(mat))
    rownames(mat) <- seq_len(nrow(mat))
    return(mat)
}
form_beta_corr <- function(n_path, d, rand_shift, a_list, b_list) {
    a <- a_list[n_path]
    b <- b_list[n_path]
    x_mat <- matrix(rep(1:d, times = d), nrow = d)
    beta <- (a * (2 ** n_path) + rand_shift) * x_mat + b
    beta[, 1] <- 1:d
    return(beta)
}
create_pc <- function(i) {
    return(paste0("PC", i))
}
test_pc_function <- function(d = 10, num_path = 0, pc_type = "HW", n_pc = 3) {
    print(paste("Test", pc_type, "with", num_path, "pathways"))
    rand_mat <- matrix(runif(d * d, 0, 1), nrow = d)
    if (num_path > 0) {
        a_list <- runif(num_path + 1, 0, 10)
        b_list <- runif(num_path + 1, 0, 10)
        mat_list <- lapply(1:(num_path + 1),
                           form_beta_corr,
                           d = d,
                           rand_shift = rand_mat,
                           a_list = a_list,
                           b_list = b_list)
        dummy_beta <- Reduce(rbind, mat_list)
    } else {
        dummy_beta <- rand_mat
    }
    dummy_se <- matrix(runif((num_path + 1) * d * d, 0, 1),
                             nrow = (num_path + 1) * d)
    dummy_p <- matrix(runif((num_path + 1) * d * d),
                            nrow = (num_path + 1) * d)
    dummy_beta <- number_row_col_names(dummy_beta)
    dummy_se <- number_row_col_names(dummy_se)
    dummy_p <- number_row_col_names(dummy_p)
    # Compute the Sample SE with na_rm
    num_axis <- ncol(dummy_beta)
    if (pc_type == "HW") {
        out_list <- find_principal_components(dummy_beta, dummy_se, dummy_p, np = n_pc)
        out_list$se <- calc_scale_mat(out_list$se)
        pc_cols <- lapply(1:n_pc, create_pc)
        colnames(out_list$se) <- pc_cols
    } else if (pc_type == "prcomp") {
    pca_beta <- stats::prcomp(dummy_beta,
      center = TRUE,
      scale = TRUE,
      rank = n_pc)
    t_mat <- pca_beta$rotation
    b_pc_mat <- pca_beta$x
    p_pc_mat <- transform_coords_in_mat(dummy_p, t_mat)
    se_pc_mat <- transform_coords_in_mat(dummy_se, t_mat)
    # Rescale se so it is within (0,1)
    se_pc_mat <- calc_scale_mat(se_pc_mat)
    pc_cols <- lapply(1:n_pc, create_pc)
    colnames(se_pc_mat) <- pc_cols
    out_list <- list("beta" = b_pc_mat,
                   "pval" = p_pc_mat,
                   "se" = se_pc_mat,
                   "transform" = t_mat)
    }
    print("Cluster")
    print(out_list)
    cluster_out <- cluster_kmeans_min(out_list, nclust = 3, space_typ = "angle")
    cluster_df <- cluster_out$clusters
    cluster_df["num_axis"] <- num_axis
    cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
    iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test",
                         "ndim_typ" = "test",
                         "res_dir" = paste0("PC_TestResults/", pc_type, "_", num_path + 1, "paths")
    )
    c1 <- colnames(dummy_beta)[1]
    c2 <- colnames(dummy_beta)[2]
    print("Plot")
    plot_clust_scatter(cluster_df, out_list$beta, out_list$se, iter_traits,
                   num_axis = num_axis)
    plot_clust_exp_out_scatter(cluster_df, dummy_beta, dummy_se, iter_traits,
                   exp_pheno = c1,
                   out_pheno = c2,
                   num_axis = num_axis)
}
pc_type_1 <- "HW"
pc_type_2 <- "prcomp"
num_paths_list <- 0:3
pc_list <- list(pc_type_1, pc_type_2)
d <- 30
for (p_type in pc_list){
    for (np in num_paths_list){
        test_pc_function(d,
                        num_path = np,
                        pc_type = p_type)
    }
}
print("Complete")