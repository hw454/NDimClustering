devtools::install(("ndimclusteringR"))
library("ndimclusteringR")
d <- 60
M2 <- matrix(runif(d * d, 0, 10), nrow = d)
M3 <- matrix(runif(d * d, 0, 1), nrow = d)
M4 <- matrix(runif(d * d), nrow = d)
colnames(M2) <- seq_len(ncol(M2))
colnames(M3) <- seq_len(ncol(M3))
colnames(M4) <- seq_len(ncol(M4))
rownames(M2) <- seq_len(nrow(M2))
rownames(M3) <- seq_len(nrow(M3))
rownames(M4) <- seq_len(nrow(M4))
narm <- TRUE
xbar <- apply(M2, 2, mean, na.rm = narm)
se <- apply(M2, 2, stats::sd, na.rm = narm)
# Compute the Sample SE with na_rm
out_list <- lapply(colnames(M2), calc_col_scale,
data = M2, mu = xbar, se = se)
out_mat <- Reduce(cbind, out_list)
M2h <- calc_scale_mat(M2)
num_axis <- ncol(M2)
stats::cor(M2h, method = "pearson", use = "pairwise.complete.obs")
out_list <- find_principal_components(M2, M3, M4)
cluster_out <- cluster_kmeans_min(out_list, 3, space_typ = "angle")
cluster_df <- cluster_out$clusters
cluster_df["num_axis"] <- num_axis
cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test",
                         "ndim_typ" = "test")
c1 <- colnames(M2)[1]
c2 <- colnames(M2)[2]
plot_clust_scatter(cluster_df, out_list$beta, out_list$se, iter_traits,
                   num_axis = num_axis)
plot_clust_exp_out_scatter(cluster_df, M2, M3, iter_traits,
                   exp_pheno = c1,
                   out_pheno = c2,
                   num_axis = num_axis)
print("Repeat test with correlated data. PCs should still be uncorrelated")
d <- 10
a <- 5
b <- 2
X <- matrix(rep(1:d, times = d), nrow = d)
M2 <- matrix(runif(d * d, 0, 1), nrow = d)
M3 <- matrix(runif(d * d, 0, 1), nrow = d)
M4 <- matrix(runif(d * d), nrow = d)
print(X)
print(M2)
M1 <- (a + M2) * X + b
M1[, 1] <- 1:d
colnames(M1) <- seq_len(ncol(M1))
colnames(M3) <- seq_len(ncol(M3))
colnames(M4) <- seq_len(ncol(M4))
rownames(M1) <- seq_len(nrow(M1))
rownames(M2) <- seq_len(nrow(M2))
rownames(M3) <- seq_len(nrow(M3))
rownames(M4) <- seq_len(nrow(M4))
narm <- TRUE
xbar <- apply(M1, 2, mean, na.rm = narm)
se <- apply(M1, 2, stats::sd, na.rm = narm)
# Compute the Sample SE with na_rm
out_list <- lapply(colnames(M1), calc_col_scale,
                    data = M1,
                    mu = xbar,
                    se = se)
out_mat <- Reduce(cbind, out_list)
M1h <- calc_scale_mat(M1)
num_axis <- ncol(M1)
stats::cor(M1h, method = "pearson", use = "pairwise.complete.obs")
out_list <- find_principal_components(M1, M3, M4)
cluster_out <- cluster_kmeans_min(out_list, 3, space_typ = "angle")
cluster_df <- cluster_out$clusters
cluster_df["num_axis"] <- num_axis
cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test2",
                         "ndim_typ" = "test2")
plot_clust_scatter(cluster_df, out_list$beta, out_list$se, iter_traits,
                   num_axis = num_axis)
plot_clust_exp_out_scatter(cluster_df, M1, M3, iter_traits,
                   exp_pheno = c1,
                   out_pheno = c2,
                   num_axis = num_axis)
print("Repeat test with two sets of correlated data. PCs should still be uncorrelated")
d <- 50
a1 <- 10
b1 <- 0
a2 <- 0.5
b2 <- 0
X <- matrix(rep(1:d, times = d), nrow = d)
M1 <- matrix(runif(d * d, 0, 1), nrow = d)
M2 <- matrix(runif(d * d, 0, 1), nrow = d)
M3 <- matrix(runif(2 * d * d, 0, 1), nrow = 2 * d)
M4 <- matrix(runif(2 * d * d), nrow = 2 * d)
M1a <- (a1 + M1) *  X + b1
M1b <- (a2 + M2) *  X + b2
M1a[, 1] <- 1:d
M1b[, 1] <- 1:d
M1 <- rbind(M1a, M1b)
colnames(M1) <- seq_len(ncol(M1))
colnames(M3) <- seq_len(ncol(M3))
colnames(M4) <- seq_len(ncol(M4))
rownames(M1) <- seq_len(nrow(M1))
rownames(M3) <- seq_len(nrow(M3))
rownames(M4) <- seq_len(nrow(M4))
narm <- TRUE
xbar <- apply(M1, 2, mean, na.rm = narm)
se <- apply(M1, 2, stats::sd, na.rm = narm)
# Compute the Sample SE with na_rm
out_list <- lapply(colnames(M1), calc_col_scale,
                    data = M1,
                    mu = xbar,
                    se = se)
out_mat <- Reduce(cbind, out_list)
M1h <- calc_scale_mat(M1)
num_axis <- ncol(M1)
stats::cor(M1h, method = "pearson", use = "pairwise.complete.obs")
out_list <- find_principal_components(M1, M3, M4)
cluster_out <- cluster_kmeans_min(out_list, 3, space_typ = "angle")
cluster_df <- cluster_out$clusters
cluster_df["num_axis"] <- num_axis
cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test3",
                         "ndim_typ" = "test3")
plot_clust_scatter(cluster_df, out_list$beta, out_list$se, iter_traits,
                   num_axis = num_axis)
plot_clust_exp_out_scatter(cluster_df, M1, M3, iter_traits,
                   exp_pheno = c1,
                   out_pheno = c2,
                   num_axis = num_axis)
# Results shuld show scatter plot with points not in a straight line.

print("Repeat with prcomp()")
d <- 52
a1 <- 10
b1 <- 0
a2 <- 0.5
b2 <- 0
X <- matrix(rep(1:d, times = d), nrow = d)
M1 <- matrix(runif(d * d, 0, 1), nrow = d)
M2 <- matrix(runif(d * d, 0, 1), nrow = d)
M3 <- matrix(runif(2 * d * d, 0, 1), nrow = 2 * d)
M4 <- matrix(runif(2 * d * d), nrow = 2 * d)
M1a <- (a1 + M1) *  X + b1
M1b <- (a2 + M2) *  X + b2
M1a[, 1] <- 1:d
M1b[, 1] <- 1:d
M1 <- rbind(M1a, M1b)
colnames(M1) <- seq_len(ncol(M1))
colnames(M3) <- seq_len(ncol(M3))
colnames(M4) <- seq_len(ncol(M4))
rownames(M1) <- seq_len(nrow(M1))
rownames(M3) <- seq_len(nrow(M3))
rownames(M4) <- seq_len(nrow(M4))
narm <- TRUE
xbar <- apply(M1, 2, mean, na.rm = narm)
se <- apply(M1, 2, stats::sd, na.rm = narm)
# Compute the Sample SE with na_rm
out_list <- lapply(colnames(M1), calc_col_scale,
                    data = M1,
                    mu = xbar,
                    se = se)
out_mat <- Reduce(cbind, out_list)
M1h <- calc_scale_mat(M1)
num_axis <- ncol(M1)
stats::cor(M1h, method = "pearson", use = "pairwise.complete.obs")
pca_beta <- stats::prcomp(M1,
      center = TRUE,
      scale = TRUE,
      rank = 3)
t_mat <- pca_beta$rotation
b_pc_mat <- pca_beta$x
p_pc_mat <- transform_coords_in_mat(M3, t_mat)
se_pc_mat <- transform_coords_in_mat(M4, t_mat)
pca_list <- list("beta" = b_pc_mat,
                   "pval" = p_pc_mat,
                   "se" = se_pc_mat,
                   "transform" = t_mat)
cluster_pca <- cluster_kmeans_min(pca_list, 3, space_typ = "angle")
cluster_df <- cluster_pca$clusters
cluster_df["num_axis"] <- num_axis
cluster_df <- tibble::rownames_to_column(cluster_df, "snp_id")
iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test4",
                         "ndim_typ" = "test4")
plot_clust_scatter(cluster_df, pca_list$beta, pca_list$se, iter_traits,
                   num_axis = num_axis)
plot_clust_exp_out_scatter(cluster_df, M1, M3, iter_traits,
                   exp_pheno = c1,
                   out_pheno = c2,
                   num_axis = num_axis)