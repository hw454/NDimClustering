devtools::install(("ndimclustering"))
library("ndimclustering")
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
stats::cor(M2h, method = "pearson", use = "pairwise.complete.obs")
out_list <- find_principal_components(M2, M3, M4)
cluster_df <- cluster_kmeans_min(out_list, 3, space_typ = "angle")
iter_traits <- data.frame("bp_on" = TRUE,
                         "clust_prob_on" = TRUE,
                         "clust_typ" = "test",
                         "ndim_typ" = "test")
plot_clust_scatter(cluster_df$clusters, out_list$beta, out_list$se, iter_traits)

# Results shuld show scatter plot with points not in a straight line.