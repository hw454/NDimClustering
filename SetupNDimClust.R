#' Data inputs needed: unstanderdised association matrix between IV and traits,
#' SE matrix of association matrix, t-stat matrix of association matrix,
#' P-value matrix of association matrix,
#' trait info (trait, description, sample size)

# variable set-up
threshmul <- 5.0
clust_threshold <- 1e-8
thresh_norm <- "F"
clust_norm <- "F"
nr <- 5 # Number of clusters
np <- 2 # Number of PCs
na_percent <- 0.25 # The percentage of a column that can acceptably be not NaN

# Testing dimensions
test <- 0 # testing switch
num_trait0 <- 380
num_trait1 <- 410
num_rows <- 50

# Create the data variables from the inputs
setup_algorithm_data(threshmul = threshmul,
                    clust_threshold = clust_threshold,
                    na_percent = na_percent,
                    nr = nr,
                    np = np,
                    clust_norm = clust_norm,
                    thresh_norm = thresh_norm)

setup_matrices(data_dir = data_dir,
              test = test,
              num_rows = num_rows,
              num_trait0 = num_trait0,
              num_trait1 = num_trait1)