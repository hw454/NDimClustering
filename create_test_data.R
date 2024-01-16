number_row_col_names <- function(mat) {
  colnames(mat) <- seq_len(ncol(mat))
  rownames(mat) <- seq_len(nrow(mat))
  return(mat)
}
form_beta_corr <- function(n_path, d, rand_shift, a_list, b_list) {
  a <- a_list[n_path]
  b <- b_list[n_path]
  x_mat <- matrix(rep(1:d, times = d), nrow = d)
  beta <- (a * x_mat + b) + rand_shift
  beta[, 1] <- 1:d
  return(beta)
}
create_pc <- function(i) {
  return(paste0("PC", i))
}
create_test_data <- function(d = 50, np = 0) {
  print(paste("Create test data with", np, "pathways"))
  iter_dir <- paste0("paths", np, "/")
  rand_mat <- matrix(runif(d * d, -1, 1), nrow = d)
  if (np > 0) {
    a_list <- runif(np, 0, 3)
    b_list <- runif(np, 0, 0)
    mat_list <- lapply(1:np,
                       form_beta_corr,
                       d = d,
                       rand_shift = rand_mat,
                       a_list = a_list,
                       b_list = b_list)
    dummy_beta <- Reduce(rbind, mat_list)
    dummy_se <- matrix(runif(np * d * d),
                       nrow = np * d)
    dummy_p <- matrix(runif(np * d * d, 0, 1),
                      nrow = np * d)
  } else {
    dummy_beta <- rand_mat
    dummy_se <- rand_mat
    dummy_p <- rand_mat
  }
  dummy_beta <- number_row_col_names(dummy_beta)
  dummy_se <- number_row_col_names(dummy_se)
  dummy_p <- number_row_col_names(dummy_p)
  trait_list <- rep("Trait", ncol(dummy_beta))
  for (i in seq_len(ncol(dummy_beta))){
    trait_list[i] <- paste0(trait_list[i], i)
  }
  colnames(dummy_beta) <- trait_list
  colnames(dummy_se) <- trait_list
  colnames(dummy_p) <- trait_list
  cat_list <- rep("Exposure", ncol(dummy_beta))
  cat_list[1] <- "Outcome"
  trait_info <- list("pheno_category" = cat_list,
                     "phenotype" = trait_list)
  maindir <- "TestData/"
  betafilename <- "unstdBeta_df.csv"
  sefilename <- "unstdSE_df.csv"
  pvalfilename <- "pval_df.csv"
  traitfilename <- "trait_info_nfil.csv"
  if (!file.exists(maindir)) {
    dir.create(maindir)
  }
  if (!file.exists(paste0(maindir, iter_dir))) {
    dir.create(paste0(maindir, iter_dir))
  }
  write.csv(dummy_beta,
    paste0(maindir, iter_dir, betafilename),
    row.names = TRUE,
    col.names = TRUE,
    eol = "\n",
    sep = ","
  )
  write.csv(dummy_se,
    paste0(maindir, iter_dir, sefilename),
    row.names = TRUE,
    col.names = TRUE,
    eol = "\n",
    sep = ","
  )
  write.csv(dummy_p,
    paste0(maindir, iter_dir, pvalfilename),
    row.names = TRUE,
    col.names = TRUE,
    eol = "\n",
    sep = ","
  )
  write.csv(trait_info,
    paste0(maindir, iter_dir, traitfilename),
    row.names = FALSE,
    quote = FALSE,
    eol = "\n",
    sep = ","
  )
}

num_path_list <- 0:5
for (np in num_path_list){
  create_test_data(np = np)
}