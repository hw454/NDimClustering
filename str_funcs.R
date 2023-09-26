# Directory function
set_directory <- function(res_dir0, iter_traits) {
  #' Set the directory for the results using the base directory
  #' and the iteration parameters
  res_dir <- paste0(res_dir0, method_str(iter_traits), "/")
  # Create results directory if it doesn't exist
  if (!dir.exists(res_dir)) {
    dir.create(res_dir)
  }
  return(res_dir)
}
#- String functions
method_str <- function(iter_traits) {
  #' Create the string that describes the type of method for filenames
  if (iter_traits$bp_on) {
    bp_str <- "_bpON"
  } else {
    bp_str <- "_bpOFF"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "_clustprobON"
  } else {
    clust_prob_str <- "_clustprobOFF"
  }
  return(paste0(iter_traits$clust_typ, bp_str, clust_prob_str))
}

desc_str <- function(iter_traits) {
  #' Create the string that describes the type of method for titles and captions
  if (iter_traits$bp_on) {
    bp_str <- "bp on"
  } else {
    bp_str <- "bp off"
  }
  if (iter_traits$clust_prob_on) {
    clust_prob_str <- "ClustProb on"
  } else {
    clust_prob_str <- "ClustProb off"
  }
  return(paste(iter_traits$clust_typ, "and", bp_str, "and", clust_prob_str))
}