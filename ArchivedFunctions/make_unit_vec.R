#' Make a unit vector
#'
#' @param i the position for 1 to be located.
#' @param tot_n the total number of terms to be in the vector.
#'
#' @description
#' Create vector of zeros "u",
#' allocate "i"th position of "u" to 1.
#'
#' @return u
#'
#' @export
make_unit_vec <- function(i, tot_n) {
  u <- integer(tot_n)
  u[i] <- 1
  return(u)
}