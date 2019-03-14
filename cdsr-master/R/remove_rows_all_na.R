#' Remove rows from a matrix that are entirely NA
#'
#' @param x a matrix
#' @param indices logical, return indices to remove instead
#' @param verbose logical, display number of rows removed
#' @return the matrix with rows removed, unless indices is TRUE in which case the
#' indices of the rows with all NA entries are returned
#' @examples
#' mat <- matrix(c(NA, 2, NA, 4), nrow=2, ncol=2)
#' remove_rows_all_na(mat)
#' @export remove_rows_all_na
remove_rows_all_na <- function(x, indices=FALSE, verbose=TRUE) {
  i <- aaply(x, 1, function(r) {all(is.na(r))})
  if (indices) {
    return(which(i))
  } else {
    if(verbose){
      cat(paste0("\nRemoved ", sum(i), " row(s) from matrix...\n"))
    }
    return(x[ ! i, ])
  }
}