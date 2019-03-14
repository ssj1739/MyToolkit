#' Remove rows from a matrix that have any NA entries
#'
#' @inheritParams remove_rows_all_na
#' @return the matrix with rows removed, unless indices is TRUE
#' @examples
#' mat <- matrix(c(NA, 2, 3, 4), nrow=2, ncol=2)
#' remove_rows_any_na(mat)
#' @export remove_rows_any_na
remove_rows_any_na <- function(x, indices=FALSE, verbose=TRUE) {
  i <- aaply(x, 1, function(r) {any(is.na(r))})
  if (indices) {
    return(which(i))
  } else {
    if(verbose){
      cat(paste0("\nRemoved ", sum(i), " row(s) from matrix...\n"))
    }
    return(x[ ! i, ])
  }
}