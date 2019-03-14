#' Remove columns from a matrix that contain any NAs
#'
#' @inheritParams remove_rows_all_na
#' @return the matrix with columns removed, unless indices is TRUE
#' @examples
#' mat <- matrix(c(NA, 2, 3, 4), nrow=2, ncol=2)
#' remove_cols_any_na(mat)
#' @export remove_cols_any_na
remove_cols_any_na <- function(x, indices=FALSE, verbose=TRUE) {
  i <- aaply(x, 2, function(r) {any(is.na(r))})
  if (indices) {
    return(which(i))
  } else {
    if(verbose){
      cat(paste0("\nRemoved ", sum(i), " cols(s) from matrix...\n"))
    }
    return(x[, ! i])
  }
}