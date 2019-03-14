#' Append columns to matrix x corresponding to columns of matrix y
#'
#' @param x a matrix (missing rows corresponding to matrix y)
#' @param colnames_y vector of colnames of matrix y 
#' @return matrix x with missing cols appended (entries are NA)
#' @examples
#' 
#' @export add_missing_cols
#' 
add_missing_cols <- function(x, colnames_y) {
  colnames_y <- colnames_y[! colnames_y %in% colnames(x)]
  x <- cbind(x, matrix(NA, nrow=nrow(x), ncol=length(colnames_y),
                       dimnames=list(rownames(x), colnames_y)))
}