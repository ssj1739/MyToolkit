#'Append rows to matrix x corresponding to rows of matrix y
#'
#' @param x a matrix (missing rows corresponding to matrix y)
#' @param rownames_y vector of rownames of matrix y 
#' @return matrix x with missing rows appended (entries are NA)
#' @examples
#' 
#' @export add_missing_rows
#' 
add_missing_rows <- function(x, rownames_y) {
  rownames_y <- rownames_y[! rownames_y %in% rownames(x)]
  x <- rbind(x, matrix(NA, nrow=length(rownames_y), ncol=ncol(x),
                       dimnames=list(rownames_y, colnames(x))))
}