#' Calculate correlations between matrices, with significance measures
#'
#' @param x matrix of values
#' @param y second matrix of values
#' @param adjust p-value correction method (fed into p.adjust).
#'
#' @return: list of matrices with elements r: correlations, n: number observations, t: t-stats, p: p-values, se: standard errors, q: q-values
#' @export
#'
#' @examples
corr_test <- function(x, y, adjust) {
  #take from Nick Clarks' post: https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
  n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
  r <- cor(x, y, use = "pairwise.complete.obs")[,1] # MUCH MUCH faster than corr.test()
  cor2pvalue = function(r, n) {
    t <- (r*sqrt(n-2))/sqrt(1-r^2)
    p <- 2*(1 - pt(abs(t),(n-2)))
    se <- sqrt((1-r*r)/(n-2))
    out <- list(r, n, t, p, se)
    names(out) <- c("r", "n", "t", "p", "se")
    return(out)
  }
  # get a list with matrices of correlation, pvalues, standard error, etc.
  result = cor2pvalue(r,n)
  if (adjust != 'none') {
    result$q <- p.adjust(result$p, method = adjust)
  }
  return(result)
}
