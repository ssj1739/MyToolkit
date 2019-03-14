#' Measure normality likelihood ratio statistic, as described in MacDonald et al Cell 2017
#'
#' @param vec vector of values, at least 10 data points required (much more probably better)
#'
#' @return twice log likelihood ratio of skewed t-dist fit compared to guassian fit
#' @export
#'
#' @examples
normLRT_test <- function(vec) {
  library(MASS); library(sn)
  min_length <- 10
  stopifnot(is.vector(vec))
  vec <- vec[!is.na(vec)]
  if (length(vec) < min_length) {return(NA)}
  st_LL <- tryCatch({
    data.frame(data = vec) %>%
      selm(data ~ 1, family = "ST", data = .) %>%
      logLik %>%
      as.numeric()
  }, error = function(e) {
      starting_nu <- 3
      print(sprintf('Initial fit failed, trying initial nu of %f', starting_nu))
      init_mod <- data.frame(data = vec) %>%
        selm(data ~ 1, family = "ST", fixed.param = list(nu=starting_nu), data = .)
      st_LL <- tryCatch({
        data.frame(data = vec) %>%
          selm(data ~ 1, family = "ST", data = ., start = c(coef(init_mod, param.type = 'DP'), list(nu = starting_nu))) %>%
          logLik %>%
          as.numeric()
      }, error = function(e) {
        print(sprintf('Second fit failed, fixing nu at %f', starting_nu))
        init_mod %>%
          logLik %>%
          as.numeric()
      })
      return(st_LL)
   })
  n_LL <- fitdistr(vec, 'normal')$loglik
  return(2*(st_LL - n_LL))
}
