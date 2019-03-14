#' ggduo module of scatter plot with linear fit 
#' @export lm_with_cor
lm_with_cor <- function(data, mapping, ..., method = "pearson") {
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  cor <- cor(x, y, method = method, use = 'p')
  ggplot(data, mapping, ...) +
    geom_point(size = .6) +
    stat_smooth(method = 'lm', color = 'blue', alpha = 0.1, size = 0.3) +
    annotate(geom = 'label', 
             x = median(x, na.rm = TRUE), y = 1.2,
             label = round(cor, 3), alpha = 0.3,
             hjust = 0, vjust = 0.5) +
    theme_bw()
}