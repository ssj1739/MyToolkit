
#' Conduct t-tests on data within a data.frame
#'
#' @importFrom magrittr "%>%"
#' @param df A data.frame
#' @param group_vars Vector of string names for columns by which to group
#' @param condition_var String name of column indicating population conditions
#' @param value_var String name of column containing value to be tested
#' @return A tibble containing a plot for each distinct combination of group_vars
#' @examples
#' plot_tbl <- t_test(mtcars, group_vars="vs", condition_var="am", value_var="mpg")
#' plot_tbl %>% dplyr::filter(vs == 1) %>% magrittr::extract2("plot")
#' @export t_test

t_test <- function(df, group_vars, condition_var, value_var){
  
  conditions <- unique(df[[condition_var]]) %>%
                  as.character()
  if(length(conditions) != 2){
    cat("\nt_test currently does not support pairwise tests among >2 conditions")
    stop()
  }
  
  df %>%
    dplyr::group_by_(group_vars) %>%
    dplyr::group_by_(condition_var, add=T) %>%
    dplyr::mutate(sample_id=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by_(group_vars) %>%
    tidyr::spread_(condition_var, value_var) %>%
    dplyr::do(test_result=t.test(.[[conditions[1]]], .[[conditions[2]]]),
              plot=ggplot2::ggplot(reshape2::melt(., measure.vars=conditions), 
                                   aes(x=variable, y=value, fill=variable)) +
                    geom_boxplot() +
                    stat_summary(fun.data = function(x){return(data.frame(y=mean(x), 
                                                                 label=paste0("n=",sum(!is.na(x)))))},
                                 geom="text", colour="white") +
                    labs(x=condition_var, y=value_var, title=paste(unique(magrittr::extract2(., group_vars))), sep="_") +
                    scale_fill_brewer(name=condition_var, palette = "Set1")) %>%
    dplyr::do(statistic=.[["test_result"]][["statistic"]], 
              p=.[["test_result"]][["p.value"]],
              diff_means=.[["test_result"]][["estimate"]][1] - .[["test_result"]][["estimate"]][2],
              plot=.[["plot"]]) %>%
    dplyr::do(plot=.[["plot"]] + annotate(geom="segment", 
                                          x=conditions[1], xend=conditions[2], 
                                          y=max(df[[value_var]]) + 0.05*max(df[[value_var]]), 
                                          yend=max(df[[value_var]]) + 0.05*max(df[[value_var]])) +
                                  annotate(geom="text",
                                           x=conditions[2],
                                           y=max(df[[value_var]]) + 0.06*max(df[[value_var]]),
                                           label=ifelse(.[["p"]] > 0.05, "",
                                                        "***")) +
                                  annotate(geom="text",
                                            x=conditions[1],
                                            y=max(df[[value_var]]) + 0.1*max(df[[value_var]]),
                                            label=paste0("p = ", format(round(.[["p"]], 2), scientific = T)))) %>%
    dplyr::ungroup() %>%
    magrittr::extract2("plot") %>%
    tibble::add_column(tibble::as_tibble(df[, group_vars] %>% as.data.frame() %>% magrittr::set_colnames(group_vars) %>% dplyr::distinct()),
                       plot=.)
  
}