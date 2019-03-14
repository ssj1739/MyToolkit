#' Collapse rows of matrix according to grouping variable
#'
#' @importFrom magrittr "%>%"
#' @param mat input matrix
#' @param group_df mapping from sample_var (row names of matrix) to specified group_var
#' @param collapse_fun function to collapse samples into groups
#' @param group_var grouping variable
#' @param sample_var sample variable; values must match row names of matrix
#' @param do_parallel T or F; parallelize collapse operation
#' @param ... additional arguments to `collapse_fun`
#' @return matrix with collapsed values
#' @examples
#' nrow <- 100
#' ncol <- 20
#' x <- rnorm(nrow*ncol) %>% matrix(nrow=nrow, ncol=ncol) %>%
#'   magrittr::set_rownames(str_c("Guide_", 1:nrow)) %>%
#'   magrittr::set_colnames(str_c("CL", 1:ncol))
#' df <- data_frame(GuideID=rownames(x),
#'                  GeneID=sample(str_c("Gene_", LETTERS), nrow, replace=T))
#' res <- cdsr::collapse_rows_of_matrix(x, df, group_var=GeneID, sample_var=GuideID,
#'                                collapse_fun=Matrix::colMeans, na.rm=T)
#' @export collapse_rows_of_matrix

collapse_rows_of_matrix <- function(mat, group_df,
                                    collapse_fun=Matrix::colMeans,
                                    group_var=Group,
                                    sample_var=Sample,
                                    do_parallel=F, ...) {
  
  GroupVar <- rlang::enquo(group_var)
  SampleVar <- rlang::enquo(sample_var)
  
  GroupStr <- deparse(substitute(group_var))
  SampleStr <- deparse(substitute(sample_var))
  
  group_df <- group_df %>%
    dplyr::filter(!is.na(!!GroupVar), !is.na(!!SampleVar)) %>%
    dplyr::filter((!!SampleVar) %in% rownames(mat)) %>%
    dplyr::distinct(GroupVar, SampleVar)
  
  single_group_df <- group_df %>% dplyr::group_by(!!GroupVar) %>% dplyr::filter(n() == 1)
  group_df <- group_df %>% dplyr::group_by(!!GroupVar) %>% dplyr::filter(n() > 1)
  
  mat_1 <- mat[single_group_df[[SampleStr]],, drop=F] %>%
    magrittr::set_rownames(single_group_df[[GroupStr]])
  
  mat_n <-
    plyr::daply(group_df, GroupVar, function(g) {
      m <- mat[unique(g[[SampleStr]]),, drop=F]
      do.call(collapse_fun, list(m, ...))
    }, .parallel=do_parallel)
  
  rbind(mat_1, mat_n) %>% {.[sort(rownames(.)), ]}
  
}