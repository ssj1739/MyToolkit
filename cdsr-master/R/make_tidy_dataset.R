#' Make a tidy dataset
#' 
#' @inheritParams align_matrix_margins
#' @param dim_names character vector of length two -
#' these become the column labels in the tidy dataset
#' for the rownames and colnames of the original matrices
#' @return a tidy dataframe
#' @description Takes a list of matrices with similar row and column names
#' and runs tidyr::gather across each to make a "long" dataset
#' with a column of values for each matrix in the list.
#' @examples 
#' CCLE_GE <-  taigr::load.from.taiga(
#'                 data.name="ccle-rnaseq-expression-genes",
#'                 data.version = 3) %>%
#'              set_colnames(stringr::str_match(colnames(.), '(.+) ')[,2]) %>%
#'              t()
#'              
#' MUT_MIS <- taigr::load.from.taiga(
#'                 data.name = 'ccle-mis-mut-binary-matrix',
#'                 data.version = 1,
#'                 transpose = T)
#' CN <- taigr::load.from.taiga(
#'                 data.name="ccle-copy-number-variants-hgnc-mapped",
#'                 data.version = 1)             
#'              
#' gene_feature_df <- cdsr::make_tidy_dataset(list(exp=CCLE_GE,
#'                                                 mut=MUT_MIS,
#'                                                 cn=CN),
#'                                            dim_names=c("Gene", "CellLine"))
#' head(gene_feature_df)
#'              
#' 
#' @export make_tidy_dataset
#' 
make_tidy_dataset <- function(mat_list, na_rows=TRUE, na_cols=FALSE,
                                use_dims=NULL,
                                dim_names=c("Gene", "Sample")) {
  
  datasets <- cdsr::align_matrix_margins(mat_list,
                                    na_rows=na_rows, na_cols=na_cols,
                                    use_dims=use_dims)
  
  tmp_df <- as.data.frame(datasets[[1]])

  rows <- rownames(tmp_df)
  cols <- colnames(tmp_df)
  
  tmp_df[, dim_names[1]] <- rows
  tmp_df <- tidyr::gather_(tmp_df, dim_names[2], "Temp", cols)
  tmp_df <- dplyr::select(tmp_df, -Temp)

  data_names <- names(datasets)
  names(data_names) <- data_names
  
  data_columns <- plyr::llply(data_names, function(d) {
    gath_df <- tidyr::gather_(as.data.frame(datasets[[d]]), dim_names[2], d, cols)
    return(gath_df[,d])
  })
  return(cbind(tmp_df, dplyr::as_data_frame(data_columns)))
  
}