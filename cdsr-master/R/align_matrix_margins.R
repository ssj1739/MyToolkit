#' Align matrix margins
#'
#' @importFrom magrittr "%<>%"
#' @param mat_list list of matrices (columns are genes, rows are samples)
#' @param na_rows if TRUE, add NAs for missing rows, otherwise they are removed
#' @param na_cols if TRUE, add NAs for missing columns, otherwise they are removed
#' @param use_dims either an integer or character vector of length 1
#' @return new list of matrices
#' @examples
#' 
#' @export align_matrix_margins
#' 
align_matrix_margins <- function(mat_list, na_rows=TRUE, na_cols=FALSE,
                                  use_dims=NULL) {
  
  data_rownames <- plyr::llply(mat_list, rownames)
  data_colnames <- plyr::llply(mat_list, colnames)
  
  if (is.null(use_dims)) {
    if (na_rows) {
      row_collapse <- union
    } else {
      row_collapse <- intersect
    }
    
    if (na_cols) {
      col_collapse <- union
    } else {
      col_collapse <- intersect
    }
    
    data_rownames <- Reduce(row_collapse, data_rownames)
    data_colnames <- Reduce(col_collapse, data_colnames)
    
    if (length(data_rownames) < 1) stop("Error: no common rownames")
    if (length(data_colnames) < 1) stop("Error: no common colnames")
    
    if (na_rows) {
      mat_list <- plyr::llply(mat_list, cdsr::add_missing_rows, data_rownames)
    }
    if (na_cols) {
      mat_list <- plyr::llply(mat_list, cdsr::add_missing_cols, data_colnames)
    }
  } else {
    stopifnot(!is.null(data_rownames[[use_dims]]))
    stopifnot(!is.null(data_colnames[[use_dims]]))
    
    data_rownames <- data_rownames[[use_dims]]
    data_colnames <- data_colnames[[use_dims]]
    
    mat_list <- plyr::llply(mat_list, cdsr::add_missing_rows, data_rownames)
    mat_list <- plyr::llply(mat_list, cdsr::add_missing_cols, data_colnames)
    
  }
  
  mat_list <- plyr::llply(mat_list,
                    function(d) d[data_rownames, data_colnames])
  
  if(is.null(names(mat_list))){
    mat_list %<>% magrittr::set_names(paste0("value_", 1:length(.)))
  }
  
  return(mat_list)
  
}