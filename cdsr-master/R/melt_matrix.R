
#' Convert a matrix into a data.frame
#'
#' @importFrom magrittr "%>%"
#' @param mat A matrix.
#' @param row_name String label for matrix row names
#' @param col_name String label for matrix column names
#' @param value_name String label for matrix entries
#' @return A data.frame containing columns for \code{row_name}, \code{col_name}, and \code{value_name}
#' @examples
#' MUT.DAM <- load.from.taiga(data.name = 'ccle-mut-data-binary-matrix', data.version = 1)
#' MUT.DAM.DF <- melt_matrix(MUT.DAM, row_name = "ccle_name", col_name = "gene", value_name="mutation")
#' @export melt_matrix

melt_matrix <- function(mat, row_name, col_name, value_name){
  
  require(reshape2)
  require(dplyr)
  require(magrittr)
  
  df <- mat %>%
          as.data.frame() %>%
          dplyr::mutate(row=row.names(.)) %>%
          reshape2::melt(id.vars="row", variable.name=col_name, value.name=value_name)

  colnames(df)[colnames(df) == "row"] <- row_name
  
  return(df)
  
}