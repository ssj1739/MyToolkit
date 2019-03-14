#' Convert gene symbols from one format to another
#'
#' @param genes A list of genes in either Ensembl, HGNC, or Entrez Gene format.
#' @param input_format 'ensembl', 'hgnc', or 'entrez'
#' @param output_format 'ensembl', 'hgnc', or 'entrez'
#' @param unique If TRUE, the mapping is one-to-one, returning the first occurrence. If FALSE, the mapping is one-to-many.
#'
#' @return A gene mapping from input to output format and a mapping of all repeats found
#'
#' @examples
#' my_gene_mappings <- map_genes(my_ensembl_genes, 'ensembl', 'hgnc', unique = F)
#' my_corresponding_hgnc_genes <- my_gene_mappings$mapping$hgnc
#' repeats_df <- my_gene_mappings$repeats
#' @export map_genes

map_genes <- function(genes, input_format, output_format, unique = TRUE) {
  require(taigr)
  require(dplyr)
  require(magrittr)
  map <- load.from.taiga(data.name='ensembl-to-hugo-gene-map-4550', data.version=3, data.file='gene_map')
  df <- data.frame(genes) %>%
    set_colnames(input_format) %>%
    left_join(map, by = input_format) %>%
    dplyr::select(input_format, output_format) %>%
    unique()

  repeats <- table(df[[input_format]])[table(df[[input_format]]) > 1] %>% names()
  repeats_df <- df[which(df[[input_format]] %in% repeats),]

  if (unique == T && nrow(repeats_df) > 0) {
    print('This function is returning a one-to-many mapping, and the repeats will be returned. The mapping with include only the first occurence.')
    # keep only first occurrence
    df <- df[match(unique(df[[input_format]]), df[[input_format]]),]
  }
  n <- sum(is.na(df[[output_format]]))
  if (n > 0) print(paste0('Note: ', n, ' genes were not mapped to ', output_format))
  return(list(mapping = df, repeats = repeats_df))
}
