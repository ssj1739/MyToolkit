% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_mapping.R
\name{map_genes}
\alias{map_genes}
\title{Convert gene symbols from one format to another}
\usage{
map_genes(genes, input_format, output_format, unique = TRUE)
}
\arguments{
\item{genes}{A list of genes in either Ensembl, HGNC, or Entrez Gene format.}

\item{input_format}{'ensembl', 'hgnc', or 'entrez'}

\item{output_format}{'ensembl', 'hgnc', or 'entrez'}

\item{unique}{If TRUE, the mapping is one-to-one, returning the first occurrence. If FALSE, the mapping is one-to-many.}
}
\value{
A gene mapping from input to output format and a mapping of all repeats found
}
\description{
Convert gene symbols from one format to another
}
\examples{
my_gene_mappings <- map_genes(my_ensembl_genes, 'ensembl', 'hgnc', unique = F)
my_corresponding_hgnc_genes <- my_gene_mappings$mapping$hgnc
repeats_df <- my_gene_mappings$repeats
}
