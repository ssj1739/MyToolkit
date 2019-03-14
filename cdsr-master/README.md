# cdsr: the R package for the Cancer Data Science team

This package is an R toolbox for the Cancer Data Science (CDS) team at The Broad Institute of MIT and Harvard. The scripts included here are meant to address common-use operations and synchronize the team's R vocabulary.

# Naming convention for functions

Functions follow the `underscore_separated` naming convention. All function names must consist of lowercase characters separated by the character `_`.

# Naming conventions for function arguments

Arguments to functions also follow the `underscore_separated` naming convention.

+ Matrices: `mat`
+ Data frames: `df`
+ Lists: `li`

# Contribute to the package

Functions must be added according to the following steps:

1. Create your `.R` file and move it to the `R` directory
2. Add `roxygen` style comments to your `.R` file ([details](http://r-pkgs.had.co.nz/man.html))
3. If your `.R` file contains a function called `your_new_function`, then run `?your_new_function` for a preview of the documentation
4. If you like how your documentation looks, navigate to the `cdsr` parent directory and run `devtools::document()`. This will generate a `man/your_new_function.Rd` containing your function documentation.

See below for an example of a commented function ready to be added to `cdsr`:

```

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

```