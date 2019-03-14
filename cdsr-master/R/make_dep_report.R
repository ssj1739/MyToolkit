
#' Create html report assessing genetic and chemical dependencies related to specified cell line features
#'
#' @param CL_var vector of cell line features. Can be either binary/logical specifying two groups, or continuous values. Must be named with CCLE_IDs
#' @param covars Optional named vector of covariates to be 'controlled for' in analysis
#' @param test_dir Direction of association testing. Default is left-tailed tests, which corresponds to stronger (more negative) dependencies when e.g. the 'in-group' is of interest in a two-group comparison. Options are ['left', 'right', and 'both']
#' @param report_title Name of report
#' @param used_gscs Name of gene set collections to use. Defaults to include GO_biological_process, canonical (from MSIGDB) and CORUM complex co-membership
#' @param out_dir Directory to generate report
#' @param save_pref_deps Should the table of preferential dep stats be saved to an rds file (default FALSE)
#' @param nproc Number of processes to use (0 for no multithreading). Default 0. Only affects gsa calcs with covariates at this stage
#'
#' @return
#' @export
#'
#' @examples
make_dep_report <- function(
  CL_var, covars = NULL, test_dir = 'left', report_title = 'dep_report',
  used_gscs = c('GO_biological_process', 'canonical', 'CORUM'),
  out_dir = NULL, save_pref_deps = FALSE, nproc = 0) {
    library(rmarkdown)
    if (is.null(out_dir)) {
      out_dir <- getwd()
    }
    if (report_title == 'dep_report') {
      report_full_title <- 'Dependency report'
    } else {
      report_full_title <- paste0(report_title, ' dependency report')
    }
  if (save_pref_deps) {
    pref_deps_file <- file.path(out_dir, sprintf('%s_pref_deps.rds', report_title))
  } else {
    pref_deps_file <- NULL
  }
    rmarkdown::render(input = system.file("rmd/depCon_template.Rmd", package="cdsr"),
                      output_file = file.path(out_dir, paste0(report_title, '.html')),
                      params = list(CL_var = CL_var,
                                    covars = covars,
                                    report_title = report_full_title,
                                    used_gscs = used_gscs,
                                    test_dir = test_dir,
                                    pref_deps_file = pref_deps_file,
                                    nproc = nproc))
}
