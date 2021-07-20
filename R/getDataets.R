#' Get Data Sets
#'
#' \code{getDataSets} returns a list containing datasets that have been specified by index according
#' to the output of \code{\link{showDatasets}}.
#'
#' @param x
#' A numeric vector containing the indexes of the datasets to be loaded.
#'
#' @return
#' A list of TreeSummarizedExperiments.
#'
#' @export
#'
#' @seealso
#' \code{\link{showDatasets}}
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#' datasets <- getDatasets(x = c(1,3))
#' datasets
#'
getDatasets <- function(x) {
  datasets <- .datasets()
  function_calls <- datasets[x,][["function_call"]]
  list_of_datasets <- lapply(function_calls, function(x) eval(parse(text = x)))
  names(list_of_datasets) <- datasets[x,][["dataset"]]
  return(list_of_datasets)
}

#' Show Datasets
#'
#' \code{showDatasets} shows the list of available datasets in the \code{MicrobiomeBenchmarkData}
#' package.
#'
#' @return
#' A tibble with available datasets.
#'
#' @export
#'
#' @seealso
#' \code{\link{getDatasets}}
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#' showDatasets()
#'
showDatasets <- function() {
  .datasets() %>%
    dplyr::select(-function_call)
}

#' Data Sets
#' \code{dataSets} imports the table containing the list of available datasets in the
#' \code{MicrobiomeBenchmarkData} package.
#'
#' @return
#' A tibble with all available datasets.
#'
#' @importFrom readr read_tsv
#'
#' @keywords internal
#'
.datasets <- function() {
  datasets_fname <- system.file("extdata/datasets.tsv", package = "MicrobiomeBenchmarkData")
  datasets <- readr::read_tsv(datasets_fname, col_types = "dccccc")
  return(datasets)
}

