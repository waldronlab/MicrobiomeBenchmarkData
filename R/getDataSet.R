#' Get Data Sets
#'
#' \code{getDataSets} returns a list containing datasets that have been specified by index according
#' to the output of \code{\link{dataSets}}.
#'
#' @param x
#' A numeric vector containing the indexes of the datasets to be loaded.
#'
#' @return
#' A list of TreeSummarizedExperiments.
#'
#' @export
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#' datasets <- getDataSets(x = 1)
#' datasets
#'
getDataSets <- function(x) {
  datasets <- dataSets()
  functions <- datasets[x,][["function_name"]]
  list_of_datasets <- lapply(functions, function(x) eval(call(x)))
  names(list_of_datasets) <- datasets[x,][["dataset"]]
  return(list_of_datasets)
}

#' Data Sets
#' \code{dataSets} returns the available datasets in the  `MicrobiomeBenchmarkData` package.
#'
#' @return
#' A tibble with all available datasets.
#'
#' @importFrom readr read_tsv
#'
#' @export
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#' datasets <- dataSets()
#'
dataSets <- function() {
  datasets_fname <- system.file("extdata/datasets.tsv", package = "MicrobiomeBenchmarkData")
  datasets <- readr::read_tsv(datasets_fname, col_types = "dccccc")
  return(datasets)
}
