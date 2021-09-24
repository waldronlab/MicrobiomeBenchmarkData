utils::globalVariables(c("dataset_name", "."))
#' Get datasets
#'
#' \code{getDataSets} imports one or several datasets as
#' TreeSummarizedExperiment(s) stored in a list vector. The datasets
#' can accessed according to their index or name as specified in the output of
#' the \code{\link{.listDatasets}} function. The datasets can be loaded
#' according to index or name, but not both.
#'
#' @param dat_name
#' A character vector containing the names of the datasets to be loaded. If
#' 'NULL' (default), a message indicating how to get detailed information
#' about the datasets is printed along with the names of available datasets;
#' furthermore, the list of available datasets is invisibly returned as a
#' character vector.
#'
#' @return
#' A list containing TreeSummarizedExperiments. If `NULL`, a character vector
#' with available datasets (invisibly returned).
#'
#' @export
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#'
#' ## Show datasets and save datasets names as a character vector
#' datasets <- getDatasets()
#'
#' ## Download datasets 1 and two in the list
#' x <- getDatasets(datasets[1:2])
#'
#' length(x)
#' vapply(x, class, character(1))
#'
getDatasets <- function(dat_name = NULL) {

  datasets_names <- .listDatasets() %>%
    dplyr::arrange(dataset_name) %>%
    dplyr::pull(dataset_name)

  if (is.null(dat_name)) {

    message(
            "Please use ",
            "vignette('datasets', package = 'MicrobiomeBenchmarkData')",
            " to get a detailed description of the datasets included in the",
            " MicrobiomeBenchmarkData package.",
            " Currently, there are ", length(datasets_names),
            " datasets alvailable: ", "\n\n",
            paste0(paste0(1:length(datasets_names), " ", datasets_names), collapse = "\n")
    )
    return(invisible(datasets_names))
  }

  if (!is.character(dat_name))
    stop("Invalid argument type. Argument 'dat_name' must be a character vector, not an object of class '", class(dat_name), "'.", call. = FALSE)

  if (any(is.na(dat_name))) {
    na_pos <- which(is.na(dat_name))
    stop("NA values. Argument 'dat_name' must not contain NA values. Remove NA values at",
         " position(s) ", paste0(na_pos, collapse = ","), call. = FALSE)
  }

  if (any(duplicated(dat_name))) {
    dup_pos <- which(duplicated(dat_name))
    stop("Duplicate values. Remove duplicate values at position(s) ", paste0(dup_pos, collapse = ","), ".", call. = FALSE)
  }

  if (any(!dat_name %in% datasets_names)) {
    no_name <- dat_name[which(!dat_name %in% datasets_names)] %>%
      paste0("'", ., "'")
    stop("Invalid dataset name. Dataset(s) ", paste0(no_name, collapse = ","),
         " not found in MicrobiomeBenchmarkData. See the list of alvailable datasets with getDatasets().", call. = FALSE)
  }

  selected_datasets <- datasets_names[datasets_names %in% dat_name]

  selected_datasets %>%
    lapply(function(x) eval(parse(text = paste0(x, "()")))) %>%
      magrittr::set_names(selected_datasets)
}

#' List datasets
#'
#' \code{.listDatasets} displays the list of available datasets.
#'
#' @return
#' A tibble with available datasets.
#'
#' @keywords internal
#'
.listDatasets <- function() {
  datasets_fname <- system.file("extdata/datasets.tsv", package = "MicrobiomeBenchmarkData")
  readr::read_tsv(datasets_fname, show_col_types = FALSE)
}
