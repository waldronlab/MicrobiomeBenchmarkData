utils::globalVariables(c("dataset_name"))

#' Get datasets
#'
#' \code{getDataSets} imports one or several datasets as
#' TreeSummarizedExperiment(s) stored in a list vector. The datasets
#' can accessed according to their index or name as specified in the output of
#' the \code{\link{listDatasets}} function. The datasets can be loaded
#' according to index or name, but not both.
#'
#' @param index
#' A numeric vector containing the indexes of the datasets to be loaded.
#'
#' @param name
#' A character vector containing the names of the datasets to be loaded.
#'
#' @return
#' A list vector containing TreeSummarizedExperiments.
#'
#' @export
#'
#' @seealso
#' \code{\link{listDatasets}}
#'
#' @examples
#'
#' library(MicrobiomeBenchmarkData)
#'
#' ## Load datasets by index
#' datasets1 <- getDatasets(index = c(1,3))
#' datasets1
#'
#' ## Load datasets by name
#' datasets_names <- c("calgaro_2020_16S_gingival_healthy",
#'                     "beghini_2019_16S_saliva_smoking" )
#' datasets2 <- getDatasets(name = datasets_names)
#' datasets2
#'
getDatasets <- function(index = NULL, name = NULL) {

  if (is.null(index) && is.null(name))
    stop("No argument values provided. You must specified either an 'index' or 'names' argument value.", call. = FALSE)

  if (!is.null(index) && !is.null(name))
    stop("Both 'index' and 'name' arguments provided. You must provide only one of them.", call. = FALSE)

  datasets <- listDatasets()

  if (!is.null(index)) {

    index <- suppressWarnings(as.integer(index))

    if (any(is.na(index)))
      stop("Invalid index value(s) or NAs at position(s) ", paste0(which(is.na(index)), collapse = ", "),
           ". Please remove invalid value(s) or NAs before importing the dataset(s)",
           call. = FALSE)

    if (any(duplicated(index)))
      stop("Index values duplicated at position(s) ", paste0(which(duplicated(index)), collapse = ","),
           ". Please remove duplicated values before importing the dataset(s)",
           call. = FALSE)

    if (any(index < 0) || any(index > nrow(datasets)))
      stop("Index out of boundaries. Thera are ", nrow(datasets), " datasets available.",
           " Please use an integer number between 1 and ", nrow(datasets), " as index argument value.",
           call. = FALSE)

    datasets_names <- datasets[index,][["dataset_name"]]
    output <- lapply(datasets_names, function(x) eval(parse(text = paste0(x, "()")))) %>%
      magrittr::set_names(datasets_names)
    return(output)

  } else if (!is.null(name)) {

    if (!is.character(name))
      stop("'name' argument must be a character vector.", call. = FALSE)

    if (any(is.na(name)))
      stop("NA value(s) at position(s) ", paste0(which(is.na(name)), collapse = ", "),
           ". Please remove all NA values before importing the dataset(s).",
           call. = FALSE)

    if (any(duplicated(name)))
      stop("Names duplicated at position(s) ", paste0(which(duplicated(name)), collapse = ", "),
           ". Please remove duplicated values before importing the dataset(s).",
           call. = FALSE)

    if (any(!name %in% datasets[["dataset_name"]]))
      stop("Dataset(s) ", paste0(name[!name %in% datasets[["dataset_name"]]], collapse = ", "),
           " not found in the list of datasets.",
           " Please check the available datasets with listDatasets().",
           call. = FALSE)

    datasets_names <- datasets %>%
      dplyr::filter(dataset_name %in% name) %>%
      dplyr::pull(dataset_name)
    output <- lapply(datasets_names, function(x) eval(parse(text = paste0(x, "()")))) %>%
      magrittr::set_names(datasets_names)
    return(output)

  }

}

#' List datasets
#'
#' \code{listDatasets} displyas the list of available datasets.
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
#' listDatasets()
#'
listDatasets <- function() {
  datasets_fname <- system.file("extdata/datasets.tsv", package = "MicrobiomeBenchmarkData")
  readr::read_tsv(datasets_fname, show_col_types = FALSE)
}
