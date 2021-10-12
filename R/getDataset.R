utils::globalVariables(c("sampleMetadata", "dataset", "sample_name"))
#' Get dataset
#'
#' \code{getDataset} imports the datasets as TreeSummarizedExperiment objects.
#'
#' @param x A character vector with the name(s) of the dataset(s). If empty
#' and dryrun = TRUE, returns a message with the names of the available
#' datasets. If empty and dryrun = FALSE, returns a list of
#' TreeSummarizedExperiments with all of the datasets.
#'
#' @param dryrun If TRUE, only returns a message and invisibly returns the
#' names of the datasets as a character vector. If FALSE, it returns the
#' TreeSummarizedExperiment datasets indicated in the argument 'x'.
#'
#' @return A list of TreeSummarizedExperiments.
#'
#' @export
#'
#' @importFrom purrr map2
#' @importFrom utils read.csv
#'
#' @examples
#'
#' datasets_names <- getDataset()
#' datasets <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)
#'
getDataset <- function(x, dryrun = TRUE) {

  metadata_csv <- system.file("extdata/metadata.csv",
                              package = "MicrobiomeBenchmarkData")
  metadata <- utils::read.csv(metadata_csv)
  titles <- sort(metadata[["Title"]])

  if (missing(x)) {

    if (isTRUE(dryrun)) {

      n_titles <- seq_along(titles)
      message(
        paste0(n_titles, " ", titles, collapse = "\n"),
        "\n\nUse vignette('datasets', package = 'MicrobiomeBenchmarkData')",
        " for a detailed description of the datasets.",
        "\n\nUse getDataset(dryrun = FALSE) to import all of the datasets."
      )
      return(invisible(titles))
    } else if (isFALSE(dryrun)) {
      x <- titles
    }
  }

  dataset_names <- x[x %in% titles]

  if (!length(dataset_names))
    stop("No datasets were found for your search.", call. = FALSE)

  dataset_names <- sort(dataset_names)

  if (isTRUE(dryrun)) {

    message(paste0(dataset_names, collapse = "\n"))
    return(invisible(dataset_names))

  } else if (isFALSE(dryrun)) {

    urls <- metadata[match(dataset_names, metadata[["Title"]]), "SourceUrl"]
    output <- purrr::map2(
      dataset_names, urls, ~ .assembleTreeSummarizedExperiment(.x, .y)
    )
    names(output) <- dataset_names
    return(output)
  }

}

#' Assemble TreeSummarizedExperiment
#'
#' \code{.assembleTreeSummarizedExperiment} assembles a TreeSummarizedDataset
#' taking as input the name of the dataset and the URL. This is a helper
#' function for the \code{\link{getDataset}} function.
#'
#' @param dat_name A character string with the name of the dataset.
#' @param dat_url A character string with the URL from Zenodo.
#'
#' @return A TreeSummarizedExperiment
#'
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.table
#' @importFrom dplyr filter
#' @importFrom purrr keep
#'
#' @keywords internal
#'
.assembleTreeSummarizedExperiment <- function(dat_name, dat_url) {

  if (dat_name %in% c("HMP_2012_16S_gingival_V13", "HMP_2012_16S_gingival_V35"))
    dat_name <- "HMP_2012_16S_gingival"

  count_matrix_file <- .getResource(dat_name, dat_url)

  count_matrix <- utils::read.table(
    file = count_matrix_file, header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE
  ) %>%
    as.matrix()

  samples <- colnames(count_matrix)
  col_data <- sampleMetadata %>%
    dplyr::filter(dataset == dat_name, sample_name %in% samples) %>%
    purrr::keep(~all(!is.na(.x)))
  rownames(col_data) <- col_data[["sample_name"]]
  col_data <- col_data[, 3:ncol(col_data)]
  col_data <- col_data[colnames(count_matrix), ]

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = S4Vectors::SimpleList(counts = count_matrix),
    colData = S4Vectors::DataFrame(col_data)
  )

  message("Finished ", dat_name, ".")

  tse

}
