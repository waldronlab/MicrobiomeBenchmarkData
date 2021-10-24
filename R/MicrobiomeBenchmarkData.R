
#' Get dataset
#'
#' \code{getDataset} imports datasets as TreeSummarizedExperiment objects.
#'
#' @param x A character vector with the name(s) of the dataset(s). If empty
#' and dryrun = TRUE, it returns a message with the names of the available
#' datasets. If empty and dryrun = FALSE, it returns a list of
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
#' @examples
#'
#' datasets_names <- getDataset()
#' datasets <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)
#'
getDataset <- function(x, dryrun = TRUE) {

    if (missing(x)) {

        if (isTRUE(dryrun)) {

            n_titles <- seq_along(titles)
            message(
                paste0(n_titles, " ", titles, collapse = "\n"),
                "\n\nUse",
                " vignette('datasets', package = 'MicrobiomeBenchmarkData')",
                " for a detailed description of the datasets.",
                "\n\nUse getDataset(dryrun = FALSE)",
                " to import all of the datasets."
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

        output <- lapply(dataset_names, .assembleTreeSummarizedExperiment)
        names(output) <- dataset_names
        return(output)
    }

}

#' Remove cache
#'
#' \code{removeCache} removes all files saved in the cache. It will ask for
#' confirmation before removing the cache.
#'
#' @export
#'
#' @return NULL The cache and all of its contents are removed.
#'
#' @examples
#'
#' ## Remove cache
#' \dontrun{
#' removeCache()
#' }
#'
removeCache <- function() {
    cache <- .getCache()
    BiocFileCache::removebfc(cache, ask = interactive())
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
#' @importFrom dplyr filter
#' @importFrom purrr keep
#' @importFrom tibble column_to_rownames
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors SimpleList
#' @importFrom utils read.table
#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang .env
#'
#' @keywords internal
#'
.assembleTreeSummarizedExperiment <- function(x) {

    if (x %in% c("HMP_2012_16S_gingival_V13", "HMP_2012_16S_gingival_V35")) {
        dat_name <- "HMP_2012_16S_gingival"
    } else {
        dat_name <- x
    }

    col_data <- MicrobiomeBenchmarkData::sampleMetadata %>%
        dplyr::filter(
            .data[["dataset"]] == .env[["dat_name"]]
        ) %>%
        purrr::keep(~all(!is.na(.x))) %>%
        tibble::column_to_rownames(var = "sample_name") %>%
        as.data.frame() %>%
        S4Vectors::DataFrame()

    count_matrix <- .getResourcePath(x, "_count_matrix") %>%
        utils::read.table(
            header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
        ) %>%
    as.matrix()

    row_data <- .getResourcePath(x, "_taxonomy_table") %>%
        utils::read.table(
            header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
        ) %>%
    S4Vectors::DataFrame()

    row_tree <- .getResourcePath(x, "_taxonomy_tree") %>%
        ape::read.tree()

    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = count_matrix),
        colData = S4Vectors::DataFrame(col_data),
        rowData = row_data,
        rowTree = row_tree
    )

    message("Finished ", x, ".")

    tse

}

#' Get cache
#'
#' \code{.getCache} creates or loads a cache to store files downloaded through
#' the \code{MicrobiomeBenchmarkData} package.
#'
#' @keywords internal
#'
#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache
#' @return A BiocFileCache object.
#'
.getCache <- function() {
    cache_path <- tools::R_user_dir(
        package = "MicrobiomeBenchmarkData", which = "cache"
    )
    BiocFileCache::BiocFileCache(cache_path, ask = FALSE)
}

#' Get resource path
#'
#' \code{.getResource} downloads the count matrix and store it in the cache.
#'
#' @param resource_name A character string with the name of the dataset.
#' @param resource_url A character string with the URL from Zenodo.
#'
#' @return A character string containing the path to the count matrix in the
#' cache.
#'
#' @importFrom BiocFileCache bfcquery
#' @importFrom BiocFileCache bfcremove
#' @importFrom BiocFileCache bfcadd
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
.getResourcePath <- function(resource, suffix) {

    resource_name <- paste0(resource, suffix)

    resource_url <- metadata %>%
        dplyr::filter(.data[["Title"]] == .env[["resource_name"]]) %>%
        dplyr::pull(.data[["SourceUrl"]])

    cache <- .getCache()

    resources <- BiocFileCache::bfcquery(
        x = cache, query = resource_name, field = "rname", exact = TRUE
    )

    if (nrow(resources) > 1) {

        rids <- dplyr::pull(resources, "rid")
        BiocFileCache::bfcremove(x = cache, rids = rids)
        resource_path <- BiocFileCache::bfcadd(
            x = cache, rname = resource_name, fpath = resource_url,
            download = TRUE
        )
        return(resource_path)

    } else if (nrow(resources) == 0) {

        resource_path <- BiocFileCache::bfcadd(
            x = cache, rname = resource_name, fpath = resource_url,
            download = TRUE
        )
        return(resource_path)

    } else if (nrow(resources) == 1) {

        resource_path <- BiocFileCache::bfcpath(
            x = cache, rids = resources[["rid"]]
        )
        return(resource_path)
    }
}

