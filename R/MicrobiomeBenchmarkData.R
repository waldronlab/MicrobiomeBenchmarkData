
#' Get dataset
#'
#' \code{getBenchmarkData} imports datasets as TreeSummarizedExperiment objects.
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
#' ## Example 1
#' datasets_names <- getBenchmarkData()
#' datasets_names
#'
#' ## Example 2
#' dataset <- getBenchmarkData("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)
#' dataset[[1]]
#'
getBenchmarkData <- function(x, dryrun = TRUE) {
    if (missing(x)) {
        if (isTRUE(dryrun)) {
            n_titles <- seq_along(titles)
            message(
                paste0(n_titles, " ", titles, collapse = "\n"),
                "\n\nUse",
                " vignette('datasets', package = 'MicrobiomeBenchmarkData')",
                " for a detailed description of the datasets.",
                "\n\nUse getBenchmarkData(dryrun = FALSE)",
                " to import all of the datasets."
            )
            return(invisible(titles))
        } else if (isFALSE(dryrun)) {
            x <- titles
        }
    }

    dataset_names <- x[x %in% titles]

    if (!length(dataset_names)) {
        stop("No datasets were found for your search.", call. = FALSE)
    }

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
#' \code{removeCache} removes all files saved in the cache.
#'
#' @param ask If TRUE, a prompt will appear asking the user to confirm removal
#' of cache. Default value is given by the \code{interactive} function.
#'
#' @export
#'
#' @return NULL The cache and all of its contents are removed.
#'
#' @examples
#'
#' ## Remove cache
#' removeCache()
#'
removeCache <- function(ask = interactive()) {

    cache <- .getCache()
    cache_info <- BiocFileCache::bfcinfo(cache)

    prompt_msg <- paste0(
        "Remove cache and ", nrow(cache_info), " resources?",
        " (yes/no): "
    )

    if (ask) {
        answer <- readline(prompt = prompt_msg)
        if (answer == 'yes') {
            message('Removing cache.')
            BiocFileCache::removebfc(cache, ask = FALSE)
            return(invisible(NULL))
        } else if (answer == 'no') {
            message('Cache was not removed.')
            return(invisible(NULL))
        } else {
            message('Not a valid option. Please enter yes or no.')
            return(invisible(NULL))
        }
    } else {
        BiocFileCache::removebfc(cache, ask = FALSE)
        return(invisible(NULL))
    }
}


#' Assemble TreeSummarizedExperiment
#'
#' \code{.assembleTreeSummarizedExperiment} assembles a TreeSummarizedDataset
#' taking as input the name of the dataset and the URL. This is a helper
#' function for the \code{\link{getBenchmarkData}} function.
#'
#' @param dat_name A character string with the name of the dataset.
#' @param dat_url A character string with the URL from Zenodo.
#'
#' @return A TreeSummarizedExperiment
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors SimpleList
#' @importFrom utils read.table
#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
#' @keywords internal
#'
.assembleTreeSummarizedExperiment <- function(x) {

    dat_name <- x

    col_data <- MicrobiomeBenchmarkData::sampleMetadata |>
        {\(y)  y[y$dataset == dat_name, ]}() |>
        {\(y) y[,vapply(y, \(x) !all(is.na(x)), logical(1)), drop = FALSE]}() |>
        S4Vectors::DataFrame()
    rownames(col_data) <- col_data$sample_id
    col_data <- col_data[, colnames(col_data) != "sample_id"]

    count_matrix <- .getResourcePath(x, "_count_matrix") |>
        utils::read.table(
            header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
            quote = ""
        ) |>
        as.matrix()

    row_data <- .getResourcePath(x, "_taxonomy_table") |>
        utils::read.table(
            header = TRUE, row.names = 1, sep = "\t", check.names = FALSE,
            quote = ""
        ) |>
        S4Vectors::DataFrame()

    row_tree_path <- .getResourcePath(x, "_taxonomy_tree")

    if (!length(row_tree_path)) {
        row_tree <- NULL
    } else {
        row_tree <- ape::read.tree(row_tree_path)
    }

    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = count_matrix),
        colData = col_data,
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
#'
#' @keywords internal
#'
.getResourcePath <- function(resource, suffix) {
    ## Inspiration for this code: https://github.com/Bioconductor/AnnotationForge/blob/3e01d4f3620396578eee7783849589fe8d4b23aa/vignettes/MakingNewAnnotationPackages.Rnw#L196-L201
    resource_name <- paste0(resource, suffix)
    resource_url <- metadata[metadata$Title == resource_name,]$SourceUrl
    if (!length(resource_url)) {
        warning(
            "No ", sub("^_", "", suffix), " available for ", resource, ".",
            call. = FALSE
        )
        return(NULL)
    }
    cache <- .getCache()
    BiocFileCache::bfcrpath(
        x = cache, rname = resource_url, exact = TRUE,
        download = TRUE, rtype = "web"
    )
}
