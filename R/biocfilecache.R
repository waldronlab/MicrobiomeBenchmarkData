#' Get cache for the MicrobiomeBenchmarkData package
#'
#' \code{.getCache} creates or loads a cache to store files downloaded through the
#' \code{MicrobiomeBenchmarkData} package.
#'
#' @keywords internal
#'
.getCache <- function() {
  tools::R_user_dir(package = "MicrobiomeBenchmarkData", which = "cache") %>%
    BiocFileCache::BiocFileCache(ask = FALSE)
}

#' Get resource path from the cache
#'
#' \code{getResourcePathFromCache} gets the path to a resource in the cache generated with
#' \code{\link{.getCache}}.
#'
#' @param resource_name
#' A character string indicating the name of a resource.
#'
#' @return
#' A character string indicating the full path to the resource in the cache.
#'
#' @keywords internal
#'
.getResourcePathFromCache <- function(resource_name) {
  cache <- .getCache()
  resources <- BiocFileCache::bfcquery(cache, query = resource_name, field = "rname", exact = TRUE)
  if (nrow(resources) > 1) {
    rids <- dplyr::pull(resources, "rid")
    BiocFileCache::bfcremove(cache, rids)
    resource_path <- BiocFileCache::bfcnew(cache, rname = resource_name, ext = ".rds", rtype = "local")
    return(resource_path)
  } else if (nrow(resources) == 0) {
    resource_path <- BiocFileCache::bfcnew(cache, rname = resource_name, ext = ".rds", rtype = "local")
    return(resource_path)
  } else if (nrow(resources) == 1) {
    resource_path <- BiocFileCache::bfcpath(cache, resources[["rid"]])
    return(resource_path)
  }
}

#' Get resource
#'
#' \code{.getResource} loads a resource (dataset) from the
#' `MicrobiomeBenchmarkData`'s cache; if the resource doesn't exist in the cache,
#' it will be downloaded first from the internet with a helper function
#' specified with the `FUN` argument.
#'
#' The \code{.getResource} function is meant to be used inside the functions
#' used to download the datasets, e.g.
#' \code{\link{calgaro_2020_16S_gingival_healthy}}.
#'
#' @param resource_name A single character string with the name of the
#' dataset.
#' @param FUN Function to download the dataset.
#'
#' @return A TreeSummarizedExperiment.
#'
#' @keywords internal
#'
.getResource <- function(resource_name, FUN) {
  resource_path <- .getResourcePathFromCache(resource_name)
  if (isTRUE(file.exists(resource_path))) {
    tse <- readRDS(resource_path)
  } else if (isFALSE(file.exists(resource_path))) {
    message("Resource ", resource_name, " not in cache. Downloading and adding to cache.")
    tse <- FUN()
    saveRDS(tse, file = resource_path)
  }
  message("Loading ", resource_name, " from cache.")
  return(tse)
}
