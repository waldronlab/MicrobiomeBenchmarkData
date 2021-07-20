#' Get MicrobiomeBenchmarkData's Cache
#'
#' \code{.getCache} creates or loads a cache to store files downloaded through the
#' \code{MicrobiomeBenchmarkData} package.
#'
#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache
#'
#' @keywords internal
#'
.getCache <- function() {
  cache_path <- tools::R_user_dir(package = "MicrobiomeBenchmarkData", which = "cache")
  cache <- BiocFileCache::BiocFileCache(cache = cache_path, ask = FALSE)
  return(cache)
}

#' Get Resource From Cache
#'
#' \code{getResourceFromCache} gets the path to a resource in the cache generated with
#' \code{\link{.getCache}}.
#'
#' \code{.getResourceFromCache} first searches for a resource in the cache. If the resource is not
#' found, then the resource is downloaded from the web to the cache.
#'
#' @param fpath
#' A character vector of length 1 indicating the full path to the resource that must be downloaded.
#'
#' @return
#' A character vector of length 1 containing the full path to the resource in the cache.
#'
#' @keywords internal
#'
.getResourceFromCache <- function(fpath) {
  cache <- .getCache()
  resources <- BiocFileCache::bfcquery(x = cache, query = fpath, field = "fpath", exact = TRUE)
  if (nrow(resources) == 0) {
    resource_path <- BiocFileCache::bfcadd(x = cache, rname = fpath, fpath = fpath, download = TRUE)
  } else if (nrow(resources) == 1) {
    resource_path <- resources[["rpath"]]
  }
  return(resource_path)
}
