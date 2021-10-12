#' Get cache
#'
#' \code{.getCache} creates or loads a cache to store files downloaded through the
#' \code{MicrobiomeBenchmarkData} package.
#'
#' @keywords internal
#'
#' @importFrom tools R_user_dir
#' @importFrom BiocFileCache BiocFileCache
#' @return A BiocFileCache object.
#'
.getCache <- function() {
  tools::R_user_dir(package = "MicrobiomeBenchmarkData", which = "cache") %>%
    BiocFileCache::BiocFileCache(ask = FALSE)
}

#' Get resource
#'
#' \code{.getResource} downloads the count matrix and store it in the cache.
#'
#' @param dat_name A character string with the name of the dataset.
#' @param dat_url A character string with the URL from Zenodo.
#'
#' @return A character string containing the path to the count matrix in the
#' cache.
#'
#' @importFrom BiocFileCache bfcquery
#' @importFrom BiocFileCache bfcremove
#' @importFrom BiocFileCache bfcadd
#' @importFrom dplyr pull
#'
#' @keywords internal
#'
.getResource <- function(dat_name, dat_url) {

  resource_name <- paste0(dat_name, "_count_matrix.tsv")

  cache <- .getCache()

  resources <- BiocFileCache::bfcquery(
    x = cache, query = resource_name, field = "rname", exact = TRUE
  )

  if (nrow(resources) > 1) {

    rids <- dplyr::pull(resources, "rid")
    BiocFileCache::bfcremove(cache, rids)
    resource_path <- BiocFileCache::bfcadd(
      x = cache, rname = resource_name, fpath = dat_url, download = TRUE
    )
    return(resource_path)

  } else if (nrow(resources) == 0) {

    resource_path <- BiocFileCache::bfcadd(
      x = cache, rname = resource_name, fpath = dat_url, download = TRUE
    )
    return(resource_path)

  } else if (nrow(resources) == 1) {

    resource_path <- BiocFileCache::bfcpath(cache, resources[["rid"]])

    return(resource_path)
  }
}

#' Remove cache
#'
#' \code{removeCache} removes all files saved in the cache. It will ask for
#' confirmation before removing the cache.
#'
#' @export
#'
#' @return The cache and all of its contents are removed.
#'
#' @examples
#'
#' ## Remove cache
#' removeCache()
#'
removeCache <- function() {
    cache <- .getCache()
    BiocFileCache::removebfc(cache, ask = interactive())
}
