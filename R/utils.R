utils::globalVariables(c("ps_list_16S", "HMP_BODY_SUBSITE", "SEX", "SAMPLE", "n"))
#' HMP16S Data from Calgaro 2020
#'
#' \code{.hmp16SDataCalgaro2020} imports the HMP16S filtered data from Calgaro, 2020,  which was
#' used for enrichment analysis in that paper.
#'
#' @param x
#' A character vector of length 1; the name of the dataset to be imported.
#'
#' @return
#' A TreeSummarizedExperiment
#'
#' @importFrom mia makeTreeSummarizedExperimentFromphyloseq
#'
#' @import phyloseq
#'
.hmp16SDataCalgaro2020 <- function(x) {
  load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
  ps <- ps_list_16S[[x]]
  tse <- mia::makeTreeSummarizedExperimentFromphyloseq(ps)
  return(tse)
}

#' Subset HMP16SData
#'
#' \code{.hmp16SDataSubset} subsets a dataset, either V13 or V35, from the \code{\link{HMP16SData}}
#' based on common samples and subjects between both datasets.
#'
#' @param se
#' A SummarizedExperiment imported with V13 or V35 functions.
#'
#' @param intersect
#' A list of intersects of samples and subjects (internal data of the \code{MicrobiomeBenchmarkData}).
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors metadata
#'
.hmp16SDataSubset <- function(se, intersect) {

  names(SummarizedExperiment::assays(se))[1] <- "counts"
  SummarizedExperiment::assays(se)$relative_abundance <- apply(SummarizedExperiment::assays(se)$counts, 2, FUN = function(x) x / sum(x) * 100)
  subsites <- c("Subgingival Plaque", "Supragingival Plaque")
  common_samples <- intersect[["samples"]]
  common_subjects <- intersect[["subjects"]]
  se <- se[, colnames(se) %in% common_samples & se$HMP_BODY_SUBSITE %in% subsites & se$VISITNO == 1 & se$RSID %in% common_subjects]
  se <- .sampleSEBySex(se)
  se <- se[rowSums(SummarizedExperiment::assays(se)$counts) > 0,]

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = S4Vectors::SimpleList(counts = SummarizedExperiment::assays(se)$counts,
                                   relative_abundance = SummarizedExperiment::assays(se)$relative_abundance),
    colData = SummarizedExperiment::colData(se),
    rowData = SummarizedExperiment::rowData(se),
    rowTree = S4Vectors::metadata(se)$phylogeneticTree,
    metadata = S4Vectors::SimpleList(experimentData = S4Vectors::metadata(se)$experimentData)
  )

  return(tse)

}

#' Sample By Sex
#'
#' \code{.sampleBySex} is a helper function to randomly select equal number of males and females per
#' HMP_BODY_SUBSITE. It implements \code{set.seed} to obtain the same random samples. Currently,
#' this only works with datasets from the MP16SData package.
#'
#' @param x
#' A SummarizedExperiment from the HMP16SData package
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr count
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_sample
#' @importFrom dplyr ungroup
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#'
#' @return
#' A SummarizedExperiment with equal number of males and females in samples.
#'
#'
.sampleSEBySex <- function(x) {
  min_n <- SummarizedExperiment::colData(x) %>%
    tibble::as_tibble() %>%
    dplyr::count(HMP_BODY_SUBSITE, SEX) %>%
    dplyr::pull(n) %>%
    min()

  set.seed(12345)
  samples <- SummarizedExperiment::colData(x) %>%
    tibble::as_tibble(rownames = "SAMPLE") %>%
    dplyr::group_by(HMP_BODY_SUBSITE, SEX) %>%
    dplyr::slice_sample(n = min_n) %>%
    dplyr::ungroup() %>%
    dplyr::pull(SAMPLE)

  output <- x[, samples]
  return(output)
}

#' HMP16SData Intersects
#'
#' \code{.hmp16SDataIntersects} creates a list of common samples and subjects of the V13 and V35
#' datasets from the \code{HMP16SData} package based on two body subsites.
#'
#' @param subsites
#' A character vector of length 2; the names of the body subsites.
#'
#' @return
#' A list of common samples and subjects between the V13 and V35 datasets.
#'
.hmp16SDataIntersects <- function(subsites) {
  if (!(length(subsites) == 2 & class(subsites) == "character"))
    stop("You must provide the names of two body sites in the HMP_BODY_SUBSITE column of HMP16SData")

  v13 <- HMP16SData::V13()
  v35 <- HMP16SData::V35()

  v13_subjects <- .intersectSubjects(v13, subsites[1], subsites[2], visit = 1)
  v35_subjects <- .intersectSubjects(v35, subsites[1], subsites[2], visit = 1)

  common_subjects <- intersect(v13_subjects, v35_subjects)
  common_samples <- intersect(colnames(v13), colnames(v35))

  output <- list(subjects = as.character(common_subjects), samples = as.character(common_samples))

  return(output)
}

#' Intersect Subjects
#'
#' \code{.intersectSubjects} creates a vector with the intersect of subjects present in
#' two body subsites of an HMP16SData (imported with V13 or V35).
#'
#' @param x
#' A SummarizedExperiment imported by the V13 of V35 functions.#'
#' @param subsite1
#' A character vector of length 1.
#' @param subsite2
#' A character vector of length 1.
#' @param visit
#' An integer; visit number.
#'
#' @return
#' Common subjects between the V13 and V35 datasets.
#'
.intersectSubjects <- function(x, subsite1, subsite2, visit = 1) {

  x1 <- x$RSID[x$HMP_BODY_SUBSITE == subsite1 & x$VISITNO == visit]
  x1 <- x1[!is.na(x1)]

  x2 <- x$RSID[x$HMP_BODY_SUBSITE == subsite2 & x$VISITNO == visit]
  x2 <- x2[!is.na(x2)]

  output <- intersect(x1, x2)
  return(output)
}
