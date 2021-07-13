utils::globalVariables(c("body_subsite", "subject_id"))
#' Subgingival vs Supragingival Plaque Dataset 1
#'
#' \code{subgingivalVsSupragingivalDataset1} loads a TreeSummarizedExperiment of Subgingival vs
#' Supragingival Plaque from the HMP16S project used in Calgaro, 2020 for enrichment analysis.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
subgingivalVsSupragingivalDataset1 <- function() {
  output <- .hmp16SDataCalgaro2020("Subgingival_Supragingival")
  return(output)
}

#' Subgingival vs Supragingival Plaque Dataset 2
#'
#' \code{subgingivalVsSupragingivalDataset2} loads a TreeSummarizedExperiment of Subgingival vs
#' Supragingival Plaque from the V13 dataset from the \code{\link{HMP16SData}} package. The samples
#' have been filtered to only include subjects present in both the V13 and V35 datasets of the
#' HMP16SData package and include only the visit number 1 to the research center. Furthermore, the samples
#' include equal numbers of males and females.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#'
subgingivalVsSupragingivalDataset2 <- function() {
  v13 <- HMP16SData::V13()
  v13tse <- .hmp16SDataSubset(v13, subgingival_supragingival_hmp16s_intersects)
  return(v13tse)

}

#' Subgingival vs Supragingival Plaque Dataset 3
#'
#' \code{subgingivalVsSupragingivalDataset3} loads a TreeSummarizedExperiment of Subgingival vs
#' Supragingival Plaque from the V35 dataset from the \code{\link{HMP16SData}} package. The samples
#' have been filtered to only include subjects present in both the V13 and V35 datasets of the
#' HMP16SData package and include only the visit number 1 to the research center. Furthermore, the samples
#' include equal numbers of males and females.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#' @importFrom HMP16SData V35
#'
subgingivalVsSupragingivalDataset3 <- function() {
  v35 <- HMP16SData::V35()
  v35tse <- .hmp16SDataSubset(v35, subgingival_supragingival_hmp16s_intersects)
  return(v35tse)

}

#' Stool vs Tongue Dorsum Dataset 1
#'
#' \code{stoolVsTongueDorsumDataset1} loads a TreeSummarizedExperiment of Stool vs
#' Tongue Dorsum from the HMP16S project used in Calgaro, 2020 for enrichment analysis.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
stoolVsTongueDorsumDataset1 <- function() {
  output <- .hmp16SDataCalgaro2020("Stool_TongueDorsum")
  return(output)
}

#' Stool vs Tongue Dorsum Dataset 2
#'
#' \code{stoolVsTongueDorsumDataset2} loads a TreeSummarizedExperiment of Stool vs Tongue Dorsum
#' from the V35 dataset from the \code{\link{HMP16SData}} package. The samples
#' have been filtered to only include subjects present in both the V13 and V35 datasets of the
#' HMP16SData package and include only the visit number 1 to the research center. Furthermore, the samples
#' include equal numbers of males and females.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#' @importFrom HMP16SData V35
#'
stoolVsTongueDorsumDataset2 <- function() {
  v13 <- HMP16SData::V13()
  v13tse <- .hmp16SDataSubset(v13, stool_tongue_dorsum_hmp16s_intersects)
  return(v13tse)
}

#' Stool vs Tongue Dorsum Dataset 3
#'
#' \code{stoolVsTongueDorsumDataset3} loads a TreeSummarizedExperiment of Stool vs Tongue Dorsum
#' from the V35 dataset from the \code{\link{HMP16SData}} package. The samples
#' have been filtered to only include subjects present in both the V13 and V35 datasets of the
#' HMP16SData package and include only the visit number 1 to the research center. Furthermore, the samples
#' include equal numbers of males and females.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#' @importFrom HMP16SData V35
#'
stoolVsTongueDorsumDataset3 <- function() {
  v35 <- HMP16SData::V35()
  v35tse <- .hmp16SDataSubset(v35, stool_tongue_dorsum_hmp16s_intersects)
  return(v35tse)
}

#' Stool vs Tongue Dorsum Dataset 4
#'
#' \code{stoolVsTongueDorsumDataset4} imports stool and tongue dorsum samples from the HMP_2012
#' study in curatedMetagenomicData (a single sample per subject).
#'
#' @importFrom curatedMetagenomicData curatedMetagenomicData
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_sample
#' @importFrom dplyr ungroup
#' @importFrom dplyr pull
#' @importFrom magrittr %>%
#'
#' @return
#' A TreeSummarizedExperiment
#'
stoolVsTongueDorsumDataset4 <- function() {
  tse <- curatedMetagenomicData::curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE)[[1]]

  relative_abundance <- SummarizedExperiment::assays(tse)[[1]]
  counts <- round(t(t(relative_abundance) * SummarizedExperiment::colData(tse)[["number_reads"]] / 100))
  mode(counts) <- "integer"

  SummarizedExperiment::assays(tse)[[1]] <- counts
  SummarizedExperiment::assays(tse)[[2]] <- relative_abundance
  names(SummarizedExperiment::assays(tse)) <- c("counts", "relative_abundance")

  tse <- tse[,SummarizedExperiment::colData(tse)$body_subsite %in% c("stool", "tongue_dorsum")]

  samples <- SummarizedExperiment::colData(tse) %>%
    tibble::as_tibble(rownames = "sample") %>%
    dplyr::group_by(body_subsite, subject_id) %>%
    dplyr::slice_sample() %>%
    dplyr::ungroup() %>%
    dplyr::pull("sample")

  tse <- tse[rowSums(SummarizedExperiment::assays(tse)$counts) > 0, samples]

  return(tse)

}
