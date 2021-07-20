utils::globalVariables(c("Class", "Family", "ps_list_WMS", "ps_list_16S", "function_call"))
# utils::globalVariables(c("body_subsite", "subject_id"))

#' Datasets Used in Calgaro, 2020
#'
#' \code{.calgaro2020Datasets} imports one of the six datasets used in the paper of Calgaro, 2020
#' for benchmarking of methods for differential abundance (DA) analysis on microbiome data. Three
#' of the datasets were generated with 16S rRNA sequencing and come from the V35 dataset in the
#' \code{\link{HMP16SData}} package (v1.2.0). The other three datasets were generated with Whole
#' Metagenome Sequencing and come from the HMP_2012 (Stool_TongueDorsum), Castro-NallarE_2015
#' (Schizophrenia), and ZellerG_2014 (CRC) datasets in the curatedMetagenomicData (v1.12.3)
#' package.
#'
#' Reference:
#' Calgaro, M., Romualdi, C., Waldron, L. et al. Assessment of statistical methods from single
#' cell, bulk RNA-seq, and metagenomics applied to microbiome data. Genome Biol 21, 191 (2020).
#' https://doi.org/10.1186/s13059-020-02104-1
#'
#' @param x
#' A character vector of length 1. Valid options: 16S_Stool_TongueDorsum, 16S_Gingiva_Mucosa,
#' 16S_Subgingival_Supragingival, WMS_Stool_TongueDorsum, WMS_Schizophrenia, WMS_CRC.
#'
#' @return
#' A TreeSummarizedExperiment.
#'
#' @keywords internal
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' library(MicrobiomeBenchmarkData)
#' library(SummarizedExperiment)
#'
#' x <- .calgaro2020Datasets(x = "16S_Subgingival_Supragingival")
#' table(colData(x)$HMP_BODY_SUBSITE)
#'
#' y <- .calgaro2020Datasets(x = "WMS_Stool_TongueDorsum")
#' table(colData(y)$body_subsite)
#'
#' }
#'
.calgaro2020Datasets <- function(x) {

  datasets_16s <- c("16S_Stool_TongueDorsum", "16S_Gingiva_Mucosa", "16S_Subgingival_Supragingival")
  datasets_wms <- c("WMS_Stool_TongueDorsum", "WMS_Schizophrenia", "WMS_CRC")

  valid_datasets <- c(datasets_16s, datasets_wms)
  if (!x %in% valid_datasets)
    stop("Input must be a character string. One of the following: 16S_Stool_TongueDorsum, 16S_Gingiva_Mucosa, 16S_Subgingival_Supragingival, WMS_Stool_TongueDorsum, WMS_Schizophrenia, WMS_CRC.",
         call. = FALSE)

  if (x %in% datasets_16s) {
    fpath <- urls$calgaro2020_16S
    rpath <- .getResourceFromCache(fpath)
    load(rpath)
    dataset_name <- sub("^16S_", "", x)
    ps <- ps_list_16S[[dataset_name]]
    tse <- mia::makeTreeSummarizedExperimentFromphyloseq(ps)
    return(tse)

  } else if (x %in% datasets_wms) {
    fpath <- urls$calgaro2020_WMS
    rpath <- .getResourceFromCache(fpath)
    load(rpath)
    dataset_name <- sub("^WMS_", "", x)
    ps <- ps_list_WMS[[dataset_name]]
    tse <- mia::makeTreeSummarizedExperimentFromphyloseq(ps)
    return(tse)
  }
}

#' Dataset used Beghini, 2019 (nychanesmicrobiome)
#'
#' \code{.beghini2019Nychanesmicrobiome} imports the dataset used in the nychanes paper.
#' A good portion of this code was taken directly from the github repository
#' of the \code{nychnamesmicrobiome} package at
#' \url{https://github.com/waldronlab/nychanesmicrobiome/blob/master/R/loadQiimeData.R}.
#'
#' Reference:
#' Beghini, F., Renson, A., Zolnik, C. P., Geistlinger, L., Usyk, M., Moody,
#' T. U., ... & Waldron, L. (2019). Tobacco exposure associated with oral
#' microbiota oxygen utilization in the New York City Health and Nutrition
#' Examination Study. Annals of epidemiology, 34, 18-25.
#'
#' @return
#' A TreeSummarizedExperiment object
#'
#' @importFrom sas7bdat read.sas7bdat
#' @importFrom phyloseq import_biom
#' @importFrom phyloseq parse_taxonomy_default
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq prune_samples
#' @importFrom phyloseq sample_names
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq subset_taxa
#' @importFrom phyloseq sample_sums
#' @importFrom utils read.delim
#' @importFrom mia makeTreeSummarizedExperimentFromphyloseq
#'
#' @keywords internal
#'
.beghini2019Nychanesmicrobiome <- function() {

  list_of_urls <- list(
    metadata = urls$nychanes_metadata,
    otu_table = urls$nychanes_otu_table,
    taxa_tree = urls$nychanes_taxa_tree,
    rep_set = urls$nychanes_rep_set,
    original_map = urls$nychanes_original_map,
    smokingsampleselection = urls$nychanes_smokingsampleselection
  )

  fpaths <- lapply(list_of_urls, .getResourceFromCache)

  ps <- phyloseq::import_biom(BIOMfilename = fpaths$otu_table,
                              treefilename = fpaths$taxa_tree,
                              refseqfilename = fpaths$rep_set,
                              refseqFunction = phyloseq::parse_taxonomy_default)

  colnames(phyloseq::tax_table(ps)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # remove controls and replicates
  pruned_samples <- c(
    # Controls
    "20151013CZTME1", "20151020CZE3", "20151020CZE4", "20151020TME1",
    "NC1", "NC2",
    # replicates
    "NYDH0036", "NYDH0051", "NYDH0060", "NYDH0152", "NYDH0213", "NYDH0487",
    "NYDH0492", "NYDH0522", "NYDH0527", "NYDH0545R", "NYDH0649c", "NYDH0661",
    "NYDH0691", "NYDH0893", "NYDH0931", "NYDH0988", "NYDH1042", "NYDH1460",
    "NYDH1353"
  )
  ps <- phyloseq::prune_samples(
    !(phyloseq::sample_names(ps) %in% pruned_samples),
    ps
    )

  # metadata block
  original_map <- utils::read.delim(fpaths$original_map)
  metadata <- sas7bdat::read.sas7bdat(fpaths$metadata)
  metadata$smokingstatus <- NULL # drop the smoking status variable that is coded by annotateFullDataset()
  metadata <- dplyr::left_join(original_map, metadata, by = 'KEY')

  phyloseq::sample_names(ps) <- gsub("c|R", "", phyloseq::sample_names(ps))

  sample_selection <- utils::read.delim(fpaths$smokingsampleselection)
  sample_selection$smokingstatus <- factor(as.matrix(sample_selection[,-1]) %*% 1:5,
                                           levels = 1:5,
                                           labels = c("alternativeonly","never","former","secondhand","cigarette"))

  new_metadata_smokingstatus <- dplyr::full_join(metadata, sample_selection, by = c('KEY' = 'key'))
  rownames(new_metadata_smokingstatus) <- new_metadata_smokingstatus$Burklab_ID
  phyloseq::sample_data(ps) <- new_metadata_smokingstatus

  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 1000, ps)

  #Remove OTU classified as chloroplasts and mitochondria
  ps <- phyloseq::subset_taxa(ps, !Class %in% c("D_2__Chloroplast") & !Family %in% c("D_4__Mitochondria"))

  #Merge splitted genera
  splitted_genera <- sapply(data.frame(phyloseq::tax_table(ps)[,"Genus"]), function(x)  grep(" [1-9]",x))
  phyloseq::tax_table(ps)[splitted_genera,"Genus"] <- sapply(phyloseq::tax_table(ps)[splitted_genera,"Genus"], function(x) gsub(" [1-9]","", x))
  drop_prefixes <- function(x) {
    x <- as.character(x)
    x <- strsplit(x,"__")
    x <- sapply(x, `[`, 2)
    x
  }
  phyloseq::tax_table(ps)[,"Phylum"] <- drop_prefixes(phyloseq::tax_table(ps)[,"Phylum"])
  phyloseq::tax_table(ps)[,"Genus"] <- drop_prefixes(phyloseq::tax_table(ps)[,"Genus"] )

  tse <- mia::makeTreeSummarizedExperimentFromphyloseq(ps)

  return(tse)

}
