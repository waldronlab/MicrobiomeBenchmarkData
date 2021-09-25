utils::globalVariables(
  c("Class", "Family", "ps_list_WMS", "ps_list_16S", "StudyID", "study-type",
    "human-readable-description", "function_call", "taxonomy", "submitted_ftp",
    "sra_ftp", "OTU ID")
  )
#' Calgaro 2020 - 16S (subgingival vs supragingival plaque)
#'
#' \code{calgaro_2020_16S} imports the 16S data used in Calgaro 2020 to
#' compare subgingival vs supragingival plaque samples. These samples are a
#' subset taken from the HMP16SData package (v1.2.0).
#'
#' **Biological signal:** anaerobic and aerobic bacteria are under- and over-
#' represented in supragingival plaque samples.
#'
#' @section Reference:
#' Calgaro, M., Romualdi, C., Waldron, L. et al. Assessment of statistical methods from single
#' cell, bulk RNA-seq, and metagenomics applied to microbiome data. Genome Biol 21, 191 (2020).
#' https://doi.org/10.1186/s13059-020-02104-1
#'
#' @keywords internal
#'
#' @return A TreeSummarizedExperiment.
#'
calgaro_2020_16S_gingival_healthy <- function() {

  download_calgaro_2020_16S_gingival_healthy <- function() {
    load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
    ps_list_16S[["Subgingival_Supragingival"]] %>%
      mia::makeTreeSummarizedExperimentFromPhyloseq()
  }

  .getResource(
    resource_name = "calgaro_2020_16S_gingival_healthy",
    FUN = download_calgaro_2020_16S_gingival_healthy
  )
}

#' Calgaro 2020 - WGS cMD3 (subgingival vs supragingival plaque)
#'
#' \code{calgaro_2020_WGS_cMD3} imports the WGS samples used in Calgaro 2020
#' to compare subgingival vs supragingival plaque. This is an updated dataset
#' taken from curatedMetagenomicData v3 instead of curatedMetagenomicData v1.
#'
#' **Biological signal:** anaerobic and aerobic bacteria are under- and over-
#' represented in supragingival plaque samples.
#'
#' @section Reference:
#' Calgaro, M., Romualdi, C., Waldron, L. et al. Assessment of statistical methods from single
#' cell, bulk RNA-seq, and metagenomics applied to microbiome data. Genome Biol 21, 191 (2020).
#' https://doi.org/10.1186/s13059-020-02104-1.
#'
#' @keywords internal
#'
#' @return A TreeSummarizedExperiment.
#'
HMP_WMS_cMD3_gingival_healthy <- function() {

  load_HMP_WMS_cMD3_gingival_healthy <- function() {
    samples <- c("SRS013950", "SRS014107", "SRS014477", "SRS014691", "SRS063215",
                 "SRS014578", "SRS016092", "SRS018665", "SRS023538", "SRS051378")
    hmp <- suppressMessages(curatedMetagenomicData::curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE, counts = TRUE)[[1]])
    hmp_subset <- hmp[,samples]
    mat <- SummarizedExperiment::assay(hmp_subset)
    hmp_subset[apply(mat, 1, function(x) any(x >= 10)),]
  }

  .getResource(
    resource_name = "HMP_WMS_cMD3_gingival_healthy",
    FUN = load_HMP_WMS_cMD3_gingival_healthy
  )
}

#' HMP16SData V13 - healthy gingival plaque
#'
#' \code{HMP_16S_V13_gingival_healthy}
#'
#' @section Reference:
#' TODO
#'
#' @keywords internal
#'
#' @return A TreeSummarizedExperiment.
#'
HMP_16S_V13_gingival_healthy <- function() {

  load_HMP_16S_V13_gingival_healthy <- function() {
    v13se <- HMP16SData::V13()
    v13se <- v13se[,hmp16s_samples[["v13"]]]
    counts <- SummarizedExperiment::assay(v13se)
    v13se[apply(counts, 1, function(x) any(x >= 10)),]

  }

  .getResource(
    resource_name = "HMP_16S_V13_gingival_healthy",
    FUN = load_HMP_16S_V13_gingival_healthy
  )

}

#' HMP16SData V35 - healthy gingival plaque
#'
#' \code{HMP_16S_V35_gingival_healthy}
#'
#' @keywords internal
#'
#' @return A TreeSummarizedExperiment.
HMP_16S_V35_gingival_healthy <- function() {

  load_HMP_16S_V35_gingival_healthy <- function() {
    v35se <- HMP16SData::V35()
    v35se <- v35se[,hmp16s_samples[["v35"]]]
    counts <- SummarizedExperiment::assay(v35se)
    v35se[apply(counts, 1, function(x) any(x >= 10)),]
  }

  .getResource(
    resource_name = "HMP_16S_V35_gingival_healthy",
    FUN = load_HMP_16S_V35_gingival_healthy
  )

}

#' Beghini 2019 (smoking)
#'
#' \code{beghini_2019_16S} imports the dataset used in the nychanes paper.
#' A good portion of this code was taken directly from the github repository
#' of the \code{nychnamesmicrobiome} package at
#' \url{https://github.com/waldronlab/nychanesmicrobiome/blob/master/R/loadQiimeData.R}.
#'
#' @section Reference:
#' Beghini, F., Renson, A., Zolnik, C. P., Geistlinger, L., Usyk, M., Moody,
#' T. U., ... & Waldron, L. (2019). Tobacco exposure associated with oral
#' microbiota oxygen utilization in the New York City Health and Nutrition
#' Examination Study. Annals of epidemiology, 34, 18-25.
#'
#' @keywords internal
#'
#' @return
#' A TreeSummarizedExperiment object.
#'
#'
beghini_2019_16S_saliva_smoking <- function() {

  download_beghini_2019_16S_saliva_smoking <- function() {

    list_of_urls <- list(
      metadata = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/public_v2_010518.sas7bdat",
      otu_table = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/otu_table_mc10_w_tax.biom",
      taxa_tree = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.tre",
      rep_set = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.fna",
      original_map = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/original_map.tsv",
      smokingsampleselection = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/smokingsampleselection.tsv"
    )

    fpaths <- lapply(list_of_urls, function(x) {
      temp_file <- tempfile()
      suppressMessages(utils::download.file(url = x, destfile = temp_file))
      temp_file
    })

    ps <- phyloseq::import_biom(BIOMfilename = fpaths$otu_table,
                                treefilename = fpaths$taxa_tree,
                                refseqfilename = fpaths$rep_set,
                                refseqFunction = phyloseq::parse_taxonomy_default)

    colnames(phyloseq::tax_table(ps)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

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

    ps <- phyloseq::subset_taxa(ps, !Class %in% c("D_2__Chloroplast") & !Family %in% c("D_4__Mitochondria"))

    splitted_genera <- sapply(data.frame(phyloseq::tax_table(ps)[,"Genus"]), function(x)  grep(" [1-9]",x))
    phyloseq::tax_table(ps)[splitted_genera,"Genus"] <- sapply(phyloseq::tax_table(ps)[splitted_genera,"Genus"], function(x) gsub(" [1-9]","", x))
    drop_prefixes <- function(x) {
      as.character(x) %>%
        strsplit("__") %>%
        vapply(function(x) `[`(x, 2), character(1))
    }

    phyloseq::tax_table(ps)[,"Phylum"] <- drop_prefixes(phyloseq::tax_table(ps)[,"Phylum"])
    phyloseq::tax_table(ps)[,"Genus"] <- drop_prefixes(phyloseq::tax_table(ps)[,"Genus"] )

    mia::makeTreeSummarizedExperimentFromphyloseq(ps)

  }

  .getResource(
    resource_name = "beghini_2019_16S_saliva_smoking",
    FUN = download_beghini_2019_16S_saliva_smoking
  )

}

#' Spike-in bacteria dataset
#'
#' \code{spikeInBacteria} provides a dataset with spike-in bacteria data.
#' TODO.
#'
#' @section Reference:
#' Stämmler, F., Gläsner, J., Hiergeist, A. et al. Adjusting microbiome
#' profiles for differences in microbial load by spike-in bacteria.
#' Microbiome 4, 28 (2016). https://doi.org/10.1186/s40168-016-0175-0
#'
#' @keywords internal
#'
#' @return
#' A TreeSummarizedExperiment
#'
spikeInBacteria <- function() {

  download_spikeInBacteria <- function() {

    # Download taxa table (with counts)
    taxa_table_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM10_ESM.txt"
    taxa_table_file <- tempfile()
    utils::download.file(url = taxa_table_url, destfile = taxa_table_file)
    taxa_table <- readr::read_tsv(taxa_table_file, show_col_types = FALSE, progress = FALSE, comment = "#")

    # Extract counts matrix
    counts <- taxa_table %>%
      dplyr::select(-taxonomy) %>%
      tibble::column_to_rownames(var = "OTU ID") %>%
      as.matrix()

    # Extract taxonomy for rowData
    # Most taxa are uncultured, so the taxonomy will be contained in a single table
    row_data <- taxa_table %>%
      dplyr::select(`OTU ID`, taxonomy) %>%
      tibble::column_to_rownames(var = "OTU ID") %>%
      as.data.frame() %>%
      S4Vectors::DataFrame()

    # Download sample metadata
    col_data_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11953&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true"
    col_data_file <- tempfile()
    utils::download.file(url = col_data_url, destfile = col_data_file)
    col_data <- readr::read_tsv(col_data_file, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::filter(grepl("ASCT.MID", submitted_ftp)) %>%
      dplyr::mutate(sample_name = sub("^.+ASCT\\.(MID[0-9]+)_.+$", "\\1", submitted_ftp)) %>%
      dplyr::select(-sra_ftp) %>%
      tibble::column_to_rownames(var = "sample_name") %>%
      as.data.frame() %>%
      S4Vectors::DataFrame()

    # Despite the fact that no tree is provided, the data will be assembled into
    # a TreeSummarizedExperiment object for consistency.
    TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = S4Vectors::SimpleList(counts = counts),
      colData = col_data,
      rowData = row_data
    )

  }

  .getResource(
    resource_name = "spikeInBacteria",
    FUN = download_spikeInBacteria
  )

}
