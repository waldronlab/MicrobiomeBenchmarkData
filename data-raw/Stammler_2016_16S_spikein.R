
## This scripts contains code to generate sample metadata, count matrix,
## taxonomy table, and taxonomy tree of the Stammler 2016 dataset (using
## spike-in bacteria)

## Execute this scritp as an independent file (not in an R package session)

library(magrittr)
library(S4Vectors)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidyr)
library(purrr)
library(taxizedb)

# Sample metadata ---------------------------------------------------------

## The sample metadata is obtained by combining the sample information available
## in EBI and the information provided by the article.

ebi_metadata_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11953&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true"
article_metadata_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM8_ESM.txt"

ebi_metadata <- readr::read_tsv(
    ebi_metadata_url, show_col_types = FALSE, progress = FALSE
    ) %>%
    filter(grepl("ASCT.MID", submitted_ftp)) %>%
    mutate(
        sample_name = sub("^.+ASCT\\.(MID[0-9]+)_.+$", "\\1", submitted_ftp)
    ) %>%
    select(-sra_ftp) %>%
    relocate(sample_name)

article_metadata <- readr::read_tsv(
    article_metadata_url, show_col_types = FALSE, progress = FALSE
)

col_data <- left_join(
    article_metadata, ebi_metadata, by = c("SampleID" = "sample_name")
) %>%
    dplyr::rename(sample_name = SampleID) %>%
    set_colnames(., tolower(colnames(.)))


# Count matrix ------------------------------------------------------------

## Import data (count_matrix/otu_table and taxonomy)
data_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM10_ESM.txt"
data <- utils::read.table(
    file = data_url, header = TRUE, sep = "\t", row.names = 1,
    comment.char = "#", check.names = FALSE
    )

## Count matrix / OTU table
count_matrix <- as.matrix(data[,colnames(data) != "taxonomy"])

# Taxonomy table ----------------------------------------------------------

row_data <- tibble::tibble(
    TAXA = rownames(data), taxonomy = data[["taxonomy"]]
)


# output <- vector("list", nrow(taxonomy))
# for (i in seq_along(output)) {
#     output[i] <- stringr::str_split(taxonomy[i, "taxonomy"], "; __")
# }
# table(vapply(output, length, integer(1)))
#
# vapply(output, \(x) tail(x, 1), character(1))
#
# taxonomy[vapply(output, length, integer(1)) == 6, ] %>% View()
# as.data.frame(table(flatten_chr(output))) %>% View()


# Test that things work ---------------------------------------------------
rowTree <- NULL

colData <- col_data %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    as.data.frame() %>%
    DataFrame()

rowData <- row_data %>%
    tibble::column_to_rownames(var = "TAXA") %>%
    as.data.frame() %>%
    DataFrame()

tse <- TreeSummarizedExperiment(
    assays = SimpleList(count_matrix),
    colData = colData,
    rowData = rowData,
    rowTree = rowTree
)

# Export files ------------------------------------------------------------

## Sample metadata
readr::write_tsv(
    x = col_data, file = "Stammler_2016_16S_spikein_sample_metadata.tsv"
)

## Count matrix
write.table(
    x = count_matrix,
    sep = "\t", row.names = TRUE, col.names = TRUE, quote = TRUE,
    file = "Stammler_2016_16S_spikein_count_matrix.tsv"
)

## Taxonomy table
readr::write_tsv(
    x = row_data, file = "Stammler_2016_16S_spikein_taxonomy_table.tsv"
)
