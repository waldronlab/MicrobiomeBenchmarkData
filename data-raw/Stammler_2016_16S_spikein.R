# "Stammler_2016_16S_spikein"

library(dplyr)
library(magrittr)
library(S4Vectors)
library(SummarizedExperiment)

## Import data (count_matrix/otu_table and taxonomy)
data_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM10_ESM.txt"
data <- utils::read.table(
    file = data_url, header = TRUE, sep = "\t", row.names = 1,
    comment.char = "#", check.names = FALSE
    )

## Taxonomy
taxonomy <- data.frame(OTU = rownames(data), taxonomy = data[["taxonomy"]])

## Count matrix / OTU table
count_matrix <- as.matrix(data[,colnames(data) != "taxonomy"])

## Sample metadata
ebi_metadata_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11953&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true"
article_metadata_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM8_ESM.txt"

ebi_metadata <- readr::read_tsv(ebi_metadata_url, show_col_types = FALSE, progress = FALSE) %>%
    filter(grepl("ASCT.MID", submitted_ftp)) %>%
    mutate(sample_name = sub("^.+ASCT\\.(MID[0-9]+)_.+$", "\\1", submitted_ftp)) %>%
    select(-sra_ftp) %>%
    relocate(sample_name)

article_metadata <- readr::read_tsv(article_metadata_url, show_col_types = FALSE, progress = FALSE)

sample_metadata <- left_join(article_metadata, ebi_metadata, by = c("SampleID" = "sample_name")) %>%
    rename(sample_name = SampleID) %>%
    set_colnames(., tolower(colnames(.)))

col_data <- sample_metadata %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    as.data.frame() %>%
    DataFrame()

se <- SummarizedExperiment(
    assays = SimpleList(count_matrix),
    colData = col_data
)


readr::write_tsv(sample_metadata, "inst/extdata/Stammler_2016_16S_spikein_sample_metadata.tsv")
write.table(x = count_matrix, sep = "\t", row.names = TRUE, col.names = TRUE, quote = TRUE,
            file = "inst/extdata/Stammler_2016_16S_spikein_count_matrix.tsv")

