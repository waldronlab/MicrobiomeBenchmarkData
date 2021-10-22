# Script to generate sample metadata, count matrix, taxonomy table, and
# taxonomy tree for the HMP WMS data

library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)

samples <- c("SRS013950", "SRS014107", "SRS014477", "SRS014691", "SRS063215",
             "SRS014578", "SRS016092", "SRS018665", "SRS023538", "SRS051378")
tse <- suppressMessages(curatedMetagenomicData::curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE, counts = TRUE)[[1]])
tse_subset <- tse[,samples]
mat <- SummarizedExperiment::assay(tse_subset)
tse_subset <- tse_subset[apply(mat, 1, function(x) any(x >= 10)),]

# Sample metadata ---------------------------------------------------------

col_data <- colData(tse_subset) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_name")

# Count matrix ------------------------------------------------------------

count_matrix <- assay(tse_subset)

# Taxonomy tree -----------------------------------------------------------

row_tree <- rowTree(tse)

# Taxonomy table ----------------------------------------------------------

row_data <- rowData(tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "TAXA")

# Export files ------------------------------------------------------------

readr::write_tsv(
    x = col_data,
    file = "inst/extdata/HMP_2012_WMS_gingival_sample_metadata.tsv"
)

write.table(
    x = count_matrix, sep = "\t", row.names = TRUE, col.names = TRUE,
    quote = TRUE,
    file = "inst/extdata/HMP_2012_WMS_gingival_count_matrix.tsv"
)

ape::write.tree(
    phy = row_tree,
    file = "inst/extdata/HMP_2012_WMS_gingival_taxonomy_tree.newick"
)

readr::write_tsv(
    x = row_data,
    file = "inst/extdata/HMP_2012_WMS_gingival_taxonomy_table.tsv"
)
