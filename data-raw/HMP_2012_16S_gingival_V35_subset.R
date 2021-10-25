
## Script to generate the sample metadata, count matrix, taxonomy table, and
## taxonomy tree of the subset of the HMP 2012 data used in Calgaro 2020.

library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)

import_calgaro_2020 <- function() {
    load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
    mia::makeTreeSummarizedExperimentFromPhyloseq(
        ps_list_16S[["Subgingival_Supragingival"]]
    )
}

tse <- import_calgaro_2020()

# Sample metadata ---------------------------------------------------------

col_data <- colData(tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_name")

# Count matrix ------------------------------------------------------------

count_matrix <- assay(tse)

# Taxonomy tree -----------------------------------------------------------

row_tree <- rowTree(tse)

# Taxonomy table ----------------------------------------------------------

fname <- "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/biosis.tsv"
biosis <- readr::read_tsv(fname, col_names = c("GENUS", "BIOSIS"),
    col_types = "cc")

row_data <-
    rowData(tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "TAXA_NAME") %>%
    left_join(biosis, by = "GENUS")

# Export files ------------------------------------------------------------

readr::write_tsv(x = col_data,
                 file = "HMP_2012_16S_gingival_V35_subset_sample_metadata.tsv")

write.table(x = count_matrix,
            file = "HMP_2012_16S_gingival_V35_subset_count_matrix.tsv",
            sep = "\t", row.names = TRUE, col.names = TRUE)

ape::write.tree(
    row_tree,
    file = "HMP_2012_16S_gingival_V35_subset_taxonomy_tree.newick"
)

readr::write_tsv(
    x = row_data,
    file = "HMP_2012_16S_gingival_V35_subset_taxonomy_table.tsv"
)
