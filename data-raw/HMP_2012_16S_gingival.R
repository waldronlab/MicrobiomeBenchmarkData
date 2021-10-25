
# This script is to generate the sample metadata, count matrix, taxonomy table,
# and taxonomy tree of the HMP_2012_16S_gingival_V35 and the
# HMP_2012_16S_gingival_V13 datasets.

library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(S4Vectors)
library(dplyr)

v35se <- HMP16SData::V35()
v13se <- HMP16SData::V13()

v35tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = assay(v35se)),
    colData = colData(v35se),
    rowData = rowData(v35se),
    rowTree = metadata(v35se)[["phylogeneticTree"]]
)

v13tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = assay(v13se)),
    colData = colData(v13se),
    rowData = rowData(v13se),
    rowTree = metadata(v13se)[["phylogeneticTree"]]
)

rm(v35se, v13se) # free memory
gc()

intersect_samples <- intersect(colnames(v35tse), colnames(v13tse))

# Sample metadata ---------------------------------------------------------

## Sample metadata must be identical for both datasets since the purpose is
## to provide the sample samples with V13 and V35

select_samples <- c("Subgingival Plaque", "Supragingival Plaque")

v35_col_data <- colData(v35tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_name") %>%
    filter(HMP_BODY_SUBSITE %in% select_samples,
           sample_name %in% intersect_samples,
           # RUN_CENTER == "WUGC",
           # one run center might be incorrectly annotated as "0"
           RUN_CENTER != "0",
           !is.na(SRS_SAMPLE_ID)) %>%
    group_by(HMP_BODY_SUBSITE, RSID) %>%
    slice_min(RSID)

v13_col_data <- colData(v13tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_name") %>%
    filter(HMP_BODY_SUBSITE %in% select_samples,
           sample_name %in% intersect_samples,
           # RUN_CENTER == "WUGC",
           # one run center might be incorrectly annotated as "0"
           RUN_CENTER != "0",
           !is.na(SRS_SAMPLE_ID)) %>%
    group_by(HMP_BODY_SUBSITE, RSID) %>%
    slice_min(RSID)

if (nrow(v35_col_data) > nrow(v13_col_data)) {

    v35_col_data <-
        v35_col_data[v35_col_data$sample_name %in% v13_col_data$sample_name,]

} else if (nrow(v13_col_data) < nrow(v35_col_data)) {

    v13_col_data <-
        v13_col_data[v13_col_data$sample_name %in% v35_col_data$sample_name,]

}

if (identical(v13_col_data, v35_col_data)) {
    hmp_gingival_samples <-
        intersect(v13_col_data$sample_name, v35_col_data$sample_name)
    # generate sample metadata
    col_data <- v13_col_data # Either v13 or v35 works
} else {
    stop("Metadata are not identical.")
}

# TSE subsets -------------------------------------------------------------

v35_subset <- v35tse[,hmp_gingival_samples]
v35_subset <- v35_subset[rowSums(assay(v35_subset)) > 0, ]

v13_subset <- v13tse[,hmp_gingival_samples]
v13_subset <- v13_subset[rowSums(assay(v35_subset)) > 0, ]

rm(v35tse, v13tse)
gc()

# Count matrix ------------------------------------------------------------

v35_count_matrix <- assay(v35_subset)
v13_count_matrix <- assay(v13_subset)

# Taxonomy tree -----------------------------------------------------------

v35_row_tree <- rowTree(v35_subset)
v13_row_tree <- rowTree(v13_subset)

# Taxonomy table ----------------------------------------------------------

v35_row_data <- rowData(v35_subset) %>%
    as.data.frame() %>%
    as_tibble(rownames = "TAXA")

v13_row_data <- rowData(v13_subset) %>%
    as.data.frame() %>%
    as_tibble(rownames = "TAXA")


# Export files ------------------------------------------------------------

## Save sample metadata
readr::write_tsv(
    x = col_data,
    file = "HMP_2012_16S_gingival_sample_metadata.tsv"
)

## Save count matrix
write.table(
    x = v35_count_matrix,
    file = "HMP_2012_16S_gingival_V35_count_matrix.tsv",
    sep = "\t", row.names = TRUE, col.names = TRUE
)

write.table(
    x = v13_count_matrix,
    file = "HMP_2012_16S_gingival_V13_count_matrix.tsv",
    sep = "\t", row.names = TRUE, col.names = TRUE
)

## Save phylogenetic tree
ape::write.tree(
    phy = v35_row_tree,
    file = "HMP_2012_16S_gingival_V35_taxonomy_tree.newick"
)
ape::write.tree(
    phy = v13_row_tree,
    file = "HMP_2012_16S_gingival_V13_taxonomy_tree.newick"
)

## Save taxonomy table
readr::write_tsv(
    x = v35_row_data,
    file = "HMP_2012_16S_gingival_V35_taxonomy_table.tsv"
)

readr::write_tsv(
    x = v13_row_data,
    file = "HMP_2012_16S_gingival_V13_taxonomy_table.tsv"
)
