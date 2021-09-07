## code to prepare `sysdata` dataset goes here

# Packages
library(SummarizedExperiment)
library(dplyr)


# Beghini 2019 -----------------------------------------------------------

urls <- list(
  nychanes_metadata = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/public_v2_010518.sas7bdat",
  nychanes_otu_table = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/otu_table_mc10_w_tax.biom",
  nychanes_taxa_tree = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.tre",
  nychanes_rep_set = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.fna",
  nychanes_original_map = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/original_map.tsv",
  nychanes_smokingsampleselection = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/smokingsampleselection.tsv"
)


# HMP16S samples ---------------------------------------------------------

v13se <- HMP16SData::V13()
v35se <- HMP16SData::V35()

v13se_subgingival <- v13se[,v13se$HMP_BODY_SUBSITE == "Subgingival Plaque"]
V13se_supragingival <- v13se[,v13se$HMP_BODY_SUBSITE == "Supragingival Plaque"]
v35se_subgingival <- v35se[,v35se$HMP_BODY_SUBSITE == "Subgingival Plaque"]
v35se_supragingival <- v35se[,v35se$HMP_BODY_SUBSITE == "Supragingival Plaque"]

list_of_data <- list(colData(v13se_subgingival)$RSID,
                 colData(V13se_supragingival)$RSID,
                 colData(v35se_subgingival)$RSID,
                 colData(v35se_supragingival)$RSID)

intersect_subjects <- purrr::reduce(list_of_data, ~ intersect(.x, .y))

x <- colData(v13se) %>%
  as.data.frame %>%
  as_tibble(rownames = "SAMPLE") %>%
  filter(RSID %in% intersect_subjects, HMP_BODY_SUBSITE %in% c("Subgingival Plaque", "Supragingival Plaque")) %>%
  group_by(RSID, HMP_BODY_SUBSITE) %>%
  slice_sample() %>%
  ungroup()

y <- colData(v35se) %>%
  as.data.frame %>%
  as_tibble(rownames = "SAMPLE") %>%
  filter(RSID %in% intersect_subjects, HMP_BODY_SUBSITE %in% c("Subgingival Plaque", "Supragingival Plaque")) %>%
  group_by(RSID, HMP_BODY_SUBSITE) %>%
  slice_sample() %>%
  ungroup()

hmp16s_samples <- list(v13 = x$SAMPLE, v35 = y$SAMPLE)


# save data --------------------------------------------------------------

usethis::use_data(
  urls,
  hmp16s_samples,
  internal = TRUE, overwrite = TRUE
)
