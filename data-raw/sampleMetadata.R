## code to prepare `sampleMetadata` dataset goes here

# HMP_2012_16S_gingival_sample_metadata <-
#   read.table(file = "inst/extdata/HMP_2012_16S_gingival_sample_metadata.tsv",
#              sep = "\t", header = TRUE)

# sampleMetadata <- HMP_2012_16S_gingival_sample_metadata

# sampleMetadata <- read.table(
#   file = "inst/extdata/sampleMetadata.tsv", sep = "\t", header = TRUE,
#   row.names = FALSE, check.names = FALSE, comment.char = "")

sampleMetadata <- readr::read_tsv(file = "inst/extdata/sampleMetadata.tsv") %>%
  as.data.frame()

usethis::use_data(sampleMetadata, overwrite = TRUE)
