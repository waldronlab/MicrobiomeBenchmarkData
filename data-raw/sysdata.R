
metadata_csv <- system.file(
    "extdata/metadata.csv",
    package = "MicrobiomeBenchmarkData"
)
metadata <- utils::read.csv(metadata_csv)
titles <- sub("_[a-z]+_[a-z]+$", "", metadata$Title)
titles <- unique(titles[titles != "sampleMetadata"])
usethis::use_data(
    metadata, titles,
    internal = TRUE, overwrite = TRUE
)

