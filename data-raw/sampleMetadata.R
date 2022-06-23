
fname <- system.file(
    'extdata/metadata.csv', package = 'MicrobiomeBenchmarkData'
)
metadata_csv <- utils::read.csv(fname)
url <- metadata_csv[metadata_csv$Title == 'sampleMetadata',]$SourceUrl
sampleMetadata <- readr::read_tsv(url)
usethis::use_data(sampleMetadata, overwrite = TRUE)

