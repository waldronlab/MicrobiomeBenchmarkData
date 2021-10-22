## code to prepare `sysdata` dataset goes here

load("inst/extdata/taxaData.RData")
usethis::use_data(taxa_data, internal = TRUE, overwrite = TRUE)
