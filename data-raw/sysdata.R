## code to prepare `sysdata` dataset goes here

subgingival_supragingival_hmp16s_intersects <- MSEBenchmarkData:::.hmp16SDataIntersects( c("Subgingival Plaque", "Supragingival Plaque"))
stool_tongue_dorsum_hmp16s_intersects <- MSEBenchmarkData:::.hmp16SDataIntersects( c("Stool", "Tongue Dorsum"))

usethis::use_data(subgingival_supragingival_hmp16s_intersects,
                  stool_tongue_dorsum_hmp16s_intersects,
                  internal = TRUE, overwrite = TRUE)
