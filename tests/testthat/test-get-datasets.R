test_that("getDatasets works.", {

  datasets <- getDatasets()

  vector_ok <- names(getDatasets(datasets[1:2]))

  vector_not_chr <- 1:3
  vector_duplicated <- c("calgaro_2020_16S_gingival_healthy", "calgaro_2020_16S_gingival_healthy", "calgaro_2020_16S_gingival_healthy")
  vector_na <- c("calgaro_2020_16S_gingival_healthy", NA, NA)
  vector_not <- c("calgaro_2020_16S_gingival_healthy", "calgaro_2020", "calgaro_2019")

  output <- getDatasets(vector_ok)
  expect_s4_class(output[[1]], "SummarizedExperiment")

  expect_error(getDatasets(vector_not_chr), regexp = "character")
  expect_error(getDatasets(vector_duplicated), regexp = "duplicate")
  expect_error(getDatasets(vector_na), regexp = "NA")

})
