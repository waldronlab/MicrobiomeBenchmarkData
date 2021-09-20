test_that("getDatasets works with datasets indexes.", {

  index_vector_ok <- 1:3
  index_vector_duplicates <- c(1, 2, 2, 3, 4)
  index_vector_dbl_duplicates <- c(1.1, 1, 2, 3, 4)
  index_vector_na <- c(1, 2, NA, 3, NA)
  index_vector_invalid_value <- c(1, 2, 3, "calgaro_2020_16S_gingival_healthy")

  expect_error(getDatasets(index = index_vector_duplicates), regexp = "duplicated")
  expect_error(getDatasets(index = index_vector_dbl_duplicates), regexp = "duplicated")
  expect_error(getDatasets(index = index_vector_na),  regexp = "Invalid")
  expect_error(getDatasets(index = index_vector_invalid_value), regexp = "Invalid")

  index_output <- getDatasets(index = index_vector_ok)
  expect_s4_class(index_output[[1]], "TreeSummarizedExperiment")

})

test_that("getDatasets works with datasets names.", {

  name_vector_ok <- c("calgaro_2020_16S_gingival_healthy", "HMP_WMS_cMD3_gingival_healthy")
  name_vector_not_chr <- 1:3
  name_vector_duplicated <- c("calgaro_2020_16S_gingival_healthy", "calgaro_2020_16S_gingival_healthy", "calgaro_2020_16S_gingival_healthy")
  name_vector_na <- c("calgaro_2020_16S_gingival_healthy", NA, NA)
  name_vector_not <- c("calgaro_2020_16S_gingival_healthy", "calgaro_2020", "calgaro_2019")

  expect_error(getDatasets(name = name_vector_not_chr), regexp = "character")
  expect_error(getDatasets(name = name_vector_duplicated), regexp = "duplicated")
  expect_error(getDatasets(name = name_vector_na), regexp = "NA")
  expect_error(getDatasets(name = name_vector_not), regexp = "list")

  name_output <- getDatasets(name = name_vector_ok)
  expect_s4_class(name_output[["calgaro_2020_16S_gingival_healthy"]], "TreeSummarizedExperiment")

})
