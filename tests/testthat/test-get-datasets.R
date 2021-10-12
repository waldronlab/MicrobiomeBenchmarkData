test_that("getDataset works.", {

  dataset_names <- getDataset()
  expect_true(is.character(dataset_names))
  expect_message(getDataset(), "Use vignette")

  tse <- getDataset("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]
  expect_s4_class(tse, "TreeSummarizedExperiment")

})
