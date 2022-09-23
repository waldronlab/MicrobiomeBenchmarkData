test_that("getBenchmarkData works.", {
    dataset_names <- getBenchmarkData()
    # expect_true(is.character(dataset_names))
    expect_true(is.data.frame(dataset_names))
    expect_message(getBenchmarkData(), "Use vignette")

    tse <- getBenchmarkData("HMP_2012_16S_gingival_V35_subset", dryrun = FALSE)[[1]]
    expect_s4_class(tse, "TreeSummarizedExperiment")
})
