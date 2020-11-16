test_that("load_to_deseq2 works", {
  dds <- load_to_deseq2(object, designs_run[["design_1"]])
  expect_true(is(dds, "DESeqDataSet"))
})
