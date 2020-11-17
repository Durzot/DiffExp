test_that("load_to_deseq2 works", {
  dds <- load_to_deseq2(object_example, designs_run[["design_1"]])
  expect_true(is(dds, "DESeqDataSet"))
})
