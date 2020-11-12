test_that("load_to_deseq2 works", {
  dds <- load_to_deseq2(object, design_a)
  expect_true(is(dds, "DESeqDataSet"))
})
