test_that("load_to_deseq2 works", {
  dds <- load_to_deseq2(sexample$object, design=~1)
  expect_true(is(dds, "DESeqDataSet"))

  # rounding is required
  object <- rexample$object
  assays(object)$counts <- round(assays(object)$counts)
  dds <- load_to_deseq2(object, design=~1)
  expect_true(is(dds, "DESeqDataSet"))
})

test_that("load_to_edgeR works", {
  dge <- load_to_edgeR(sexample$object)
  expect_true(is(dge, "DGEList"))

  dge <- load_to_edgeR(rexample$object)
  expect_true(is(dge, "DGEList"))
})
