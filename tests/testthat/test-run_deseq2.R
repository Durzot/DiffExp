test_that("run_deseq2 works", {
  cat("Running DESEQ2 ...\n")

  for (design_name in names(designs_run)){
    design <- designs_run[[design_name]]
    contrasts <- contrasts_run[[design_name]]
    object <- preprocess_object(object_example, opts$prepro)
    table_deseq2 <- run_deseq2(object, 
                               design=design,
                               contrasts=contrasts,
                               opts_algo=opts$diffexp$deseq2,
                               opts_comm=opts$diffexp$comm)

    expect_true(unique(table_deseq2$design) == c(design))
    expect_true(all(c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in% colnames(table_deseq2)))
  }
})
