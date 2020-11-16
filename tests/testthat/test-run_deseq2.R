test_that("run_deseq2 works", {
  cat("Running DESEQ2 ...\n")

  for (design_name in names(designs_run)){
    design <- designs_run[[design_name]]
    contrasts <- contrasts_run[[design_name]]
    object_prepro <- preprocess_object(object, min_count=0, min_total_count=15, norm_factors_method="TMM")
    table_deseq2 <- run_deseq2(object_prepro, 
                               design=design,
                               contrasts=contrasts,
                               opts_algo=opts_diffexp$deseq2,
                               opts_comm=opts_diffexp$comm)

    expect_true(unique(table_deseq2$design) == c(design))
    expect_true(all(c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in% colnames(table_deseq2)))
  }
})
