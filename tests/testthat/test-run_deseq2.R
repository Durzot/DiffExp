test_that("run_deseq2 works", {
  cat("Running DESEQ2 ...\n")

  for (design_name in names(designs_run)){
    df_results <- run_deseq2(object=preprocess_object(object_example, opts$prepro), 
                             design=designs_run[[design_name]],
                             contrasts=contrasts_run[[design_name]],
                             opts_algo=opts$diffexp$deseq2,
                             opts_comm=opts$diffexp$comm)

    expect_true(all(c("beta", "stat", "pval", "padj") %in% colnames(df_results)))
  }
})
