test_that("run edgeR works", {
  cat("Running edgeR ...\n")

  for (design_name in names(designs_run)){
    df_results <- run_edgeR(object=preprocess_object(object_example, opts$prepro), 
                            design=designs_run[[design_name]],
                            contrasts=contrasts_run[[design_name]],
                            opts_algo=opts$diffexp$edgeR,
                            opts_comm=opts$diffexp$comm)

    expect_true(all(c("beta", "stat", "pval", "padj") %in% colnames(df_results)))
  }
})
