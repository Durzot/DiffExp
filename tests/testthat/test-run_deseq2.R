test_that("run_deseq2 synthetic works", {
  object <- sexample$object
  designs_run <- sexample$designs_run
  contrasts_run <- sexample$contrasts_run
  opts <- sexample$opts
  for (design_name in names(designs_run)){
    df_results <- run_deseq2(object=preprocess_object(object, opts$prepro), 
                             design=designs_run[[design_name]],
                             contrasts=contrasts_run[[design_name]],
                             opts_algo=opts$diffexp$deseq2,
                             opts_comm=opts$diffexp$comm)

    expect_true(all(c("beta", "stat", "pval", "padj") %in% colnames(df_results)))
  }
})

test_that("run_deseq2 real works", {
  object <- rexample$object
  designs_run <- rexample$designs_run
  contrasts_run <- rexample$contrasts_run
  opts <- rexample$opts
  for (design_name in names(designs_run)){
    df_results <- run_deseq2(object=preprocess_object(object, opts$prepro), 
                             design=designs_run[[design_name]],
                             contrasts=contrasts_run[[design_name]],
                             opts_algo=opts$diffexp$deseq2,
                             opts_comm=opts$diffexp$comm)

    expect_true(all(c("beta", "stat", "pval", "padj") %in% colnames(df_results)))
  }
})
