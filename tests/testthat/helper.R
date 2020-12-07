# simulation example
synthetic_example <- function(n=5000, m=12, betaSD=1, seed=25){
  set.seed(seed)
  dds_example <- DESeq2::makeExampleDESeqDataSet(n=n, m=m, betaSD=betaSD)
  dds_example$genotype <- factor(rep(rep(c("I","II"),each=3),2))

  object <- SummarizedExperiment::SummarizedExperiment(assays=SummarizedExperiment::assays(dds_example),
                                                       colData=SummarizedExperiment::colData(dds_example),
                                                       rowData=SummarizedExperiment::rowData(dds_example),
                                                       metadata=list(name="example data"))
  # examples designs and contrasts
  designs_run <- list("design_1"= ~ genotype + condition + genotype:condition, 
                      "design_2"= ~ genotype + genotype:condition)

  # genotype_I.condition_B_vs_A: the condition effect for genotype I (the main effect)
  contrasts_run <- list("design_1"=list("genotype_I.condition_B_vs_A"="conditionB"),
                        "design_2"=list("genotype_I.condition_B_vs_A"="genotypeI.conditionB"))

  # options
  opts <- list()
  opts[["diffexp"]] <- opts_diffexp(alpha=0.1, ncores=6, save_table=T, only_significant=T, 
                                    folder_results="./out/results_synthetic")
  opts[["prepro"]] <- opts_prepro(min_count=0, min_total_count=15)

  list(object=object, designs_run=designs_run, contrasts_run=contrasts_run, opts=opts)
}

# real example
real_example <- function(){
  # examples designs and contrasts
  designs_run <- list("design_1"= ~ Subtype_Ihc + Age_At_Diagnosis + Age_At_Diagnosis:Subtype_Ihc, 
                      "design_2"= ~ Age_At_Diagnosis + Age_At_Diagnosis:Subtype_Ihc)

  # Age.IHC_TNBC_vs_HER2_plus: the SubType_Ihc effect difference between TNBC and HER2_plus interacting with Age
  contrasts_run <- list("design_1"=list("Age.IHC_TNBC_vs_HER2_plus"="Age_At_Diagnosis"),
                        "design_2"=list("Age.IHC_TNBC_vs_HER2_plus"="Age_At_Diagnosis.Subtype_IhcTNBC"))

  # load options
  opts <- list()
  opts[["diffexp"]] <- opts_diffexp(ncores=6, save_table=T, only_significant=T,
                                    folder_results="./out/results_real")
  opts[["prepro"]] <- opts_prepro(min_count=0, min_total_count=15)

  list(object=tcga_brca, designs_run=designs_run, contrasts_run=contrasts_run, opts=opts)
}

sexample <- synthetic_example()
rexample <- real_example()
