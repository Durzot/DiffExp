# example data
dds_example <- DESeq2::makeExampleDESeqDataSet(n=1000, m=12, betaSD=1)
dds_example$genotype <- factor(rep(rep(c("I","II"),each=3),2))
object_example <- SummarizedExperiment::SummarizedExperiment(assays=SummarizedExperiment::assays(dds_example),
                                                             colData=SummarizedExperiment::colData(dds_example),
                                                             rowData=SummarizedExperiment::rowData(dds_example),
                                                             metadata=list(name="example data"))

# examples designs and contrasts
designs_run <- list("design_1"= ~ genotype + condition + genotype:condition, 
                    "design_2"= ~ genotype + genotype:condition)

# genotype_I.condition_B_vs_A: the condition effect for genotype I (the main effect)
contrasts_run <- list("design_1"=list("genotype_I.condition_B_vs_A"="condition_B_vs_A"),
                      "design_2"=list("genotype_I.condition_B_vs_A"="genotypeI.condition_B_vs_A"))

# load options
opts <- list()
opts[["diffexp"]] <- opts_diffexp(ncores=6, save_table=T, only_significant=T)
opts[["prepro"]] <- opts_prepro(min_count=0, min_total_count=15)

# clean previous run results
unlink(opts$diffexp$comm$folder_results, recursive=T)
