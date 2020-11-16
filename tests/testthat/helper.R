# example data
dds_example <- DESeq2::makeExampleDESeqDataSet(n=1000, m=12, betaSD=1)
dds_example$genotype <- factor(rep(rep(c("I","II"),each=3),2))
object <- SummarizedExperiment::SummarizedExperiment(assays=SummarizedExperiment::assays(dds_example),
                                                     colData=SummarizedExperiment::colData(dds_example),
                                                     rowData=SummarizedExperiment::rowData(dds_example),
                                                     metadata=list(name="example data"))

# examples designs and contrasts
designs_run <- list("design_1"= ~ genotype + condition + genotype:condition, 
                    "design_2"= ~ genotype + genotype:condition)

# contrast_1: the condition effect for genotype I (the main effect)
contrasts_run <- list("design_1"=list("contrast_1"=list("condition_B_vs_A")),
                      "design_2"=list("contrast_1"=list("genotypeI.conditionB")))

# load default options
opts_diffexp <- opts_diffexp_default(ncores=6, save_table=T, only_significant=T)

# clean previous run results
unlink(opts_diffexp$comm$folder_results, recursive=T)
