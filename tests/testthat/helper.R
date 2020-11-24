# example data
dds_example <- DESeq2::makeExampleDESeqDataSet(n=5000, m=12, betaSD=1)
dds_example$genotype <- factor(rep(rep(c("I","II"),each=3),2))
object_example <- SummarizedExperiment::SummarizedExperiment(assays=SummarizedExperiment::assays(dds_example),
                                                             colData=SummarizedExperiment::colData(dds_example),
                                                             rowData=SummarizedExperiment::rowData(dds_example),
                                                             metadata=list(name="example data"))

# examples designs and contrasts
designs_run <- list("design_1"= ~ genotype + condition + genotype:condition, 
                    "design_2"= ~ genotype + genotype:condition)

# genotype_I.condition_B_vs_A: the condition effect for genotype I (the main effect)
contrasts_run <- list("design_1"=list("genotype_I.condition_B_vs_A"="conditionB"),
                      "design_2"=list("genotype_I.condition_B_vs_A"="genotypeI.conditionB"))

# load options
opts <- list()
opts[["diffexp"]] <- opts_diffexp(ncores=6, save_table=T, only_significant=T)
opts[["prepro"]] <- opts_prepro(min_count=0, min_total_count=15)

# clean previous run results
unlink(opts$diffexp$comm$folder_results, recursive=T)

# # tcga_brca data
# object_example <- tcga_brca
# 
# # examples designs and contrasts
# designs_run <- list("design_1"= ~ Subtype_Ihc + Age_At_Diagnosis + Age_At_Diagnosis:Subtype_Ihc, 
#                     "design_2"= ~ Age_At_Diagnosis + Age_At_Diagnosis:Subtype_Ihc)
# 
# # Age.IHC_TNBC_vs_HER2_plus: the SubType_Ihc effect difference between TNBC and HER2_plus interacting with Age
# contrasts_run <- list("design_1"=list("Age.IHC_TNBC_vs_HER2_plus"="Subtype_IhcTNBC"),
#                       "design_2"=list("Age.IHC_TNBC_vs_HER2_plus"="Age_At_Diagnosis.Subtype_IhcTNBC"))
# 
# # load options
# opts <- list()
# opts[["diffexp"]] <- opts_diffexp(ncores=6, save_table=T, only_significant=T)
# opts[["prepro"]] <- opts_prepro(min_count=0, min_total_count=15)
# 
# # clean previous run results
# unlink(opts$diffexp$comm$folder_results, recursive=T)
