# example data
dds_example <- DESeq2::makeExampleDESeqDataSet(n=100, m=12)
dds_example$genotype <- factor(rep(rep(c("I","II"),each=3),2))
object <- SummarizedExperiment::SummarizedExperiment(assays=SummarizedExperiment::assays(dds_example),
                                                     colData=SummarizedExperiment::colData(dds_example),
                                                     rowData=SummarizedExperiment::rowData(dds_example),
                                                     metadata=list(name="example data"))
design_a <- ~ genotype + genotype:condition
design_b <- ~ genotype + condition + genotype:condition
