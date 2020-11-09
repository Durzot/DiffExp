# Differential Expression

Sets of methods used to perform differential expression analysis between medical cohorts. The code for the original
methods is included as git submodules in the `tools` folder of this repository. Annotations and/or presentations about
the methods may be found in the `docs` folder.

For now, the following methods have been studied:

- `DESeq2`
    1. **Done**
    - Read the paper <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8>
    - Read the code in details
    - Used it on a project.
    - Skimmed through the vignette
      <https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>
    2. **TODO**
    - Understand the usage of cook's distances, the filter threshold selection for independent filtering and
      the refitting at the end of `estimateDispersionGeneEst` and some other refittings.

- `limmaVoom`
    1. **TODO**
    - Read the paper <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29>
    - Read the code in details
    - Used it on a project.
    - Skimmed through the vignette
      <https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html>
