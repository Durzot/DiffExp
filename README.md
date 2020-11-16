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
    - Skim through the vignette
      <https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>
    2. **TODO**
    - Understand the usage of cook's distances, the filter threshold selection for independent filtering and
      the refitting at the end of `estimateDispersionGeneEst` and some other refittings.

- `edgeR`
    1. **Done**
    - Read the short paper <https://doi.org/10.1093/bioinformatics/btp616>
    - Read the extended paper <https://doi.org/10.1093/nar/gks042>
    - Skim through the user guide.

- `limmaVoom`
    1. **Done**
    - Read the paper <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29>
    - Read the paper <https://academic.oup.com/nar/article/43/7/e47/2414268>
    2. **TODO**
    - Read the code in details
    - Used it on a project.
    - Skimmed through the vignette
      <https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html>


## The `DiffExp` package.

The package is under development. Details about functions may be found in the manual package in the root folder.
