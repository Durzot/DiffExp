# @modified: 12 Nov 2020
# @created: 28 Oct 2020
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
# 
#     Institut Gustave Roussy
#     Prism Center
#     114 rue Edouard Vaillant, Villejuif, 94800 France
# 
# Apply the DESeq 2 R library on the RNA-seq data.
# 
# Reference
# ---------
# Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with
# DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

#' Define the parameters specific to each differential analysis method.
#'
#' @param alpha fdr level when adjusting for multiple testing.
#' @param minreads the minimum number of reads a variable must have to be differentially analyzed.
#' @param ncores number of cores available for doing parallel computations. Used in DESeq.
#' @param run_deseq2 boolean to choose to run DESeq2.
#' @param run_edgeR boolean to choose to run edgeR.
#' @param run_limma boolean to choose to run limma.
#' @param ... extra parameters added to the configuration list
#'
#' @importFrom utils modifyList
#'
#' @author Yoann Pradat
#'
#' @export
config_default <- function(alpha=0.1,
                           minreads=10,
                           ncores=4,
                           run_deseq2=T,
                           run_edgeR=T,
                           run_limma=T, ...){

  common <- list(alpha=alpha,
                 minreads=minreads,
                 ncores=ncores)
  deseq2 <- list(run=run_deseq2, lfcShrink_type="apeglm")
  edgeR <- list(run=run_edgeR)
  limma <- list(run=run_limma)

  default <- list(common=common,
                  deseq2=deseq2,
                  edgeR=edgeR,
                  limma=limma)

  modifyList(default, list(...))
}
