# @modified: 16 Nov 2020
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

#' Define the parameters specific to each differential analysis method.
#'
#' @param alpha fdr level when adjusting for multiple testing.
#' @param ncores number of cores available for doing parallel computations. Used in DESeq.
#' @param save_table boolean to decide whether to save tables in txt files or not
#' @param only_significant boolean to decide whether only significant (FDR) variables are kept in the results tables
#' or not
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
opts_diffexp_default <- function(alpha=0.1,
                                 ncores=6,
                                 save_table=T,
                                 only_significant=T,
                                 folder_results="./results",
                                 run_deseq2=T,
                                 run_edgeR=F,
                                 run_limma=F, ...){

  common <- list(alpha=alpha, 
                 ncores=ncores,
                 save_table=save_table,
                 only_significant=only_significant,
                 folder_results=folder_results)
  deseq2 <- list(run=run_deseq2, 
                 lfcShrink_type="apeglm",
                 altHypothesis="greaterAbs",
                 lfcThreshold=0)
  edgeR <- list(run=run_edgeR)
  limma <- list(run=run_limma)

  default <- list(common=common,
                  deseq2=deseq2,
                  edgeR=edgeR,
                  limma=limma)

  modifyList(default, list(...))
}
