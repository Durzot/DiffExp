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

#' Define the options for \code{\link{preprocess_object}}
#'
#' @return a \code{SummarizedExperiment} object
#' @param design a formula object specifying the design matrix
#' @param min_count \code{min.total.count} param of \code{\link[edgeR]{filterByExpr}}
#' @param min_total_count \code{min.count} param of \code{\link[edgeR]{filterByExpr}}
#' @param large_n \code{large.n} param of \code{\link[edgeR]{filterByExpr}}
#' @param min_prop \code{min.prop} param of \code{\link[edgeR]{filterByExpr}}
#' @param norm_factors_method \code{method} argument of \code{\link[edgeR]{calcNormFactors}}.
#' @param ... extra parameters added to the options list
#'
#' @importFrom utils modifyList
#'
#' @author Yoann Pradat
#'
#' @export
opts_prepro <- function(design=NULL,
                        min_count=0,
                        min_total_count=15,
                        large_n=10,
                        min_prop=0.7,
                        norm_factors_method=c("TMM", "TMMwsp", "RLE", "upperquartile", "none"), ...){

  prepro <- list(design=design,
                 min_count=min_count,
                 min_total_count=min_total_count,
                 norm_factors_method=norm_factors_method,
                 large_n=large_n,
                 min_prop=min_prop,
                 norm_factors_method=norm_factors_method)

  modifyList(prepro, list(...))
}


#' Define the parameters specific to each differential analysis method.
#'
#' @param alpha fdr level when adjusting for multiple testing.
#' @param ncores number of cores available for doing parallel computations. Used in DESeq.
#' @param save_table boolean to decide whether to save tables in txt files or not
#' @param only_significant boolean to decide whether only significant (FDR) variables are kept in the results tables
#' or not
#' @param use_deseq2 boolean to choose to run DESeq2.
#' @param use_edgeR boolean to choose to run edgeR.
#' @param use_limma boolean to choose to run limma.
#' @param ... extra parameters added to the configuration list
#'
#' @importFrom utils modifyList
#'
#' @author Yoann Pradat
#'
#' @export
opts_diffexp <- function(alpha=0.1,
                         ncores=6,
                         save_table=T,
                         only_significant=T,
                         folder_results="./results",
                         use_deseq2=T,
                         use_edgeR=F,
                         use_limma=F, ...){

  common <- list(alpha=alpha, 
                 ncores=ncores,
                 save_table=save_table,
                 only_significant=only_significant,
                 folder_results=folder_results)

  deseq2 <- list(run=use_deseq2, 
                 lfc_fit_shrink_type="apeglm",
                 lfc_test_alt_hypothesis="greaterAbs",
                 lfc_test_threshold=0)

  edgeR <- list(run=use_edgeR,
                disp_trend_method="locfit",
                disp_robust_prior_df=TRUE,
                lfc_fit_robust_prior_ql_disp=TRUE,
                lfc_test_type=c("QL","LRT"),
                lfc_test_threshold=0)

  limma <- list(run=use_limma)

  default <- list(common=common,
                  deseq2=deseq2,
                  edgeR=edgeR,
                  limma=limma)

  modifyList(default, list(...))
}
