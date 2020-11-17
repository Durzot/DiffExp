# @modified: 16 Nov 2020
# @created: 06 Nov 2020
# @author: Yoann Pradat
#
# One function for loading the data to the object expected in each differential analysis procedure.

#' Transform \code{SummarizedExperiment} object to \code{DESeqDataSet} object.
#'
#' @return a \code{DESeqDataSet} object
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula that expresses how the counts for each variable depend on the variables in
#' \code{colData}
#'
#' @import SummarizedExperiment
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom S4Vectors metadata
#'
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_deseq2 <- function(object, design=NULL){
  # make the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData=assays(object)$counts,
                                colData=colData(object),
                                design=design,
                                metadata=metadata(object))
  # add row data
  mcols(dds) <- cbind(mcols(dds), rowData(object))

  if ("norm.factors" %in% colnames(colData(dds))){
    sizeFactors(dds) <- colData(dds)$norm.factors
    colData(dds)$norm.factors <- NULL
  }

  dds
}

#' Transform \code{SummarizedExperiment} object to \code{DGEList} object.
#'
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula that expresses how the counts for each variable depend on the variables in
#' \code{colData}
#'
#' @import SummarizedExperiment
#' @importFrom edgeR cpm
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_edgeR <- function(object){

  # norm factors if any
  if ("norm.factors" %in% colnames(colData(object))){
    norm_factors <- colData(object)$norm.factors
  } else{
    norm_factors <- NULL
  }

  # DGEList
  dgel <- DGEList(counts=assays(object)$count,
                  lib.size=NULL, # default to colSums(counts)
                  norm.factors=norm_factors,
                  samples=colData(object),
                  gene=rowData(object),
                  remove.zeros=F)

  dgel
}

run_edgeR

#' Transform \code{SummarizedExperiment} object to ????
#'
#' @param object a \code{SummarizedExperiment} object
#'
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_limma <- function(object){
  T
}
