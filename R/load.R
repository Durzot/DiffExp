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

  if ("normFactors" %in% colnames(colData(dds))){
    sizeFactors(dds) <- colData(dds)$normFactors
    colData(dds)$normFactors <- NULL
  }

  return(dds)
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
load_to_edgeR <- function(object, design){
  # counts
  count_data <- as.matrix(assays(object)$counts)
  mode(count_data) <- "integer"
  
  # design
  design_matrix <- model.matrix(design, data=colData(object))

  # DGEList
  dgel <- DGEList(counts=count_data,
                  lib.siz=colSums(count_data),
                  norm.factors=NULL,
                  gene=rowData(object))
  # method=RLE is the same scaling method as DESeq2
  # method=TMM is trimmed mean M-values method.
  dgel <- calcNormFactors(object=dgel, method="RLE") 
  dgel <- filterByExpr(dgel, design_matrix)

  return(T)
}

#' Transform \code{SummarizedExperiment} object to ????
#'
#' @param object a \code{SummarizedExperiment} object
#'
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_limma <- function(object){
  return(T)
}
