# @modified: 16 Nov 2020
# @created: 06 Nov 2020
# @author: Yoann Pradat
#
# One function for loading the data to the object expected in each differential analysis procedure.

#' Transform \code{SummarizedExperiment} object to \code{DESeqDataSet} object.
#'
#' @return a \code{DESeqDataSet} object
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula or matrix that expresses how the counts for each variable depend on the variables in
#' \code{colData}. See \code{DESeqDataSet}.
#'
#' @import SummarizedExperiment
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom S4Vectors metadata
#'
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_deseq2 <- function(object, design=NULL){
  # counts
  count_data <- as.matrix(assays(object)$counts)
  mode(count_data) <- "integer"

  # make the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData=count_data,
                                colData=colData(object),
                                design=design,
                                metadata=metadata(object))

  # add row data
  mcols(dds) <- cbind(mcols(dds), rowData(object))

  return(dds)
}

#' Transform \code{SummarizedExperiment} object to \code{DGEList} object.
#'
#' @param object a \code{SummarizedExperiment} object
#'
#' @importFrom edgeR cpm
#' @author Yoann Pradat
#'
#' @keywords internal
load_to_edgeR <- function(object){
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
