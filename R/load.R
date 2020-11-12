# @modified: 12 Nov 2020
# @created: 06 Nov 2020
# @author: Yoann Pradat
#
# One function for loading the data to the object expected in each differential analysis procedure.

#' Transform \code{SummarizedExperiment} object to \code{\link[DESeq2]{DESeq2DataSet}} object.
#'
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula or matrix that expresses how the counts for each variable depend on the variables in
#' \code{colData}. See \code{\link[DESeq2]{DESeqDataSet}}.
#'
#' @import SummarizedExperiment
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @author Yoann Pradat
#'
#' @export
load_to_deseq2 <- function(object, design=NULL){
  # counts
  count_data <- as.matrix(assays(object)$counts)
  mode(count_data) <- "integer"

  # colData
  col_data <- colData(object)

  # make the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design =  design)

  # add row data
  mcols(dds) <- cbind(mcols(dds), rowData(object))

  return(dds)
}

#' Transform \code{SummarizedExperiment} object to \code{\link[edgeR]{DGEList}} object.
#'
#' @param object a \code{SummarizedExperiment} object
#'
#' @importFrom edgeR cpm
#' @author Yoann Pradat
#'
#' @export
load_to_edgeR <- function(object){
  return(T)
}

#' Transform \code{SummarizedExperiment} object to ????
#'
#' @param object a \code{SummarizedExperiment} object
#'
#' @author Yoann Pradat
#'
#' @export
load_to_limma <- function(object){
  return(T)
}
