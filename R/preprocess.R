# @modified: 16 Nov 2020
# @created: 16 Nov 2020
# @author: Yoann Pradat
#
# Perform preprocessing options like computing normalization factors and removing low-count variables

#' Perform preprocessing steps.
#'
#' @return a \code{SummarizedExperiment} object
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula object specifying the design matrix
#' @param min_count \code{min.count} param of \code{\link[edgeR]{filterByExpr}}
#' @param min_total_count \code{min.count} param of \code{\link[edgeR]{filterByExpr}}
#' @param large_n \code{large.n} param of \code{\link[edgeR]{filterByExpr}}
#' @param min_prop \code{min.prop} param of \code{\link[edgeR]{filterByExpr}}
#' @param norm_factors_method \code{method} argument of \code{\link[edgeR]{calcNormFactors}}.
#'  
#' @importFrom edgeR filterByExpr calcNormFactors
#'
#' @author Yoann Pradat
#' @export
preprocess_object <- function(object, design=NULL, min_count=10, min_total_count=15, large_n=10, min_prop=0.7,
                              norm_factors_method=c("TMM", "TMMwsp", "RLE", "upperquartile", "none")){

  # get design matrix
  if (is.null(design)){
    design_matrix <- NULL
  } else{
    design_matrix <- model.matrix(design, data=colData(object))
  }

  # check integer counts
  counts <- as.matrix(assays(object)$counts)
  if (!all(counts==floor(counts))){
    message("some counts are not integer, transforming them to integers")
  }
  mode(counts) <- "integer"

  # remove low-counts variables
  keep <- filterByExpr(y=counts, design=design_matrix,
                       min.count=min_count, min.total.count=min_total_count,
                       large.n=large_n, min.prop=min_prop)
  message(paste0("edgeR::filterByExpr keeps ", sum(keep), "/", length(keep), " variables")) 
  object <- object[keep,]
  counts <- counts[keep,]

  # compute normalization factors
  if ("normFactors" %in% colnames(rowData(object))){
    rowData(object)[,"normFactors"] <- NULL
  }
  object$normFactors <- calcNormFactors(counts, method=norm_factors_method)

  # update counts
  assays(object)$counts <- counts

  object
}
