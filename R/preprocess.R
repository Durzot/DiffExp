# @modified: 16 Nov 2020
# @created: 16 Nov 2020
# @author: Yoann Pradat
#
# Perform preprocessing options like computing normalization factors and removing low-count variables

#' Perform preprocessing steps.
#'
#' @return a \code{SummarizedExperiment} object
#' @param object a \code{SummarizedExperiment} object
#' @param opts a named list of options. See \code{\link{opts_prepro}}
#'  
#' @importFrom edgeR filterByExpr calcNormFactors
#' @importFrom stats model.matrix
#'
#' @author Yoann Pradat
#' @export
preprocess_object <- function(object, opts){
  # options
  design <- opts$design
  min_count <- opts$min_count
  min_total_count <- opts$min_total_count
  large_n <- opts$large_n
  min_prop <- opts$min_prop
  norm_factors_method <- opts$norm_factors_method

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
  if ("norm.factors" %in% colnames(rowData(object))){
    rowData(object)[,"norm.factors"] <- NULL
  }
  object$norm.factors <- calcNormFactors(counts, method=norm_factors_method)

  # update counts
  assays(object)$counts <- counts

  object
}
