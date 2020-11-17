# @modified: 16 Nov 2020
# @created: 12 Nov 2020
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
# One function for running the main functions of each differential analysis procedure.
 

#' Run DESEQ2 algorithm.
#' 
#' @return dataframe with results
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula specifying the design for the model matrix of DESeq2. Any variable appearing should be present
#' in the \code{colData} of object
#' @param contrasts (optional) a character vector specifying the contrasts (one or multiple beta coefficient) to be used
#' for making tests and building results table. If set to NULL, one table for each element of \code{resultsNames(dds)}
#' will be built
#' @param opts_algo a named list of options passed to \code{DESeq}, \code{results} and \code{lfcShrink} functions. See
#' \code{\link[DESeq2]{DEseq}}, \code{\link[DESeq2]{results}}, \code{\link[DESeq2]{lfcShrink}}
#'
#' @import DESeq2
#' @importFrom BiocParallel MulticoreParam
#'
#' @author Yoann Pradat
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with
#' DESeq2. Genome Biology, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
run_deseq2 <- function(object, design=NULL, contrasts, opts_algo, opts_comm){

  # 0. options and settings ============================================================================================

  # Options for the function \code{DESeq}
  lfcThreshold <- opts_algo$lfcThreshold
  ncores <- opts_comm$ncores

  # Options for the function \code{results}
  alpha <- opts_comm$alpha

  # Options for the function \code{lfcShrink}
  lfcShrink_type <- opts_algo$lfcShrink_type
  altHypothesis <- opts_algo$altHypothesis 

  # Options for saving
  save_table <- opts_comm$save_table
  only_significant <- opts_comm$only_significant

  # init results folder
  folder_results <- file.path(opts_comm$folder_results, "deseq2")
  dir.create(folder_results, showWarnings=F, recursive=T)

  # used only if save_table is TRUE
  file_table_results <- file.path(folder_results, "table_results.txt")

  # for file names and plot titles
  meta_char <- paste(lapply(metadata(object), function(x) paste(x, collapse="_")), collapse="-")
  design_char <- paste(design, collapse="")

  # 1. build dds object ================================================================================================
  cat("Building the DESeqDataSet object ...\n")
  dds <- load_to_deseq2(object=object, design=design)

  # 2. fit the deseq2 models ===========================================================================================
  cat("Fitting the DESeq2 models ...\n")
  BPPARAM <- MulticoreParam(ncores)

  # fit DESeq2 model in parallel without betaPrior
  #
  # Perform the following
  # - Estimate size factors
  # - Estimate dispersions
  #   * estimate gene-wise dispersions
  #   * estimate dispersions prior
  #   * estimate dispersions MAP
  # - nBinomWaldTest with betaPrior=False, test="Wald"
  #   * fitNbinomGLMs with a wide prior (\sigma^2_r = 10^6) using IRLS to get initial estimates of the \beta_{i,r}
  dds <- DESeq(object=dds, 
               test="Wald", 
               fitType="parametric", 
               sfType="ratio", 
               quiet=F, 
               parallel=T, 
               BPPARAM=BPPARAM)

  # save plot dispersion
  filename <- paste0("deseq2_dispersion_", meta_char, "_", design_char, ".pdf")
  title <- paste("Dispersion plot on:", meta_char)
  plot_dispersion_deseq2(dds, filepath=file.path(folder_results, filename), title=title, subtitle=design_char)

  # 3. results w/wo lcf shrinkage ===========================================================================================
  
  # # Exploration of the link between designs
  # # - ~genotype + condition + genotype:condition
  # # - ~genotype + genotype:condition
  #
  # dds1 <- makeExampleDESeqDataSet(n=100,m=18)
  # dds1$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
  # design(dds1) <- ~ genotype + condition + genotype:condition
  # dds1 <- DESeq(dds1)
  # design(dds2) <- ~ genotype + genotype:condition
  # dds2 <- DESeq(dds2)
  # resultsNames(dds1)
  # resultsNames(dds2)
  # 
  # # the condition effect for genotype I (the main effect)
  # results(dds1, contrast=c("condition","B","A"))
  # results(dds2, contrast=c(list("genotypeI.conditionB")))
  # 
  # # the condition effect for genotype III.
  # # this is the main effect *plus* the interaction term
  # # (the extra condition effect in genotype III compared to genotype I).
  # results(dds1, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
  # results(dds2, contrast=c(list("genotypeIII.conditionB")))
  #  
  # # the interaction term for condition effect in genotype III vs genotype I.
  # # this tests if the condition effect is different in III compared to I
  # results(dds1, name="genotypeIII.conditionB")
  # results(dds2, contrast=list("genotypeIII.conditionB", "genotypeI.conditionB"))
  # 
  # # the interaction term for condition effect in genotype III vs genotype II.
  # # this tests if the condition effect is different in III compared to II
  # results(dds1, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))
  # results(dds2, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))

  # Extract table results with base means, log2 fold changes, standard errors, test statistics,
  # p-values, adjusted p-values for different priors over LFC.

  if (is.null(contrasts)){
    contrasts <- resultsNames(dds)
    contrasts <- contrasts[contrasts != "Intercept"]
    contrasts <- sapply(contrasts, list)
  }

  cat("Computing results tables for each for the following contrasts...\n\t")
  cat(paste(contrasts, collapse="\n\t"), "\n")

  table_results_all <- data.frame()

  for (contrast_name in names(contrasts)){
    # No LFC shrinkage
    # The statistic is based on the MLE of LFC
    resLFC_MLE <- results(object=dds,
                          contrast=contrasts[[contrast_name]],
                          lfcThreshold=lfcThreshold,         # if specifying lfcThreshold, tests are Wald tests
                          altHypothesis=altHypothesis,
                          alpha=alpha,                       # for p-value adjustment level
                          filter=mcols(dds)$baseMean,
                          independentFiltering=T,            # for choosing genes whose p-values will be ajdusted
                          pAdjustMethod="BH",
                          format="DataFrame",
                          parallel=T,
                          BPPARAM=BPPARAM)

    lfc_description <- attr(resLFC_MLE, "elementMetadata")[2,"description"]

    # WARNING: the p-values are unchanged, only lfc and lfcSE change (with lfcSE being a posterior SD here).
    resLFC_shrink <- NULL

    if (is.null(lfcShrink_type)){
      resLFC_shrink <- NULL
    } else if (lfcShrink_type=="normal"){
      if (any(attr(terms.formula(design(dds)),"order") > 1)){
        cat("Interactions with lfcShrink='normal' is not implemented. Skipping.\n")
      } else {
      # Normal LFC prior is a zero-centered Gaussian (original in DESeq2). Not available with interactions.
      resLFC_shrink <- lfcShrink(dds=dds, 
                                 contrast=contrasts[[contrast_name]],
                                 type=lfcShrink_type,
                                 format="DataFrame",
                                 parallel=T,
                                 BPPARAM=BPPARAM)
      }
    } else if(lfcShrink_type=="apeglm"){
      if (length(contrasts[[contrast_name]])>1){
          cat("Contrasts that are not coefficient names for lfcShrink='apeglm' is not implemented. Skipping.\n")
      } else {
        # Apeglm sets a Cauchy prior with null location
        resLFC_shrink <- lfcShrink(dds=dds, 
                                   coef=unlist(contrasts[[contrast_name]]),
                                   type=lfcShrink_type,
                                   format="DataFrame",
                                   parallel=T,
                                   BPPARAM=BPPARAM)
      }
    }
  
    plot_MA_deseq2(resLFC_MLE=resLFC_MLE,
                   resLFC_shrink=resLFC_shrink,
                   filepath=file.path(folder_results, paste0("deseq2_ma_", meta_char, "_", design_char, ".pdf")), 
                   title=paste("LFC over mean norm counts on:", meta_char),
                   subtitle=paste(design_char, lfc_description))

    # summary of results
    #
    # if "svalue" not in names of the results object (when is it the case?), outliers and low counts
    #   outlier <- sum(object$baseMean > 0 & is.na(object$pvalue)) (value after "outliers [1] :")
    #   filt <- sum(!is.na(object$pvalue) & is.na(object$padj))    (value after "low counts [2] :")
    #
    # and the value after "(mean count <") is
    #
    #   if (is.null(metadata(object)$filterThreshold)) {
    #     ft <- 0
    #   } else {
    #     ft <- round(metadata(object)$filterThreshold)
    #   }

    # default filter value for independent filtering is the baseMean (mean of normalized counts) of the dds object
    # the number of low counts variable should be
    #   sum(mcols(dds)$baseMean <= round(metadata(resLFC_shrink)[["filterThreshold"]]))

    # save MA plot
    # "In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of 
    # normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is
    # less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down."

    # table of results
    if (only_significant){
      keep <- !is.na(resLFC_MLE$padj) & resLFC_MLE$padj < alpha
    } else {
      keep <- rep(T, nrow(resLFC_MLE))
    }

    if (grepl("contrast_ref", contrast_name)){
      contrast_char <- sub("contrast_ref_", "\\2", contrast_name)
    } else {
      contrast_char <- sub("(.*)(?=:\\s):\\s", "\\2", lfc_description, perl=T)
    }

    tags <- list("meta"=meta_char, "design"=design_char, "contrast"=contrast_char)
    extra <- list(time=as.character(Sys.time()))
    table_results <- cbind.data.frame(c(tags, extra), resLFC_MLE[keep,])
    table_results <- cbind.data.frame(rowData(object[keep,]), table_results)
    
    if (save_table){
      save_update_table(file_table_results, table_results, tags)
    }

    table_results_all <- rbind(table_results_all, table_results)
  }

  table_results_all
}

#' Run edgeR algorithm.
#' 
#' @author Yoann Pradat
#'
#' @references
#' Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of
#' digital gene expression data.” Bioinformatics, 26(1), 139-140. \url{https://doi.org/10.1093/bioinformatics/btp616}
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with
#' respect to biological variation.” Nucleic Acids Research, 40(10), 4288-4297. \url{https://doi.org/10.1093/nar/gks042}
run_edgeR <- function(object, design=NULL, opts_algo, opts_comm){
  # 0. options and settings ============================================================================================

  # Options for the function \code{estimateDisp}
  disp_trend_method <- opts_algo$disp_trend_method
  disp_robust <- opts_algo$disp_robust

  # Options for saving
  save_table <- opts_comm$save_table
  only_significant <- opts_comm$only_significant

  # init results folder
  folder_results <- file.path(opts_comm$folder_results, "edgeR")
  dir.create(folder_results, showWarnings=F, recursive=T)

  # used only if save_table is TRUE
  file_table_results <- file.path(folder_results, "table_results.txt")

  # for file names and plot titles
  meta_char <- paste(lapply(metadata(object), function(x) paste(x, collapse="_")), collapse="-")
  design_char <- paste(design, collapse="")

  # 1. build the DGEList object ========================================================================================
  cat("Building the DGEList object ...\n")
  dgel <- load_to_edgeR(object=object)

  # design matrix
  design_matrix <- model.matrix(design, data=colData(object))

  # 2. fit the dispersions and LFCs ====================================================================================

  # equivalent to running sucessively
  #   object <- estimateGLMCommonDisp(object, design)
  #   object <- estimateGLMTrendedDisp(object, design)
  #   object <- estimateGLMTagwiseDisp(object, design)
  #
  # - the trend method may be one of "locfit", "none", "movingave", "loess" and "locfit.mixed"
  # - robust specifies the method used in estimation of prior.df
  dgel <- estimateDisp(dgel, design_matrix,
                       prior.df=NULL,
                       trend.method=disp_trend_method,
                       robust=disp_robust)

  return(T)
}

#' Run limma algorithm.
#' 
#' @author Yoann Pradat
#'
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for
#' RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47. \url{https://doi.org/10.1093/nar/gkv007}
#'
#' Law, C.W., Chen, Y., Shi, W. et al. voom: precision weights unlock linear model analysis tools for RNA-seq read
#' counts. Genome Biol 15, R29 (2014). \url{https://doi.org/10.1186/gb-2014-15-2-r29}
run_limma <- function(object){
  return(T)
}

