# @modified: 17 Nov 2020
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
#' @return a dataframe of results
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula specifying the design for the model matrix of DESeq2. Any variable appearing should be 
#' present in the \code{colData} of object
#' @param contrasts a character vector specifying the contrasts (one or multiple beta coefficient) to 
#' be used for making tests and building results table. 
#' @param opts_algo a named list of options specific to DESeq2
#' @param opts_comm a named list of options common to all methods
#'
#' @import DESeq2
#' @importFrom BiocParallel MulticoreParam
#'
#' @author Yoann Pradat
#'
#' @export
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with
#' DESeq2. Genome Biology, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
run_deseq2 <- function(object, design=NULL, contrasts, opts_algo, opts_comm){

  # 0. options and settings ============================================================================================

  # Options for the function \code{DESeq}
  ncores <- opts_comm$ncores

  # Options for the function \code{results}
  alpha <- opts_comm$alpha
  lfc_test_threshold <- opts_algo$lfc_test_threshold

  # Options for the function \code{lfcShrink}
  lfc_fit_shrink_type <- opts_algo$lfc_fit_shrink_type
  lfc_test_alt_hypothesis <- opts_algo$lfc_test_alt_hypothesis

  # Options for saving
  only_significant <- opts_comm$only_significant

  # init results folder
  folder_results <- file.path(opts_comm$folder_results, "deseq2")
  dir.create(folder_results, showWarnings=F, recursive=T)

  if (opts_comm$save_table) {
    file_table_results <- file.path(folder_results, "table_results.txt")
  } else {
    file_table_results <- NULL
  }

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
  plot_dispersion_deseq2(dds, 
                         filepath=file.path(folder_results, 
                                            paste0("deseq2_dispersion_", meta_char, "_", design_char, ".pdf")),
                         title=paste("Dispersion plot on:", meta_char),
                         subtitle=design_char)

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

  cat("Computing results tables for each for the following contrasts...\n\t")
  cat(paste(lapply(names(contrasts), function(x) paste0(x," = ",contrasts[[x]])), collapse="\n\t"), "\n")
  df_results_all <- data.frame()

  for (contrast_name in names(contrasts)){
    contrast <- contrasts[[contrast_name]]
    contrast_vec <- get_contrast_vector(contrast=contrast, 
                                        design=design, data=as.data.frame(colData(dds)))

    # No LFC shrinkage
    # The statistic is based on the MLE of LFC
    resLFC_MLE <- results(object=dds,
                          contrast=contrast_vec,
                          lfcThreshold=lfc_test_threshold,         # if specifying lfcThreshold, tests are Wald tests
                          altHypothesis=lfc_test_alt_hypothesis,
                          alpha=alpha,                             # for p-value adjustment level
                          filter=mcols(dds)$baseMean,
                          independentFiltering=T,            # for choosing genes whose p-values will be ajdusted
                          pAdjustMethod="BH",
                          format="DataFrame",
                          parallel=T,
                          BPPARAM=BPPARAM)

    # WARNING: the p-values are unchanged, only lfc and lfcSE change (with lfcSE being a posterior SD here).
    resLFC_shrink <- NULL

    if (is.null(lfc_fit_shrink_type)){
      resLFC_shrink <- NULL
    } else if (lfc_fit_shrink_type=="normal"){
      if (any(attr(terms.formula(design(dds)),"order") > 1)){
        cat("Interactions with lfcShrink='normal' is not implemented. Skipping.\n")
      } else {
      # Normal LFC prior is a zero-centered Gaussian (original in DESeq2). Not available with interactions.
      resLFC_shrink <- lfcShrink(dds=dds, 
                                 contrast=contrast_vec,
                                 type=lfc_fit_shrink_type,
                                 format="DataFrame",
                                 parallel=T,
                                 BPPARAM=BPPARAM)
      }
    } else if(lfc_fit_shrink_type=="apeglm"){
      if (sum(contrast_vec==1) != 1){
          cat("Contrasts that are not coefficient names for lfcShrink='apeglm' is not implemented. Skipping.\n")
      } else if (resultsNames(dds)[contrast_vec==1] != rownames(contrast_vec)[contrast_vec==1]){
          cat("Contrasts that are not coefficient names for lfcShrink='apeglm' is not implemented. Skipping.\n")
      } else {
        # Apeglm sets a Cauchy prior with null location
        resLFC_shrink <- lfcShrink(dds=dds, 
                                   coef=resultsNames(dds)[contrast_vec==1],
                                   type=lfc_fit_shrink_type,
                                   format="DataFrame",
                                   parallel=T,
                                   BPPARAM=BPPARAM)
      }
    }

    # NOTE: on results summary
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
  
    # save MA plot
    # "In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of 
    # normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is
    # less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down."

    plot_MA_deseq2(resLFC_MLE=resLFC_MLE,
                   resLFC_shrink=resLFC_shrink,
                   filepath=file.path(folder_results, paste0("deseq2_ma_", meta_char, "_", design_char, ".pdf")), 
                   title=paste("LFC over mean norm counts on:", meta_char),
                   subtitle=paste("design:", design_char, "     contrast:", contrast_name))

    # 4 build contrast-specific table of results =======================================================================

    # statistic (Wald or LRT)
    stat <- resLFC_MLE$stat

    # unadjusted pvalues
    pval <- resLFC_MLE$pvalue
 
    # adjusted pvalues
    padj <- resLFC_MLE$padj

    # beta
    betas <- resLFC_MLE$log2FoldChange

    df_results_all <- make_table_results(data=list(beta=betas, stat=stat, pval=pval, padj=padj),
                                         tags=list(meta=meta_char, design=design_char, contrast=contrast_name),
                                         only_significant=only_significant, alpha=alpha,
                                         df_results_all=df_results_all,
                                         df_row=rowData(object),
                                         file=file_table_results)

  }

  df_results_all 
}

#' Run edgeR algorithm.
#'
#' @return a dataframe of results
#' @param object a \code{SummarizedExperiment} object
#' @param design a formula specifying the design for the model matrix of DESeq2. Any variable appearing should be 
#' present in the \code{colData} of object
#' @param contrasts a character vector specifying the contrasts (one or multiple beta coefficient) to 
#' be used for making tests and building results table. 
#' @param opts_algo a named list of options specific to edgeR
#' @param opts_comm a named list of options common to all methods
#'
#' @import edgeR
#' @importFrom stats p.adjust
#'
#' @author Yoann Pradat
#'
#' @export
#'
#' @references
#' Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of
#' digital gene expression data.” Bioinformatics, 26(1), 139-140. \url{https://doi.org/10.1093/bioinformatics/btp616}
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). “Differential expression analysis of multifactor RNA-Seq experiments with
#' respect to biological variation.” Nucleic Acids Research, 40(10), 4288-4297. \url{https://doi.org/10.1093/nar/gks042}
run_edgeR <- function(object, design=NULL, contrasts, opts_algo, opts_comm){

  # 0. options and settings ============================================================================================

  # Options for dispersion fitting
  disp_trend_method <- opts_algo$disp_trend_method
  disp_robust_prior_df <- opts_algo$disp_robust_prior_df
  disp_robust_prior_ql <- opts_algo$disp_robust_prior_ql

  # Options for lfc fitting and testing
  # Either
  # - QL: then use the functions glmQLFit and glmQLTest. This is the preferred and default option.
  # - LRT: then use the functions glmFit, glmLRT
  lfc_test_type <- match.arg(opts_algo$lfc_test_type, c("QL", "LRT"))
  lfc_test_threshold <- opts_algo$lfc_test_threshold
  lfc_fit_robust_prior_ql_disp <- opts_algo$lfc_fit_robust_prior_ql_disp
  alpha <- opts_comm$alpha

  # Options for saving
  only_significant <- opts_comm$only_significant

  # init results folder
  folder_results <- file.path(opts_comm$folder_results, "edgeR")
  dir.create(folder_results, showWarnings=F, recursive=T)

  # used only if save_table is TRUE
  if (opts_comm$save_table){
    file_table_results <- file.path(folder_results, "table_results.txt")
  } else {
    file_table_results <- NULL
  }

  # for file names and plot titles
  meta_char <- paste(lapply(metadata(object), function(x) paste(x, collapse="_")), collapse="-")
  design_char <- paste(design, collapse="")

  # 1. build the DGEList object ========================================================================================
  cat("Building the DGEList object ...\n")
  dgel <- load_to_edgeR(object=object)

  # design matrix
  design_matrix <- model.matrix(design, data=colData(object))

  # 2. fit the dispersions =============================================================================================

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
                       robust=disp_robust_prior_df)

  # save plot dispersion
  plot_dispersion_edgeR(dgel=dgel, 
                        filepath=file.path(folder_results, 
                                           paste0("edgeR_dispersion_", meta_char, "_", design_char, ".pdf")),
                        title=paste("Dispersion plot on:", meta_char), subtitle=design_char)

  # 3. fit the LFCs and compute test ===================================================================================

  if (lfc_test_type == "QL"){
    # fit
    glm_fit <- glmQLFit(dgel, design_matrix, 
                        winsor.tail.p=c(0.05,0.1),
                        abundance.trend=TRUE,
                        robust=lfc_fit_robust_prior_ql_disp)

    # save plot QL dispersion
    plot_dispersion_ql_edgeR(glm_fit=glm_fit, 
                             filepath=file.path(folder_results, 
                                                paste0("edgeR_dispersion_ql_", meta_char, "_", design_char, ".pdf")),
                             title=paste("QL Dispersion plot on:", meta_char), subtitle=design_char)

  } else if (lfc_test_type == "LRT") {
    # fit
    # TODO: clarify the role of prior.count parametrs
    glm_fit <- glmFit(dgel, design_matrix, 
                      prior.count=0.125)

    # one can extract beta coefficients with
    #   predbeta <- predFC(y=dgel$counts, design_matrix, prior.count=0.125,  offset=getOffset(dgel), 
    #                      dispersion=dgel$tagwise.dispersion)
    # or individually for each beta name after running glmLRT with
    #   predbeta[[beta_name]] <- glm_test$table$logFC

  }


  cat("Computing results tables for each for the following contrasts...\n\t")
  cat(paste(lapply(names(contrasts), function(x) paste0(x," = ",contrasts[[x]])), collapse="\n\t"), "\n")
  df_results_all <- data.frame()

  for (contrast_name in names(contrasts)){
    contrast <- contrasts[[contrast_name]]
    contrast_vec <- get_contrast_vector(contrast=contrast, 
                                        design=design, data=as.data.frame(colData(object)))

    # F-test
    #   - if poisson.bound=T, the pvalue returned will not be less than the pval that would be  obtained for a LRT 
    #     with NB dispersion set to 0
    #
    # for lfc=0 with a glm fit from glmQLFit, glmTreat automatically calls glmLRT i.e
    #   glm_test <- glmQLFTest(glmfit=glm_fit,
    #                          contrast=contrast_vec,
    #                          poisson.bound=T)

    # LRT
    #
    # for lfc=0 with a glm fit from glmFit, glmTreat automatically calls glmLRT i.e
    #   glm_test <- glmLRT(glmfit=glm_fit,
    #                      contrast=contrast_vec)

    glm_test <- glmTreat(glm_fit,
                         contrast=contrast_vec,
                         lfc=lfc_test_threshold)

    # 4 build contrast-specific table of results =======================================================================

    # statistic (LRT or F)
    if (lfc_test_type == "QL"){
      stat <- glm_test$table$F
    } else if (lfc_test_type == "LRT") {
      stat <- glm_test$table$LR
    }

    # unadjusted pvalues
    pval <- glm_test$table$PValue

    # adjusted pvalues
    padj <- p.adjust(pval, method="BH")

    # beta
    #
    # NOTE:
    # in the code of DESeq2Paper, beta coeffs for edgeR are obtained using
    #   log2(exp(1))*glm_fit$coefficients
    # which is equal to
    #   glm_test$table
    # for the particular coefficients tested
    betas <- glm_test$table$logFC

    df_results_all <- make_table_results(data=list(beta=betas, stat=stat, pval=pval, padj=padj),
                                         tags=list(meta=meta_char, design=design_char, contrast=contrast_name),
                                         only_significant=only_significant, alpha=alpha,
                                         df_results_all=df_results_all,
                                         df_row=rowData(object),
                                         file=file_table_results)
  }

  df_results_all
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
