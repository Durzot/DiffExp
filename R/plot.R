# @created: 06 Nov 2020
# @modified: 06 Nov 2020
# @author: Yoann Pradat
# 
# DESEQ2 ANALYSIS: Functions for plotting intermediate and final results.

#### DESEQ2 PLOTS ======================================================================================================

#' Plot dispersion estimates from DESeq2
#'
#' @param dds a \code{DESeqDataSet} object
#' @param filepath a character specifying the path to the plot file
#' @param title a character
#' @param subtitle (optional) a character
#'
#' @importFrom graphics mtext par text
#' @importFrom grDevices dev.off pdf
#'
#' @keywords internal
plot_dispersion_deseq2 <- function(dds, filepath, title, subtitle=NULL){
  pdf(file = filepath,
      width = 6,
      height = 6)

  # plot
  ymin <- 1e-2
  ymax <- max(10, ncol(dds))
  DESeq2::plotDispEsts(object = dds,
                       ylim = c(ymin, ymax),
                       ymin = ymin,
                       cex = 0.2)
  px <- mcols(dds)$baseMean
  xmin <- min(px[px>0], na.rm=T)
  xmax <- max(px[px>0], na.rm=T)

  # add form of the fitted function
  coefs_trend <- attr(dispersionFunction(dds), "coefficients")
  text(x = (xmin*xmax)**0.6,
       y = 0.8*ymax, 
       as.expression(substitute(alpha[tr](mu) == A + frac(B, mu), 
                                list(A = as.name(coefs_trend[[1]]),
                                     B = as.name(coefs_trend[[2]])))),
       cex = 1)

  if(!is.null(subtitle)){
    mtext(subtitle, side=3, line=0.2, cex=0.8, font=1L)
  }
  title(title)
  dev.off()

  cat(paste("DESeq2 dispersion plot saved at", filepath),"\n")
}

#' Plot MA from DESeq2
#'
#' @param resLFC_MLE a \code{DESeqResults} object
#' @param resLFC_shrink (optional) a \code{DESeqResults} object
#' @param filepath a character specifying the path to the plot file
#' @param title a character
#' @param subtitle (optional) a character
#' @param alpha a float passed to \code{\link[DESeq2]{plotMA}}
#' @param ylim a length-2 float vector passed to \code{\link[DESeq2]{plotMA}}
#'
#' @importFrom graphics mtext par text
#' @importFrom grDevices dev.off pdf
#'
#' @importFrom DESeq2 plotMA
#'
#' @keywords internal
plot_MA_deseq2 <- function(resLFC_MLE, resLFC_shrink=NULL, filepath, title, subtitle=NULL, alpha=0.1, ylim=c(-2,2)){
  if (is.null(resLFC_shrink)){
    pdf(file=filepath,
        width=8,
        height=4
    )

    # 1 row, 1 column + margins
    par(mfrow=c(1,1), oma=c(0,0,2,0), mar=c(3,3,2,2))

    DESeq2::plotMA(object=resLFC_MLE ,
           alpha=alpha,
           ylim=ylim
    )
    title("No shrinkage", font=3, line=1, cex=1)

    if(!is.null(subtitle)){
      mtext(title, font=2, side=3, line=1, outer=T, cex=1.2)
      mtext(subtitle, side=3, line=2, cex=0.8, font=1L)
    } else {
      mtext(title, font=2, side=3, line=0.5, outer=T, cex=1.2)
    }
  } else {
    pdf(file=filepath,
        width=8,
        height=8
    )

    # 2 rows, 1 column + margins
    par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(3,3,2,2))

    # A (no shrink)
    DESeq2::plotMA(object=resLFC_MLE ,
           alpha=alpha,
           ylim=ylim
    )
    title("No shrinkage", font=3, line=1, cex=1)

    # title(s)
    if(!is.null(subtitle)){
      mtext(title, font=2, side=3, line=1, outer=T, cex=1.2)
      mtext(subtitle, side=3, line=2, cex=0.8, font=1L)
    } else {
      mtext(title, font=2, side=3, line=0.5, outer=T, cex=1.2)
    }

    DESeq2::plotMA(object=resLFC_shrink,
           alpha=alpha,
           ylim=c(-2,2)
    )
    title("Shrinkage", font=3, cex=1)
  }

  dev.off()
  cat(paste("DESeq2 MA plot saved at", filepath),"\n")
}

#### edgeR PLOTS ======================================================================================================

#' Plot dispersion estimates from edgeR
#'
#' @param dgel a \code{DGEList} object
#' @param filepath a character specifying the path to the plot file
#' @param title a character
#' @param subtitle (optional) a character
#'
#' @importFrom graphics mtext par text
#' @importFrom grDevices dev.off pdf
#'
#' @keywords internal
plot_dispersion_edgeR <- function(dgel, filepath, title, subtitle=NULL){
  pdf(file = filepath,
      width = 6,
      height = 6)

  xmin <- min(dgel$AveLogCPM, na.rm=T)
  xmax <- max(dgel$AveLogCPM, na.rm=T)

  ymin <- min(dgel$tagwise.dispersion, na.rm=T)
  ymax <- max(dgel$tagwise.dispersion, na.rm=T)

  # plot
  edgeR::plotBCV(y=dgel,
                 cex=0.6,
                 col.common="green",
                 col.trend="red",
                 col.tagwise="black")

  # add annotations
  text(x = 0.5*xmax + 0.5*xmin,
       y = 0.90*ymax + 0.1*ymin,
       bquote("Common BCV:" ~ hat(phi)^frac(1,2) == .(sqrt(dgel$common.dispersion))),
       cex = 0.8)

  # add annotations
  text(x = 0.5*xmax + 0.5*xmin,
       y = 0.85*ymax + 0.15*ymin,
       paste("Trend meth:", dgel$trend.method),
       cex = 0.8)

  if(!is.null(subtitle)){
    mtext(subtitle, side=3, line=0.2, cex=0.8, font=1L)
  }
  title(title)
  dev.off()
  cat(paste("edgeR dispersion plot saved at", filepath),"\n")
}

#' Plot dispersion estimates from edgeR
#'
#' @param glm_fit a \code{DGEGLM} object
#' @param filepath a character specifying the path to the plot file
#' @param title a character
#' @param subtitle (optional) a character
#'
#' @importFrom graphics mtext par text
#' @importFrom grDevices dev.off pdf
#'
#' @keywords internal
plot_dispersion_ql_edgeR <- function(glm_fit, filepath, title, subtitle=NULL){
  pdf(file = filepath,
      width = 6,
      height = 6)

  # plot
  edgeR::plotQLDisp(glmfit=glm_fit,
                    cex=0.6,
                    col.raw="black",
                    col.shrunk="dodgerblue",
                    col.trend="red")

  if(!is.null(subtitle)){
    mtext(subtitle, side=3, line=0.2, cex=0.8, font=1L)
  }
  title(title)
  dev.off()
  cat(paste("edgeR ql dispersion plot saved at", filepath),"\n")
}
