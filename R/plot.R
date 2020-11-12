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

  # add form of the fitted function
  coefs_trend <- attr(dispersionFunction(dds), "coefficients")
  text(x = 1e4,
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
}
