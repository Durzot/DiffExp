suppressMessages(library(DESeq2))

## load bottomly =======================================================================================================

idx <- c(1:3,5:7)

load("../../../tools/DESeq2paper/data/bottomly_sumexp.RData")
ddsB <- DESeqDataSetFromMatrix(
  bottomly@assays$data[["counts"]][,idx], 
  DataFrame(colData(bottomly)[idx,]), 
  ~ strain
)
ddsB <- estimateSizeFactors(ddsB)
ddsB <- estimateDispersions(ddsB)

## plot glm Gamma family ===============================================================================================

pdfGamma <- function(x, mu, phi){
  p <- x^(1/phi-1) * (1/(phi*mu))^(phi^-1) * exp(-x/(phi*mu)) / gamma(1/phi)
  return(p)
}

plotPdfGamma <- function(x, phis, cols){
  par(mar=c(4.5,4.5,1.5,1.5))

  for (i in seq(1, length(phis))){
    phi <- phis[i]
    col <- cols[i]
    print(i)
    if (i == 1){
      plot(x, pdfGamma(x, mu=0.5, phi), type="l", lty=1, col=col, lwd=3, ylab=expression(Gamma(x, mu, phi)), cex=1.5)
    } else {
      lines(x, pdfGamma(x, mu=0.5, phi), type="l", lty=1, col=col, lwd=3)
    }
  }

  text(1.25, 2, expression(Gamma(x, mu, phi) == paste(frac(1, Gamma(phi^-1)), bgroup("(", frac(1, phi * mu), ")")^phi^-1, x^(phi^-1 - 1), e^frac(-x, phi * mu))), cex = 1.75)

  legend("topright",
         legend=sapply(phis, function(x) as.expression(substitute(phi == p, list(p = as.name(x))))),
         pch=20, 
         col=cols,
         bg="white",
         cex=1.5
  )

}

x <- (0:2000)/1000
phis <- c(0.05, 0.1, 0.5, 0.8)
cols <- c("black", "red", "orange", "dodgerblue")


pdf(file = "../img/plotGamma.pdf",
    width = 8,
    height = 6
)
plotPdfGamma(x, phis, cols)
dev.off()

## plot Cauchy distribution ============================================================================================

pdfCauchy <- function(x, mu, gam){
  p <- (gam^2)/((x-mu)^2 + gam^2) * 1/(pi*gam)
  return(p)
}

plotPdfCauchy <- function(x, gams, cols){
  par(mar=c(4.5,4.5,1.5,1.5))

  for (i in seq(1, length(gams))){
    gam <- gams[i]
    col <- cols[i]
    if (i == 1){
      plot(x, pdfCauchy(x, 0, gam), type="l", lty=1, col=col, lwd=3, ylab=expression(Cauchy(x, 0, S)), cex=1.5)
    } else {
      lines(x, pdfCauchy(x, 0, gam), type="l", lty=1, col=col, lwd=3)
    }
  }

  text(0, 4, expression(Cauchy(x, mu, gamma) == paste(frac(1, pi*gamma), bgroup("(", frac(gamma^2, (x-mu)^2 + gamma^2), ")"))), cex = 1.75)

  legend("topright",
         legend=sapply(gams, function(x) as.expression(substitute(gamma == S, list(S = as.name(x))))),
         pch=20, 
         col=cols,
         bg="white",
         cex=1.5
  )

}

x <- (-1000:1000)/1000
gams <- c(0.05, 0.1, 0.5, 0.8)
cols <- c("black", "red", "orange", "dodgerblue")


pdf(file = "../img/plotCauchy.pdf",
    width = 8,
    height = 5
)
plotPdfCauchy(x, gams, cols)
dev.off()

## plot disp ===========================================================================================================

plotDispShrink <- function(dds, show.trend=F, show.MAP=F, show.outliers=F) {
  # pick 40 equally spaced genes along the base mean
  bins <- 10^seq(from=0, to=5, length=20)
  pickone <- function(x) {
    if (sum(x) == 0) return(NULL)
    if (sum(x) == 1) return(which(x))
    sample(which(x), 1)
  }
  up <- sapply(seq_along(bins[-1]), function(i) 
              pickone(mcols(dds)$dispGeneEst > 1e-4 & 
                       !mcols(dds)$dispOutlier & 
                       mcols(dds)$dispGeneEst > mcols(dds)$dispFit & 
                       mcols(dds)$baseMean > bins[i] &
                       mcols(dds)$baseMean < bins[i+1])
  )

  down <- sapply(seq_along(bins[-1]), function(i) 
              pickone(mcols(dds)$dispGeneEst > 1e-4 & 
                      !mcols(dds)$dispOutlier & 
                      mcols(dds)$dispGeneEst < mcols(dds)$dispFit & 
                      mcols(dds)$baseMean > bins[i] & 
                      mcols(dds)$baseMean < bins[i+1])
  )

  # pick 5 outliers
  bins <- 10^seq(from=1, to=4, length=6)
  outliers <- do.call(c, lapply(seq_along(bins[-1]), function(i) 
              pickone(mcols(dds)$dispGeneEst / mcols(dds)$dispMAP > 2 & 
                      mcols(dds)$dispOutlier & mcols(dds)$baseMean > bins[i] & 
                      mcols(dds)$baseMean < bins[i+1]))
  )

  s <- c(up, down, outliers)
  s <- s[!is.na(s)]

  # scatter plot gene-wise estimates
  with(mcols(dds[s,]),
       plot(baseMean, dispGeneEst,log="xy", pch=16, xlab="mean of normalized counts", ylab="dispersion estimate",
            yaxt="n",ylim=c(.001,10))
  )
  axis(2,at=10^(-3:2),label=10^(-3:2))

  legendText <- c("MLE")
  legendCol <- c("black")
  
  if (show.trend){
    # trend function
    xs <- 10^(-20:50/10)
    lines(xs,dispersionFunction(dds)(xs),col="red",lwd=2)
    legendText <- c(legendText, "prior mean")
    legendCol <- c(legendCol, "red")

  }

  if (show.MAP){
    # arrows showing change gene-wise estmites -> MAP
    with(mcols(dds[s,][!mcols(dds[s,])$dispOutlier,]),
         arrows(baseMean, dispGeneEst, baseMean, dispersion,
                length=.075, col="dodgerblue",lwd=2)
    )

    with(mcols(dds[s,][mcols(dds[s,])$dispOutlier,]),
         segments(baseMean, dispGeneEst, baseMean, dispMAP,
                  col="dodgerblue",lwd=2, lty=3)
    )

    legendText <- c(legendText, "MAP")
    legendCol <- c(legendCol, "dodgerblue")
  }

  if (show.outliers){
    # circles around outliers
    with(mcols(dds[s,][mcols(dds[s,])$dispOutlier,]),
         points(baseMean, dispersion,
                cex=2,col="dodgerblue",lwd=2)
    )
  }

  legend("topright",legendText,pch=20, col=legendCol,bg="white")
}

makePlotDispShrink <- function(dds, show.trend=F, show.MAP=F, show.outliers=F){
  set.seed(1)
  plotDispShrink(dds, show.trend, show.MAP, show.outliers)
}

line <- -.1
adj <- -.3
cex <- 1.5

# graph 1

pdf(file = "../img/plotDispShrink_gw_and_trend.pdf",
    width = 10,
    height = 4
)
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5))
makePlotDispShrink(ddsB, show.trend=F, show.MAP=F, show.outliers=F)
mtext("A",side=3,line=line,adj=adj,cex=cex)

makePlotDispShrink(ddsB, show.trend=T, show.MAP=F, show.outliers=F)
mtext("B",side=3,line=line,adj=adj,cex=cex)

dev.off()


## graph 2

pdf(file = "../img/plotDispShrink_gw_trend_and_MAP.pdf",
    width = 10,
    height = 4
)

par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5))
makePlotDispShrink(ddsB, show.trend=T, show.MAP=F, show.outliers=F)
mtext("A",side=3,line=line,adj=adj,cex=cex)

makePlotDispShrink(ddsB, show.trend=T, show.MAP=T, show.outliers=T)
mtext("B",side=3,line=line,adj=adj,cex=cex)
dev.off()

pdf(file = "../img/plotDispShrink_gw_trend_MAP_outliers.pdf",
    width = 8,
    height = 6
)
makePlotDispShrink(ddsB, show.trend=T, show.MAP=T, show.outliers=T)
dev.off()

