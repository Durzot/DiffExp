#' RNA-seq and clinical data for TCGA BRCA Stage I patients.
#'
#' Data from the reanalysis of RNA-seq files by Zheng et al. 2019. The clinical data were taken from the GDC data
#' portal and the cBio data portal. Both the RNA-seq data and the clinical data are public access. The samples were
#' chosen so that
#'  1. all samples are Stage I
#'  2. all samples are Female
#'  3. all samples have complete clinical data
#'
#' This selection was performed only to keep this example data set small and may not be suitable for meaningful
#' analyses.
#' 
#' @format A \code{SummarizedExperiment} object with 58,288 genes and 98 samples.
#' \describe{
#'   \item{assays}{A \code{SimpleList} with just one "counts" matrix}
#'   \item{colData}{A \code{DFrame} with 98 rows and 8 columns}
#'   \item{rowData}{A \code{DFrame} with 58,288 rows and 2 columns}
#'   \item{metadata}{An empty list}
#'   \item{elementMetadata}{Alias for \code{rowData}}
#' }
#'
#' @docType data
#' @usage data(tcga_brca)
#' @keywords datasets
#'
#' @references 
#' Hong Zheng, Kevin Brennan, Mikel Hernaez, Olivier Gevaert, Benchmark of long non-coding RNA quantification for RNA
#' sequencing of cancer samples, GigaScience, Volume 8, Issue 12, December 2019, giz145,
#' <https://doi.org/10.1093/gigascience/giz145>
#'
#' @source RNA-seq source: <https://stanfordmedicine.app.box.com/s/lu703xuaulfz02vgd2lunxnvt4mfvo3q>. cBioPortal
#' BRCA Pan-cancer Atlas <https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018>. GDC legacy archive
#' clinical BRCA <https://portal.gdc.cancer.gov/legacy-archive/files/735bc5ff-86d1-421a-8693-6e6f92055563>
#'
#' @examples
#' data(tcga_brca)
"tcga_brca"
