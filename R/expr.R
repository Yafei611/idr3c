#' Gene expression data for human tissues
#'
#' A dataset containing log-fold change of gene expression between 2 human tissues. The fold change were computed based on the expression data generated using Affymetrix array platform (microarray data) and Illumina Hiseq (RNA-seq data). Data from 2 different platform have been normalized and averaged across technical replicates.
#'
#' @format A data frame with 11836 rows and 2 columns:
#' \describe{
#'   \item{Array}{Expression fold change measured using Affymetrix array}
#'   \item{Seq}{Illumina Hiseq}
#' }
#' @source More detailed information can be found in \bold{Seqc/Maqc-Iii Consortium. (2014). A comprehensive assessment of RNA-seq accuracy, reproducibility and information content by the Sequencing Quality Control Consortium. Nature biotechnology, 32(9), 903-914.}. Data processing method can be found in \bold{Lyu, Y., & Li, Q. (2016). A semi-parametric statistical model for integrating gene expression profiles across different platforms. BMC Bioinformatics, 17(1), 51.}. The raw data was obtained from \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47792}
"expr"