#' srnadiff: A package for differential expression of sRNA-Seq.
#'
#' The srnadiff package provides uses four strategies to find differentially
#'     expressed loci.
#'
#' @docType package
#' @name srnadiff
#'
#' @author Matthias Zytnicki, \email{matthias.zytnicki@@inra.fr}
#'
#' @import methods
#' @import devtools
#' @import BiocStyle
#' @importFrom utils read.csv
#' @import S4Vectors
#' @import GenomeInfoDb
#' @import rtracklayer
#' @import SummarizedExperiment
#' @import IRanges
#' @import Rsamtools
#' @import Rsubread
#' @import DESeq2
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import BiocParallel
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom Rcpp evalCpp
#' @useDynLib srnadiff
NULL
