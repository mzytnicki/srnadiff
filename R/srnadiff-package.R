###############################################################################
#
# srnadiff-package organization of R files
#
# AllClasses .... class definitions and object constructors
# AllGenerics ... the generics defined in srnadiff-package
# methods ....... the S4 methods (accessors for slots and replace methods)
# readAnnot ..... reads and parses GFF/GTF files
# srnadiff ...... core function to find differentially expressed regions
# naive \
# hmm    ........ segmentation method functions
# IR    /
# pvalues ....... method that compute the p-values
# plotRegions ... plot function
# helper ........ computeNormFactors, computeLogFoldChange, cvgNormalization,
#                 reconcileRegions, checkParameters, IRlist2GR
# RcppExports ... the R wrappers for the C++ functions (auto)
#
#
# General outline of the internal function calls.
# Note: not all of these functions are exported.
#
# --------------------------------
# srnadiffExp
# |- computeNormFactors
# --------------------------------
#
# --------------------------------
# srnadiff
# |- srnadiffCore
#    |- runNaive
#       |- computeLogFoldChange
#          |- cvgNormalization
#       |- IRlist2GR
#    |- runIR
#       |- computeLogFoldChange
#          |- cvgNormalization
#       |- rcpp_ir (C++)
#    |- runHmm
#       |- buildDataHmm
#          |- rcpp_buildHmm (C++)
#       |- computePvalues
#          |- DESeq or
#          |- edgeR
#       |- hmm
#          |- rcpp_viterbi (C++)
#    |- reconcileRegions
#       |- computePvalues
# --------------------------------
#
# --------------------------------
# plotRegions
# |- cvgNormalization
# --------------------------------
#
###############################################################################

#' Finding differentially expressed unannotated genomic regions from
#' RNA-seq data
#'
#' \code{srnadiff} is a package that finds differently expressed regions
#' from RNA-seq data at \emph{base-resolution} level without relying on
#' existing annotation. To do so, the package implements the
#' \emph{identify-then-annotate} methodology that builds on the idea of
#' combining two pipelines approach: differential expressed regions detection
#' and differential expression quantification.
#'
#' The \code{srnadiff} package implements two major methods to produce
#' potential differentially expressed regions: the HMM and IR method.
#' Briefly, these methods identify contiguous base-pairs in the genome
#' that present differential expression signal, then these are regrouped
#' into genomic intervals called differentially expressed regions (DERs).
#'
#' Once DERs are detected, the second step in a sRNA-diff approach is to
#' quantify the statistic signification of these. To do so, reads (including
#' fractions of reads) that overlap each expressed region are counted to
#' arrive at a count matrix with one row per region and one column per sample.
#' Then, this count matrix is analyzed using the standard workflow of
#' \code{DESeq2} for differential expression of RNA-seq data, assigning a
#' p-value to each candidate DER. Alternatively, \code{edgeR} can be used.
#'
#' The main functions for finds differently expressed regions are
#' \code{\link{srnadiffExp}} and \code{\link{srnadiff}}. The first one
#' creats an S4 class providing the infrastructure (slots) to store the
#' input data, methods parameters, intermediate calculations and results
#' of an sRNA-diff approach. The second one implement four methods to find
#' candidate differentially expressed regions and quantify the statistic
#' signification of the finded regions. Details about the implemented methods
#' are further described in the vignette and the manual page of the
#' \code{\link{srnadiff}} function.

#'
#' @examples
#' ## A typical srnadiff session might look like the following.
#'
#' ## Here we assume that 'bamFiles' is a vector with the full
#' ## paths to the BAM files and the sample and experimental
#' ## design information are stored in a data frame 'sampleInfo'.
#'
#' \dontrun{
#'
#' #-- Data preparation
#' srnaExp <- srnadiffExp(bamFiles, sampleInfo)
#'
#' #-- Detecting DERs and quantifying differential expression
#' srnaExp <- srnadiff(srnaExp)
#'
#' #-- Visualization of the results
#' plotRegions(srnaExp, regions(srnaExp)[1])
#' }
#'
#' @docType package
#' @name srnadiff
#' @aliases srnadiff-package
#'
#' @author Matthias Zytnicki and Ignacio González
#'
#' @import stats
#' @import methods
#' @import devtools
#' @import S4Vectors
#' @import GenomeInfoDb
#' @import rtracklayer
#' @import SummarizedExperiment
#' @import IRanges
#' @import Rsamtools
#' @import DESeq2
#' @import edgeR
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import BiocParallel
#' @import GenomicRanges
#' @import Gviz
#' @import BiocStyle
#' @importFrom BiocManager version
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices col2rgb
#' @useDynLib srnadiff
#'
#' @keywords package
NULL
