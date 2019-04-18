###############################################################################
### srnadiffExp S4 class definition
###############################################################################
#' Infrastructure for sRNA-Seq experiment and differential expression
#'
#' \code{srnadiffExp} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a sRNA-diff approach.
#'
#' @details \code{srnadiffExp} load and summarize sample BAM
#' files into base-resolution coverage and estimate the size factors (the
#' effective library size) from the coverage data.
#'
#' To facilitate programming pipelines, \code{NULL} values are input
#' for \code{regions}, \code{parameters} and \code{countMatrix} slots,
#' in which case the default value is used as if the argument had been
#' missing. These slots will be updated after differential expression
#' (\code{\link{srnadiff}}) approach.
#'
#' @name srnadiffExp
#' @rdname srnadiffExp
#' @docType class
#' @aliases srnadiffExp srnadiffExp-class
#'
#' @slot bamFiles        A \code{\link[Rsamtools]{BamFileList}} object
#'                       with the full paths to the BAM files.
#' @slot sampleInfo      A \code{data.frame} with sample and experimental
#'                       design information. Each row describes one sample.
#' @slot annotReg        A \code{GRanges} with annotation information.
#' @slot chromosomeSizes A named vector with the sizes of the chromosomes.
#' @slot coverages       The sample coverages, a named
#'                       \code{\link[IRanges]{RleList}} object.
#' @slot normFactors     A vector of normalization factors.
#' @slot regions         A \code{GenomicRanges} of the candidate
#'                       differentially expressed regions.
#' @slot countMatrix     A matrix of non-negative integer count values,
#'                       one row per region and one column per sample.
#' @slot parameters      An named \code{list}. The parameters for the
#'                       segmentation methods. See \code{\link{parameters}}.
#'
setClass("srnadiffExp", slots=c(bamFiles="ANY",
                                sampleInfo="ANY",
                                annotReg="ANY",
                                chromosomeSizes="ANY",
                                coverages="ANY",
                                normFactors="ANY",
                                regions="ANY",
                                countMatrix="ANY",
                                parameters="ANY")
)


##- srnadiffExp S4 class constructor -----------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname srnadiffExp
#' @docType class
#'
#' @param bamFiles    A vector with the full paths to the BAM files.
#' @param sampleInfo  A \code{data.frame} with three columns labelled
#'                    \code{FileName}, \code{SampleName} and \code{Condition}.
#'                    The first column is the BAM file name (without extension),
#'                    the second column the sample name, and the third column
#'                    the condition to which sample belongs. Each row describes
#'                    one sample.
#' @param annotReg    Optional annotation information. Annotated regions as a
#'                    \code{GRanges} object. By example, ranges in the output
#'                    from \code{\link{readAnnotation}}.
#' @param normFactors A numeric vector, one size factor for each sample in the
#'                    data.
#'
#' @return \code{srnadiffExp} constructor returns an \code{srnadiffExp}
#'         object of class S4.
#'
#' @examples
#' basedir    <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' sampleInfo <- read.csv(file.path(basedir, "dataInfo.csv"))
#' gtfFile    <- file.path(basedir, "Homo_sapiens.GRCh38.76.gtf.gz")
#' annotReg   <- readAnnotation(gtfFile, feature="gene", source="miRNA")
#' bamFiles   <- paste(file.path(basedir, sampleInfo$FileName), "bam", sep = ".")
#'
#' srnaExp <- srnadiffExp(bamFiles, sampleInfo, annotReg)
#' srnaExp
#'
#' @export
srnadiffExp <- function(bamFiles=NULL,
                        sampleInfo=NULL,
                        annotReg=NULL,
                        normFactors=NULL) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#

    ##- bamFiles
    if (is.null(bamFiles)) {
        stop("'bamFiles' must be a vector with the full paths to the BAM",
        " files.", call.=FALSE)
    }

    if (!is.vector(bamFiles)) {
        stop("'bamFiles' must be a vector.", call.=FALSE)
    }

    ##- sampleInfo
    if (is.null(sampleInfo)) {
        stop("'sampleInfo' must be a 'data.frame'.", call.=FALSE)
    }

    if (class(sampleInfo) != "data.frame") {
        stop("'sampleInfo' must be a 'data.frame'.", call.=FALSE)
    }

    if (ncol(sampleInfo) != 3) {
        stop("'sampleInfo' must be a 'data.frame' with three columns.",
            call.=FALSE)
    }

    tmp <- (colnames(sampleInfo) == c("FileName", "SampleName", "Condition"))

    if (sum(tmp) != 3) {
        stop("'sampleInfo' must be a 'data.frame' with three columns",
            " labelled: 'FileName', 'SampleName' and 'Condition'.",
            call.=FALSE)
    }

    if (any(duplicated(sampleInfo$FileName))) {
        stop("non-unique file names in 'sampleInfo$FileName'", call.=FALSE)
    }

    if (any(duplicated(sampleInfo$SampleName))) {
        stop("non-unique sample names in 'sampleInfo$SampleName'", call.=FALSE)
    }

    if (length(unique(sampleInfo$Condition)) != 2) {
        stop("only two conditions are possible in 'srnadiff'.", call.=FALSE)
    }

    ##- compatibility between bamFiles and sampleInfo
    nBams  <- length(bamFiles)
    nReps  <- nrow(sampleInfo)

    if (nBams != nReps) {
        stop("The number of input BAM files should be equal to the number",
            " of samples (", nBams, " and ", nReps, " resp.).", call.=FALSE)
    }

    bamNames <- strsplit(basename(bamFiles), ".", fixed = TRUE)
    bamNames <- unlist(lapply(bamNames, function(x) x[1]))

    if (!(all(as.character(sampleInfo$FileName) %in% bamNames))) {
        stop("'FileName' in sampleInfo do not match BAM file names in",
            " bamFiles.", call.=FALSE)
    }

    ##- annotReg
    if (!is.null(annotReg)) {
        if (!is(annotReg, "GRanges")) {
            stop("'annotReg' must be a 'GRanges' object.", call.=FALSE)
        }
    }

    ##- normFactors
    if (!is.null(normFactors)) {
        n <- nrow(sampleInfo)

        if (length(normFactors) != n) {
            stop("'normFactors' must be a vector of length ", n, ".",
                call.=FALSE)
        }

        if (any(is.na(normFactors))) {
            stop("NA values in 'normFactors'.", call.=FALSE)
        }

        if (!all(is.finite(normFactors))) {
            stop("Infinite values in 'normFactors'.", call.=FALSE)
        }

        if (!all(normFactors > 0)) {
            stop("'normFactors' must be a vector of positive values.",
                call.=FALSE)
        }

        factorNames <- names(normFactors)
        sampleName <- sampleInfo$SampleName

        if (is.null(factorNames)) {
            names(normFactors) <- sampleName
        } else {
            if (any(factorNames != sampleName)) {
                stop("'normFactors' names must be to match the",
                     " 'SampleName' column in sampleInfo.", call.=FALSE)
            }
        }
    }

    ##- end checking ---------------------------------------------------------#

    message("Constructing object...")
    indexFileNames <- paste0(bamFiles, ".bai")
    unIndexedFiles <- bamFiles[!file.exists(indexFileNames)]

    if (length(unIndexedFiles) > 0) { indexBam(unIndexedFiles) }

    bamFiles <- BamFileList(lapply(bamFiles,
                                   function (b) {
                                        BamFile(b, yieldSize=500000,
                                                paste0(b, ".bai"))
                                    }))

    object <- new("srnadiffExp")

    object@annotReg <- annotReg
    object@bamFiles <- bamFiles
    object@sampleInfo <- sampleInfo
    object@chromosomeSizes <- seqlengths(bamFiles[[1]])
    object@coverages <- lapply(bamFiles, coverage)

    if (is.null(normFactors)) {
        normFactors <- computeNormFactors(object@coverages)
    }

    names(normFactors) <- sampleInfo$SampleName
    object@normFactors <- normFactors
    object@regions <- NULL
    object@countMatrix <- NULL
    object@parameters <- NULL

    message("... done.")

    return(invisible(object))
}


##- Example constructor ------------------------------------------------------#
##----------------------------------------------------------------------------#
#' Example constructor
#'
#' This function provides an example of a \code{srnadiffExp} object from
#' two conditions ... à completer ...
#'
#' Raw data have been downloaded from the GEO data set GSE62830, provided in
#' Viollet et \emph{al}. (2015). Adapters were removed with
#' fastx_clipper and mapped with bowtie2 (Salzberg and Langmead, 2012) on the
#' human genome version GRCh38.
#'
#' This example is restricted to a small locus on chr14.
#' It uses the whole genome annotation (with coding genes, etc.) and extracts
#' miRNAs.
#'
#' @return An \code{srnadiffExp} object called '\code{srnaExp}'.
#'
#' @references
#' Viollet, Coralie, David A. Davis, Martin Reczko, Joseph M. Ziegelbauer,
#' Francesco Pezzella, Jiannis Ragoussis, and Robert Yarchoan (2015).
#' "Next-Generation Sequencing Analysis Reveals Differential Expression
#' Profiles of Mirna-mRNA Target Pairs in Kshv-Infected Cells."
#' \emph{PLOS ONE}, 10:1–23.
#'
#' Salzberg, Steven, and Ben Langmead (2012). "Fast gapped-read alignment
#' with Bowtie 2." \emph{Nature Methods}, 9:357–59.
#'
#' @examples
#' ## The 'srnadiffExp' object in this example was constructed by:
#'
#' \dontrun{
#'
#' basedir    <- system.file("extdata", package="srnadiff", mustWork = TRUE)
#' sampleInfo <- read.csv(file.path(basedir, "dataInfo.csv"))
#' gtfFile    <- file.path(basedir, "Homo_sapiens.GRCh38.76.gtf.gz")
#' annotReg   <- readAnnotation(gtfFile, feature="gene", source="miRNA")
#' bamFiles   <- paste(file.path(basedir, sampleInfo$FileName), "bam", sep = ".")
#' srnaExp    <- srnadiffExp(bamFiles, sampleInfo, annotReg)
#' }
#'
#' srnaExp <- srnadiffExample()
#' srnaExp
#'
#' @export
srnadiffExample <- function() {
    srnaExp <- NULL
    basedir <- system.file("extdata", package="srnadiff", mustWork=TRUE)
    load(file.path(basedir, "srnadiffExample.rda"))
    ind <- setdiff(grep(".bam", dir(basedir)), grep(".bai", dir(basedir)))
    bamFiles <- file.path(basedir, dir(basedir)[ind])
    indexFileNames <- paste0(bamFiles, ".bai")
    unIndexedFiles <- bamFiles[!file.exists(indexFileNames)]

    if (length(unIndexedFiles) > 0) { indexBam(unIndexedFiles) }

    bamFiles <- BamFileList(lapply(bamFiles,
                                   function (b) {
                                       BamFile(b, yieldSize=500000,
                                               paste0(b, ".bai"))
                                   }))

    srnaExp@bamFiles <- bamFiles
    return(invisible(srnaExp))
}
