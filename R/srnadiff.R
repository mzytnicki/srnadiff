#' Find differentially expressed sRNA regions
#'
#' This is the main wrapper for running several key functions from this
#' package.  It is meant to be used after that a \code{\link{srnadiffExp}}
#' object has been created. \code{\link{srnadiff}} implement four methods to
#' produce potential DERs (see Details).
#' Once DERs are detected, the second step in \code{\link{srnadiff}} is to
#' quantify the statistic signification of these.
#'
#' Implemented methods to produce potential differentially expressed
#' regions in \code{\link{srnadiff}} are:
#'
#' \describe{
#' \item{\code{annotation:}}{This method simply provides the genomic regions
#'     corresponding to the annotation file that is optionally given by the
#'     user. It can be a set of known miRNAs, siRNAs, piRNAs, genes, or a
#'     combination thereof.}
#'
#' \item{\code{hmm:}}{This approach assumes that continuous regions of RNA
#'       along the chromosome are either "differentially expressed" or "not".
#'       This is captured with a hidden Markov model (HMM) with binary latent
#'       state of each nucleotide: \emph{differentially expressed} or
#'       \emph{not differentially expressed}. The observations of the HMM are
#'       then the empirical p-values arising from the differential expression
#'       analysis corresponding to each nucleotide position.
#'       The HMM approach normally needs emission, transition, and starting
#'       probabilities values (see \code{\link{parameters}}). They
#'       can be tuned by the user. In order to finding the most likely sequence
#'       of states from the HMM, the Viterbi algorithm is performed. This
#'       essentially segments the genome into regions, where a region is
#'       defined as a set of consecutive bases showing a common expression
#'       signature.}
#'
#' \item{\code{IR:}}{In this approach, for each base, the average from
#'       the normalized coverage is calculated across all samples into each
#'       condition. This generates a vector of (normalized) mean coverage
#'       expression per condition. These two vectors
#'       are then used to compute per-nucleotide log-ratios (in absolute value)
#'       across the genome. For the computed log-ratio expression, the
#'       method uses a sliding threshold \emph{h} that run across the log-ratio
#'       levels identifying bases with log-ratio value above of \emph{h}.
#'       Regions of contiguous bases passing this threshold are then analyzed
#'       using an adaptation of Aumann and Lindell algorithm for irreducibility
#'       property (Aumann and Lindell (2003)).}
#'
#' \item{\code{naive:}}{This method is the simplest, gived a fixed threshold
#'       \emph{h}, contiguous bases with log-ratio expression
#'       (in absolute value) passing this threshold are then considered as
#'       candidate differentially expressed regions.}
#' }
#'
#' @references
#' Aumann Y. and, Lindell Y. (2003). A Statistical Theory for Quantitative
#' Association Rules. \emph{Journal of Intelligent Information Systems},
#' 20(3):255-283.
#'
#' @param object        An \code{\link{srnadiffExp}} object.
#' @param segMethod     A character vector. The segmentation methods to use,
#'                      one of \code{'annotation'}, \code{'naive'},
#'                      \code{'hmm'}, \code{'IR'} or combinations thereof.
#'                      Default \code{'all'}, all methods are used. See
#'                      Details.
#' @param diffMethod    A character. The differential expression testing
#'                      method to use, one of \code{'DESeq2'}, \code{'edgeR'}.
#'                      See Details.
#' @param nThreads      \code{integer(1)}. Number of workers.
#'                      Defaults to all cores available as determined by
#'                      \code{\link[BiocParallel]{multicoreWorkers}}.
#' @param useParameters A named list containing the methods parameters to use.
#'                      If missing, default parameter values are supplied.
#'                      See \code{\link{parameters}} for details.
#'
#' @return An \code{srnadiffExp} object containing additional slots for:
#' \itemize{
#'     \item \code{regions}
#'     \item \code{parameters}
#'     \item \code{countMatrix}
#' }
#'
#' @seealso
#' \code{\link{regions}}, \code{\link{parameters}}, \code{\link{countMatrix}}
#' and \code{\link{srnadiffExp}}
#'
#' @examples
#' srnaExp <- srnadiffExample()
#' srnaExp <- srnadiff(srnaExp)
#' srnaExp
#'
#' @export
srnadiff <- function(object,
                     segMethod=c("hmm", "IR"),
                     diffMethod="DESeq2",
                     useParameters=srnadiffDefaultParameters,
                     nThreads=1) {

    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    ##- object
    if (!is(object, "srnadiffExp")) {
        stop("'object' must be an object of class 'srnadiffExp'.", call.=FALSE)
    }

    ##- segMethod
    choices <- c("all", "annotation", "naive", "hmm", "IR")
    segMethod <- choices[pmatch(segMethod, choices)]

    if (any(is.na(segMethod))) {
        stop("'segMethod' should be 'all' or one of 'annotation', 'naive'",
            " 'hmm', 'IR' or combinations thereof.", call.=FALSE)
    }

    if ("all" %in% segMethod) {
        segMethod <- c("annotation", "naive", "hmm", "IR")
    }

    ##- diffMethod
    if (is.null(diffMethod) || !is.character(diffMethod) ||
        (length(diffMethod) != 1)) {
        stop("Invalid declaration off 'diffMethod'. It must be a character",
             " of size 1.", call.=FALSE)
    }
    choices <- c("deseq2", "edger")
    if (!(tolower(diffMethod) %in% choices)) {
        stop("'diffMethod' should be 'DESeq2', or 'edgeR'. ",
             "Got '", diffMethod, "'.", call.=FALSE)
    }
    diffMethod <- tolower(diffMethod)

    ##- nThreads
    if (is.null(nThreads) || !is.numeric(nThreads) || (nThreads < 1) ||
        !is.finite(nThreads)) {
        stop("invalid number of threads, 'nThreads'. It must be a positive",
             " integer.", call.=FALSE)
    }

    nThreads <- nThreads - trunc(nThreads)

    if (nThreads > 0) {
        stop("invalid number of threads, 'nThreads'. It must be a positive",
            " integer.", call.=FALSE)
    }

    ##- useParameters
    defaultparnames <- names(srnadiffDefaultParameters)

    if (is.null(object@parameters)) {
        if (missing(useParameters)) {
            parameters(object) <- srnadiffDefaultParameters
        } else {
            if (!is(useParameters, "list")) {
                stop("'useParameters' must be a named list. See",
                     " help(parameters) for details.", call.=FALSE)
            }

            valueNames <- names(useParameters)

            if (any(duplicated(valueNames))) {
                stop("duplicate name parameters in 'useParameters'. See",
                     " help(parameters) for details.", call.=FALSE)
            }

            if (!all(valueNames %in% defaultparnames)) {
                stop("'useParameters' must be a named list of valid",
                     " parameters. See help(parameters) for details.",
                     call.=FALSE)
            }

            parameters(object) <- useParameters
        }
    }

    ##- end checking ---------------------------------------------------------#

    object@diffMethod <- diffMethod
    args = c(list(object     = object,
                  segMethod  = segMethod,
                  nThreads   = nThreads),
            parameters(object))

    do.call('srnadiffCore', args)
}


##- srnadiff core function ---------------------------------------------------#
##----------------------------------------------------------------------------#
srnadiffCore <- function(object,
                         segMethod=c("annotation", "naive", "hmm", "IR"),
                         diffMethod="DESeq2",
                         nThreads=1,
                         minDepth=10,
                         minSize=18,
                         maxSize=1000000,
                         minGap=100,
                         minOverlap=10,
                         noDiffToDiff=0.001,
                         diffToNoDiff=0.000001,
                         emission=0.9,
                         emissionThreshold=0.1,
                         minLogFC=0.5) {

    if (nThreads > 1) {
        register(MulticoreParam(nThreads))
    }

    allRegions <- list()

    ##- run annotation -------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if (("annotation" %in% segMethod) & !is.null(object@annotReg)) {
        allRegions <- do.call("c", list(allRegions, object@annotReg))
    }

    ##- run naive  -----------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if ("naive" %in% segMethod) {
        allRegions <- do.call("c", list(allRegions, runNaive(object)))
    }

    ##- run IR ---------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if ("IR" %in% segMethod) {
        allRegions <- do.call("c", list(allRegions, runIR(object)))
    }

    ##- run hmm --------------------------------------------------------------#
    ##------------------------------------------------------------------------#
    if ("hmm" %in% segMethod) {
        allRegions <- do.call("c", list(allRegions, runHmm(object, nThreads)))
    }

    allRegions <- unique(sort(do.call("c", allRegions)))
    res <- reconcileRegions(object, allRegions, minOverlap)
    object@regions <- res$regions
    object@countMatrix <- res$countMatrix

    return(object)
}
