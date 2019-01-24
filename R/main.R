#' Overloading the show method
#' @param object An \code{srnadiff} object.
#' @return       A description of the object.
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp
#'
#' @export
setMethod(  f         ="show",
            signature ="sRNADiff",
            definition=function (object) {
                cat("Object of class sRNADiff.\n",
                    "Annotation (length = ",
                    length(object@annotation),
                    ")\nBAM files (length = ",
                    length(object@bamFileNames),
                    "):\n",
                    paste(' ', head(object@bamFileNames, 3)),
                    "...\nChromosomes (length = ",
                    length(object@chromosomes),
                    "):\n",
                    paste(' ', head(object@chromosomes, 3)),
                    "...\nReplicates (length = ",
                    length(object@replicates),
                    "):\n",
                    paste(' ', head(object@replicates, 3)),
                    "...\nConditions (length = ",
                    length(object@conditions),
                    "):\n",
                    paste(' ', head(object@conditions, 3)),
                    "...\nDifferentially expressed regions (length = ",
                    length(object@regions),
                    ")\nSkip annotation step: ",
                    object@skipAnnotation,
                    "\nSkip naive step: ",
                    object@skipNaive,
                    "\nSkip HMM step: ",
                    object@skipHmm,
                    "\nSkip slice step: ",
                    object@skipSlice,
                    "\n", sep = "")
            }
)

#' Get the output regions
#' @rdname regions-method
#' @param  object  An \code{srnadiff} object.
#' @param  pvalue  A minimum p-value
#' @return         The selected regions
#'
#' @examples
#' exp <- sRNADiffExample()
#' regions(exp)
#'
#' @export
setGeneric( name="regions",
            def =function(object, pvalue = 0.05) {
                standardGeneric("regions")
            }
)

#' @rdname  regions-method
#' @export
setMethod(  f        ="regions",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, pvalue) {
                return(object@regions[mcols(object@regions)$padj <= pvalue])
            }
)

#' @rdname  regions-method
#' @export
setMethod(  f        ="regions",
            signature= c("sRNADiff"),
            definition=function(object) {
                return(object@regions)
            }
)


#' Set the different steps
#' @rdname setStrategies-method
#' @param  object     An \code{srnadiff} object.
#' @param  annotation The annotation step.
#' @param  naive      The naive step.
#' @param  hmm        The HMM step.
#' @param  slice      The slice step.
#' @return         The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setStrategies(exp, TRUE, FALSE, TRUE, TRUE)
#'
#' @export
setGeneric( name="setStrategies",
            def =function(object, annotation, naive, hmm, slice) {
                standardGeneric("setStrategies")
            }
)

#' @rdname  setStrategies-method
#' @export
setMethod(  f        ="setStrategies",
            signature= c("sRNADiff", "logical", "logical", "logical",
                            "logical"),
            definition=function(object, annotation, naive, hmm, slice) {
                                object@skipAnnotation <- ! annotation
                                object@skipNaive      <- ! naive
                                object@skipHmm        <- ! hmm
                                object@skipSlice      <- ! slice
                return(object)
            }
)


#' Set min and max sizes of the regions
#' @rdname setSizes-method
#' @param  object  An \code{srnadiff} object.
#' @param  minSize The minimum size.
#' @param  maxSize The maximum size.
#' @return         The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setSizes(exp, 10, 1000)
#' regions(exp)
#'
#' @export
setGeneric( name="setSizes",
            def =function(object, minSize, maxSize) {
                standardGeneric("setSizes")
            }
)

#' @rdname  setSizes-method
#' @export
setMethod(  f        ="setSizes",
            signature= c("sRNADiff", "numeric", "numeric"),
            definition=function(object, minSize, maxSize) {
                                object@minSize <- minSize
                                object@maxSize <- maxSize
                return(object)
            }
)


#' Set min minimum depth to localize regions
#' @rdname setMinDepth-method
#' @param  object An \code{srnadiff} object.
#' @param  depth  The minimum depth
#' @return        The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setMinDepth(exp, 3)
#'
#' @export
setGeneric( name="setMinDepth",
            def =function(object, depth) {
                standardGeneric("setMinDepth")
            }
)

#' @rdname  setMinDepth-method
#' @export
setMethod(  f        ="setMinDepth",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, depth) {
                                object@minDepth <- depth
                return(object)
            }
)


#' Set the threshold to merge close regions (in the naive step)
#' @rdname setMergeDistance-method
#' @param  object   An \code{srnadiff} object.
#' @param  distance The maximum distance
#' @return          The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setMergeDistance(exp, 1000)
#'
#' @export
setGeneric( name="setMergeDistance",
            def =function(object, distance) {
                standardGeneric("setMergeDistance")
            }
)

#' @rdname  setMergeDistance-method
#' @export
setMethod(  f        ="setMergeDistance",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, distance) {
                                object@mergeDistance <- distance
                return(object)
            }
)


#' Set the threshold to remove similar regions (in the slice step)
#' @rdname setMinDifferences-method
#' @param  object      An \code{srnadiff} object.
#' @param  differences The minimum number of different nt.
#' @return             The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setMinDifferences(exp, 10)
#'
#' @export
setGeneric( name="setMinDifferences",
            def =function(object, differences) {
                standardGeneric("setMinDifferences")
            }
)

#' @rdname  setMinDifferences-method
#' @export
setMethod(  f        ="setMinDifferences",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, differences) {
                                object@minDifferences <- differences
                return(object)
            }
)


#' Set transition probabilities (for the HMM step).
#' @rdname setTransitionProbabilities-method
#' @param  object       An \code{srnadiff} object.
#' @param  noDiffToDiff probability to change from the "not-differentially
#'                        expressed" state to the "differentially expressed"
#'                        state
#' @param  diffToNoDiff probability to change from the "differentially
#'                        expressed" state to the "not-differentially expressed"
#'                        state
#' @return         The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setTransitionProbabilities(exp, 0.001, 0.000001)
#'
#' @export
setGeneric( name="setTransitionProbabilities",
            def =function(object, noDiffToDiff, diffToNoDiff) {
                standardGeneric("setTransitionProbabilities")
            }
)

#' @rdname  setTransitionProbabilities-method
#' @export
setMethod(  f        ="setTransitionProbabilities",
            signature= c("sRNADiff", "numeric", "numeric"),
            definition=function(object, noDiffToDiff, diffToNoDiff) {
                                object@noDiffToDiff <- noDiffToDiff
                                object@diffToNoDiff <- diffToNoDiff
                return(object)
            }
)


#' Set emission probabilities (for the HMM step): probability to have a p-value
#'    not less than a threshold in the "not-differentially expressed" state,
#'    and a p-value not greater than this threshold in the "differentially
#'    expressed" state (supposed equal).
#' @rdname setEmissionProbabilities-method
#' @param  object       An \code{srnadiff} object.
#' @param  probability  The emission probability
#' @return              The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setEmissionProbabilities(exp, 0.9)
#'
#' @export
setGeneric( name="setEmissionProbabilities",
            def =function(object, probability) {
                standardGeneric("setEmissionProbabilities")
            }
)

#' @rdname  setEmissionProbabilities-method
#' @export
setMethod(  f        ="setEmissionProbabilities",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, probability) {
                                object@emission <- probability
                return(object)
            }
)


#' Set emission threshold (for the HMM step): the emission distribution being
#'    binomial, all the p-values less than this threshold belong to one class,
#'    and all the p-values greater than this threshold belong to the
#'    other class.
#' @rdname setEmissionThreshold-method
#' @param  object    An \code{srnadiff} object.
#' @param  threshold The emission threshold
#' @return           The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setEmissionThreshold(exp, 0.1)
#'
#' @export
setGeneric( name="setEmissionThreshold",
            def =function(object, threshold) {
                standardGeneric("setEmissionThreshold")
            }
)


#' @rdname  setEmissionThreshold-method
#' @export
setMethod(  f        ="setEmissionThreshold",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, threshold) {
                                object@emissionThreshold <- threshold
                return(object)
            }
)


#' Set minimum overlap (for the last quantification step): all the reads with
#'    at least n nucleotides shared with a feature will be used for
#'    quantification of this feature.
#' @rdname setMinOverlap-method
#' @param  object     An \code{srnadiff} object.
#' @param  minOverlap The minimum overlap
#' @return            The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setMinOverlap(exp, 10)
#'
#' @export
setGeneric( name="setMinOverlap",
            def =function(object, minOverlap) {
                standardGeneric("setMinOverlap")
            }
)

#' @rdname  setMinOverlap-method
#' @export
setMethod(  f        ="setMinOverlap",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, minOverlap) {
                                object@minOverlap <- minOverlap
                return(object)
            }
)


#' Set number of threads to use
#' @rdname setNThreads-method
#' @param  object   An \code{srnadiff} object.
#' @param  nThreads The number of threads
#' @return          The same object
#'
#' @examples
#' exp <- sRNADiffExample()
#' exp <- setNThreads(exp, 4)
#'
#' @export
setGeneric( name="setNThreads",
            def =function(object, nThreads) {
                standardGeneric("setNThreads")
            }
)

#' @rdname  setNThreads-method
#' @export
setMethod(  f        ="setNThreads",
            signature= c("sRNADiff", "numeric"),
            definition=function(object, nThreads) {
                object@nThreads <- nThreads
                if (nThreads > 1) register(MulticoreParam(workers = nThreads))
                return(object)
            }
)


#' Remove redundant regions
#'
#' @param regions An \code{GRanges} object.
#' @param padj A list of adjusted p-values.
#' @return A \code{GRanges} object.
removeRedundant <- function(regions, padj) {
    names(regions) <- paste0(seqnames(regions), "_", start(regions),
                                "_", end(regions))
    sizes          <- width(regions)
    overlaps       <- findOverlaps(regions, regions)
    overlaps       <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
    if (length(overlaps) == 0) {
        return(regions)
    }
    from      <- queryHits(overlaps)
    to        <- subjectHits(overlaps)
    dominance <- padj[from] < padj[to] |
                            (padj[from] == padj[to] & sizes[from] < sizes[to]) |
                            (padj[from] == padj[to] & sizes[from] == sizes[to]
                                & from < to)
    toBeRemoved <- c()
    while (TRUE) {
        dominated   <- table(to[dominance])
        dominator   <- table(from[dominance])
        both        <- intersect(names(dominated), names(dominator))
        if (length(both) == 0) break
        toBeRemoved <- union(toBeRemoved,
                                intersect(both, names(dominated[dominated ==
                                        min(dominated[both])])))
        dominance[(queryHits(overlaps) %in% toBeRemoved |
                    subjectHits(overlaps) %in% toBeRemoved)] <- FALSE
    }
    toBeRemoved    <- union(as.numeric(toBeRemoved), to[dominance])
    regions        <- regions[- toBeRemoved]
    return(regions)
}


#' Keep regions with best p-values.
#'
#' @param object An \code{srnadiff} object.
#' @param allSets A \code{GRanges} object.
#' @return A \code{GRanges} object.
reconcileRegions <- function(object, allSets) {
    message("Computing differential expression...")
    counts           <- summarizeOverlaps(allSets, object@bamFiles,
                            inter.feature=FALSE,
                            mode=function(features, reads, ignore.strand,
                                inter.feature)
                                countOverlaps(features, reads,
                                    minoverlap=object@minOverlap))
    colData(counts)  <- object@design
    dds              <- DESeqDataSet(counts, design =~condition)
    names(dds)       <- names(allSets)
    dds              <- dds[rowSums(counts(dds)) > 1, ]
    dds              <- DESeq(dds)
    regions          <- allSets[names(allSets) %in% names(dds),]
    mcols(regions)   <- results(dds)
    padj             <- mcols(regions)$padj
    regions          <- regions[! is.na(padj)]
    dds              <- dds[! is.na(padj)]
    padj             <- padj[! is.na(padj)]
    regions          <- removeRedundant(regions, padj)
    message("... done.")
    return(regions)
}


#' Compute normalization factors
#'
#' @param object An \code{srnadiff} object.
#' @return The normalization factor as a list of numeric.
computeNormalizationFactors <- function(object) {
    librarySize    <- vapply(lapply(object@coverages, sum), sum, integer(1))
    normalizationFactors <- rcpp_normalization(object@lengths, object@values,
                                               object@chromosomeSizes,
                                               librarySize)
    object@normalizationFactors <- normalizationFactors
    return(object)
}


#' Run the segmentation using 3 different methods, and reconcile them.
#'
#' @param object An \code{srnadiff} object.
#' @return A \code{GRanges} object.
#' @examples
#' exp        <- sRNADiffExample()
#' exp        <- runAll(exp)
#'
#' @export
runAll <- function(object) {
    if (object@nThreads > 1) {
        register(MulticoreParam(object@nThreads))
    }
    object         <- computeNormalizationFactors(object)
    setAnnotation  <- runAllAnnotation(object)
    setNaive       <- runAllNaive(object)
    setHmm         <- runAllHmm(object)
    setSlice       <- runAllSlice(object)
    allSets        <- unique(sort(do.call("c", list(setAnnotation, setNaive,
                                                    setHmm, setSlice))))
    object@regions <- reconcileRegions(object, allSets)
    return(object)
}
