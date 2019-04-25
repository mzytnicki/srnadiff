##- Creates a matrix of unique counts ----------------------------------------#
##----------------------------------------------------------------------------#
buildDataHmm <- function (object) {
    counts <- rcpp_buildHmm(object@coverages, parameters(object)$minDepth)
    colnames(counts) <- sampleInfo(object)$SampleName
    rownames(counts) <- paste("data", seq(dim(counts)[1]), sep="_")

    return(counts)
}


##- Compute p-values of the selected counts ----------------------------------#
##----------------------------------------------------------------------------#
computePvalues <- function(object, counts, nThreads) {
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                    colData=sampleInfo(object),
                                    design=~Condition)
    sizeFactors(dds) <- normFactors(object)
    dds <- suppressMessages(DESeq(dds, parallel=(nThreads > 1)))
    res <- results(dds)
    pvalues <- res$pvalue

    return(pvalues)
}


##- Initialize and run the HMM -----------------------------------------------#
##----------------------------------------------------------------------------#
hmm <- function(object, counts, pvalues) {
    param <- parameters(object)

    transitions <- t(matrix(data=c(1 - param$noDiffToDiff,
                            param$noDiffToDiff, param$diffToNoDiff,
                            1 - param$diffToNoDiff), ncol=2))
    transitions <- -log10(transitions / rowSums(transitions))

    emissions <- t(matrix(data=c(1 - param$emission, param$emission,
                                param$emission, 1 - param$emission), ncol=2))
    emissions <- -log10(emissions / rowSums(emissions))

    starts <- c(param$noDiffToDiff / param$diffToNoDiff,
                param$diffToNoDiff / param$noDiffToDiff)
    starts <- -log10(starts / sum(starts))

    ranges <- rcpp_viterbi(coverages(object), transitions,
                            emissions, param$emissionThreshold,
                            starts, counts, pvalues, param$minDepth,
                            param$minSize, param$maxSize)

    if (nrow(ranges) == 0) return(GRanges())

    return(GRanges(ranges))
}

##- Segmentation of the genome using the HMM method --------------------------#
##----------------------------------------------------------------------------#
runHmm <- function(object, nThreads) {

    message("Starting HMM step...")
    message("  Building data...")
    counts <- buildDataHmm(object)
    message("  ... data built")

    message("  Computing p-values...")
    pvalues <- computePvalues(object, counts, nThreads)
    message("  ... values computed")

    message("  Running HMM...")
    intervals <- hmm(object, counts, pvalues)
    message("  ... HMM run.")

    if (length(intervals) > 0) {
        names(intervals) <- paste("hmm", seq(length(intervals)), sep="_")
    }

    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... HMM step done.\n")

    return(intervals)
}
