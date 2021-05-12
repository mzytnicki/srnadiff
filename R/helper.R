##- Remove redundant regions -------------------------------------------------#
##----------------------------------------------------------------------------#
removeRedundant <- function(regions) {
    padj           <- mcols(regions)$padj
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


##- Statistic quantification of the DERs -------------------------------------#
##----------------------------------------------------------------------------#
reconcileRegions <- function(object, allRegions, minOverlap) {
    message("Quantifying differential expression...")

    countM <- summarizeOverlaps(allRegions, bamFiles(object),
                                inter.feature=FALSE,
                                mode=function(features, reads, ignore.strand,
                                                inter.feature) {
                                        countOverlaps(features, reads,
                                                minoverlap=minOverlap)
                                    }
                                )
    counts <- assays(countM)$counts
    padj <- computePvalues(object, counts)
    if (length(padj) != length(allRegions)) {
        stop("Error!  Cannot compute p-values: lengths differ.")
    }
    mcols(allRegions)$padj <- padj

    recRegions <- removeRedundant(allRegions)
    message("... done.")

    return(list(regions = recRegions, countMatrix = countM))
}


##- Compute normalization factors --------------------------------------------#
##----------------------------------------------------------------------------#
computeNormFactors <- function(cvg) {
    librarySize <- unlist(lapply(lapply(cvg, sum), sum))
    normFactors <- rcpp_normalization(cvg, librarySize)

    return(normFactors)
}


##- Coverage normalization ---------------------------------------------------#
##----------------------------------------------------------------------------#
cvgNormalization <- function(object) {
    librarySize <- unlist(lapply(lapply(coverages(object), sum), sum))
    md <- median(librarySize / sum(chromosomeSizes(object)))
    librarySize <- librarySize * normFactors(object)

    #-- normalization
    normCvg <- mapply('/', coverages(object), librarySize)

    normLibSize <- unlist(lapply(lapply(normCvg, sum), sum))
    normMd <- median(normLibSize / sum(chromosomeSizes(object)))
    scaleFactor <- md / normMd

    normCvg <- mapply(function(x, fs) { round(fs * x) },
                      normCvg, MoreArgs=list(fs = scaleFactor))

    normCvg <- lapply(normCvg, function(x, chrNames) {
                                    names(x) <- chrNames
                                    return(x) },
                      chrNames=names(chromosomeSizes(object)))
    return(normCvg)
}


##- Compute fold-change ------------------------------------------------------#
##----------------------------------------------------------------------------#
computeLogFoldChange <- function(object) {
    conditions <- as.character(sampleInfo(object)[["Condition"]])
    avgCounts <- lapply(split(cvgNormalization(object), conditions),
                        function(s) { round(Reduce('+', s) / length(s)) })

    lowValues <- (pmin(avgCounts[[1]], avgCounts[[2]]) <
                    parameters(object)$minDepth)
    avgCounts[[1]][lowValues] <- 0
    avgCounts[[2]][lowValues] <- 0
    logFC <- log2((avgCounts[[2]] + 1) / (avgCounts[[1]] + 1))

    return(logFC)
}


##- From IRanges list to GRanges ---------------------------------------------#
##----------------------------------------------------------------------------#
IRlist2GR <- function(IRlist) {
    res <- NULL

    for (i in seq_along(IRlist)) {
        df <- as.data.frame(IRlist[[i]])
        df$seqnames <- rep(names(IRlist)[i], nrow(df))
        res <- rbind(res, df)
    }

    res <- makeGRangesFromDataFrame(res)
    return(res)
}


##- Helper function for check input parameters -------------------------------#
##----------------------------------------------------------------------------#
checkParameters <- function(value) {

    minDepth <- value$minDepth
    minSize <- value$minSize
    maxSize <- value$maxSize
    minGap <- value$minGap
    maxDiff <- value$maxDiff
    minOverlap <- value$minOverlap
    noDiffToDiff <- value$noDiffToDiff
    diffToNoDiff <- value$diffToNoDiff
    emission <- value$emission
    emissionThreshold <- value$emissionThreshold
    cutoff <- value$cutoff
    minLogFC <- value$minLogFC

    ##- minDepth
    if (length(minDepth) != 1) {
        stop("'minDepth' must be of length 1.", call.=FALSE)
    }

    if (is.null(minDepth) || !is.numeric(minDepth) || (minDepth < 0) ||
        !is.finite(minDepth)) {
        stop("invalid value for 'minDepth'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    dec <- minDepth - trunc(minDepth)

    if (dec > 0) {
        stop("invalid value for 'minDepth'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- minSize
    if (length(minSize) != 1) {
        stop("'minSize' must be of length 1.", call.=FALSE)
    }

    if (is.null(minSize) || !is.numeric(minSize) || (minSize < 0) ||
        !is.finite(minSize)) {
        stop("invalid value for 'minSize'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    dec <- minSize - trunc(minSize)

    if (dec > 0) {
        stop("invalid value for 'minSize'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- maxSize
    if (length(maxSize) != 1) {
        stop("'maxSize' must be of length 1.", call.=FALSE)
    }

    if (is.null(maxSize) || !is.numeric(maxSize) || (maxSize < 0) ||
        !is.finite(maxSize)) {
        stop("invalid value for 'maxSize'. It must be a not negative",
             " integer.", call.=FALSE)
    }

    dec <- maxSize - trunc(maxSize)

    if (dec > 0) {
        stop("invalid value for 'maxSize'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- minGap
    if (length(minGap) != 1) {
        stop("'minGap' must be of length 1.", call.=FALSE)
    }

    if (is.null(minGap) || !is.numeric(minGap) ||
        (minGap < 0) || !is.finite(minGap)) {
        stop("invalid value for 'minGap'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    dec <- minGap - trunc(minGap)

    if (dec > 0) {
        stop("invalid value for 'minGap'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- maxDiff
    if (length(maxDiff) != 1) {
        stop("'maxDiff' must be of length 1.", call.=FALSE)
    }

    if (is.null(maxDiff) || !is.numeric(maxDiff) ||
        (maxDiff < 0) || !is.finite(maxDiff)) {
        stop("invalid value for 'maxDiff'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    dec <- maxDiff - trunc(maxDiff)

    if (dec > 0) {
        stop("invalid value for 'maxDiff'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- noDiffToDiff
    if (length(noDiffToDiff) != 1) {
        stop("'noDiffToDiff' must be a single value.", call.=FALSE)
    }

    if (is.null(noDiffToDiff) || !is.numeric(noDiffToDiff) ||
        !is.finite(noDiffToDiff)) {
        stop("'noDiffToDiff' value must be numeric.", call.=FALSE)
    }

    if ((noDiffToDiff > 1) || (noDiffToDiff < 0)) {
        stop("'noDiffToDiff' value outside the interval [0,1].", call.=FALSE)
    }

    ##- diffToNoDiff
    if (length(diffToNoDiff) != 1) {
        stop("'diffToNoDiff' must be a single value.", call.=FALSE)
    }

    if (is.null(diffToNoDiff) || !is.numeric(diffToNoDiff) ||
        !is.finite(diffToNoDiff)) {
        stop("'diffToNoDiff' value must be numeric.", call.=FALSE)
    }

    if ((diffToNoDiff > 1) || (diffToNoDiff < 0)) {
        stop("'diffToNoDiff' value outside the interval [0,1].", call.=FALSE)
    }

    ##- emission
    if (length(emission) != 1) {
        stop("'emission' must be a single value.", call.=FALSE)
    }

    if (is.null(emission) || !is.numeric(emission) ||
        !is.finite(emission)) {
        stop("'emission' value must be numeric.", call.=FALSE)
    }

    if ((emission > 1) || (emission < 0)) {
        stop("'emission' value outside the interval [0,1].", call.=FALSE)
    }

    ##- emissionThreshold
    if (length(emissionThreshold) != 1) {
        stop("'emissionThreshold' must be a single value.", call.=FALSE)
    }

    if (is.null(emissionThreshold) || !is.numeric(emissionThreshold) ||
        !is.finite(emissionThreshold)) {
        stop("'emissionThreshold' value must be numeric.", call.=FALSE)
    }

    if ((emissionThreshold > 1) || (emissionThreshold < 0)) {
        stop("'emissionThreshold' value outside the interval [0,1].",
             call.=FALSE)
    }

    ##- minOverlap
    if (length(minOverlap) != 1) {
        stop("'minOverlap' must be of length 1.", call.=FALSE)
    }

    if (is.null(minOverlap) || !is.numeric(minOverlap) || (minOverlap < 0) ||
        !is.finite(minOverlap)) {
        stop("invalid value for 'minOverlap'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    dec <- minOverlap - trunc(minOverlap)

    if (dec > 0) {
        stop("invalid value for 'minOverlap'. It must be a non-negative",
             " integer.", call.=FALSE)
    }

    ##- cutoff
    if (length(cutoff) != 1) {
        stop("'cutoff' must be of length 1.", call.=FALSE)
    }

    if (is.null(cutoff) || !is.numeric(cutoff) ||
        (cutoff < 0) || !is.finite(cutoff)) {
        stop("invalid value for 'cutoff'. It must be a non-negative",
             " number.", call.=FALSE)
    }

    ##- minLogFC
    if (length(minLogFC) != 1) {
        stop("'minLogFC' must be of length 1.", call.=FALSE)
    }

    if (is.null(minLogFC) || !is.numeric(minLogFC) ||
        (minLogFC < 0) || !is.finite(minLogFC)) {
        stop("invalid value for 'minLogFC'. It must be a non-negative",
             " number.", call.=FALSE)
    }
}
