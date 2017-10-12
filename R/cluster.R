#' Compute normalization factors using RLE values.
#'
#' @param object An \code{srnadiff} object.
#' @param counts the merged coverage
#' @param libSize the library sizes
#' @return The normalization factors
calcFactorRLE <- function (object, counts, libSize) {
    #-- Remove all zero rows
    nonzero <- (Reduce('+', counts) > 0)
    counts <- sapply(counts, function (x) { subset(x, nonzero) })
    #-- Calculate factors
    gm <- exp(Reduce('+', sapply(counts, log)) / length(counts))
    f <- sapply(counts, function(x) { median((x / gm)[gm > 0]) }) / libSize
    #-- Factors should multiple to one
    f <- f / exp(mean(log(f)))
    #-- Output
    return(f)
}

#' Normalize counts
#'
#' @param object An \code{srnadiff} object.
#' @param mergedCounts the counts: one per condition
#' @param libSize the library size
#' @param factors the normalization factors
#' @return The normalized counts
CPM <- function(object, mergedCounts, libSize, factors) {
    n       <- length(libSize)
    counts  <- object@coverages
    libSize <- libSize * factors
    md      <- median(sapply(mergedCounts, mean))
    counts  <- mapply('/', counts, libSize)
    md.norm <- median(sapply(sapply(counts, unlist), mean))
    fs      <- md / md.norm
    counts  <- lapply(counts, function(x) { round(x * fs) } )
    #-- Output
    return(counts)
}

#' Initialize and run the clustering.
#'
#' @param object An \code{srnadiff} object.
#' @param counts The normalized counts
#' @return A GRanges.
runClustering <- function(object, counts) {
    MIN_DEPTH      <- 10
    MIN_SIZE       <- 10
    MAX_SIZE       <- 10000
    MIN_DIFFERENCE <- 20
    avgCounts      <- sapply(split(counts,
                                object@conditions == object@conditions[[1]]),
                            function (s) { round(Reduce('+', s) / length(s)) })
    lengths <- lapply(avgCounts, function (x) {lapply(x, slot, "lengths") })
    values  <- lapply(avgCounts, function (x) {lapply(x, slot, "values") })
    ranges  <- rcpp_clustering(lengths, values, object@chromosomeSizes,
                                MIN_DEPTH, MIN_SIZE, MAX_SIZE, MIN_DIFFERENCE)
    return(GRanges(ranges))
}


#' Segmentation of the genome using a clustering method.
#'
#' @param object An \code{srnadiff} object.
#' @return       A \code{GRanges} object.
runAllClustering <- function(object) {
    if (object@skipClustering) {
        return(GRanges())
    }
    message("Starting clustering step...")
    message("  Normalizing data...")
    mergedCounts <- sapply(object@coverages, unlist)
    libSize      <- sapply(mergedCounts, sum)
    factors      <- calcFactorRLE(object, mergedCounts, libSize)
    normCounts   <- CPM(object, mergedCounts, libSize, factors)
    message("  ... data normalized")
    message("  Clustering...")
    intervals <- runClustering(object, normCounts)
    message("  ... clustered.")
    names(intervals) <- paste("cluster", seq(length(intervals)), sep="_")
    message(paste0(c("  ... ", length(intervals), " regions found.")))
    message("... clustering step done.")
    return(intervals)
}
