build.data.hmm <- function (object, lengths, values) {
    counts           <- rcpp_buildHmm(lengths, values, object@chromosome.sizes)
    counts           <- do.call(rbind, counts)
    colnames(counts) <- object@replicates
    rownames(counts) <- paste("data", seq(dim(counts)[1]), sep = "_")
    return(counts)
}

compute.pvalues <- function(object, counts) {
    dds     <- DESeqDataSetFromMatrix(countData = counts, colData = object@design, design = ~ condition)
    dds     <- DESeq(dds)
    res     <- results(dds)
    pvalues <- res$padj
    return(pvalues)
}

NO.DIFF.THRESHOLD <- 0.85
DIFF.THRESHOLD <- 0.15
NO.DIFF.TO.DIFF <- 0.000001
DIFF.TO.NO.DIFF <- 0.001

run.hmm <- function(object, counts, pvalues, lengths, values) {
    transitions <- t(matrix(data = c(1 - NO.DIFF.TO.DIFF, NO.DIFF.TO.DIFF, DIFF.TO.NO.DIFF, 1 - DIFF.TO.NO.DIFF), ncol = 2))
    transitions <- -log10(transitions / rowSums(transitions))
    emissions   <- t(matrix(data = c(1, 9, 9, 1), ncol = 2))
    emissions   <- -log10(emissions / rowSums(emissions))
    starts      <- c(DIFF.TO.NO.DIFF / NO.DIFF.TO.DIFF, NO.DIFF.TO.DIFF / DIFF.TO.NO.DIFF)
    starts      <- -log10(starts / sum(starts))
    ranges      <- rcpp_viterbi(object@chromosome.sizes, transitions, emissions, starts, counts, pvalues, lengths, values)
    return(GRanges(ranges))
}

run.all.hmm <- function(object) {
    lengths          <- lapply(object@chromosomes, function(chromosome) lapply(lapply(object@coverages, `[[`, chromosome), slot, "lengths"))
    values           <- lapply(object@chromosomes, function(chromosome) lapply(lapply(object@coverages, `[[`, chromosome), slot, "values"))
    counts           <- build.data.hmm(object, lengths, values)
    pvalues          <- compute.pvalues(object, counts)
    intervals        <- run.hmm(object, counts, pvalues, lengths, values)
    names(intervals) <- paste("hmm", seq(length(intervals)), sep = "_")
    return(intervals)
}
