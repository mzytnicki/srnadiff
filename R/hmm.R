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

compute.observations <- function(object, map, chromosome) {
    values                <- vector()
    pos                   <- vector()
    chunk.size            <- 100000
    values.per.chromosome <- vector()
    pos.per.chromosome    <- vector()
    chromosome.size       <- object@chromosome.sizes[chromosome]
    for (i in 1:floor(chromosome.size/chunk.size)) {
        start <- chunk.size * i + 1
        end   <- min(chunk.size * (i+1), chromosome.size)
        counts.per.chunk <- Reduce(cbind, lapply(lapply(lapply(object@coverages, `[[`, chromosome), window, start = start, end = end), as.numeric))
        colnames(counts.per.chunk) <- NULL
        counts.per.chunk.as.vector <- apply(counts.per.chunk, 1, paste, collapse = " ")
        values.per.chunk <- unlist(mget(counts.per.chunk.as.vector, envir = map, ifnotfound = 1))
        values.per.chunk[is.na(values.per.chunk)] <- 1
        pos.per.chunk <- which(values.per.chunk < 1) + start - 1
        trimmed.values.per.chunk <- values.per.chunk[values < 1]
        values.per.chromosome <- c(values.per.chromosome, trimmed.values.per.chunk)
        pos.per.chromosome    <- c(pos.per.chromosome, pos.per.chunk)
    }
    return(list(pos = pos.per.chromosome, values = values.per.chromosome))
}

NO.DIFF.CLASS <- 1
DIFF.CLASS <- 2
N.CLASSES <- 2
NO.DIFF.THRESHOLD <- 0.85
DIFF.THRESHOLD <- 0.15
NO.DIFF.TO.DIFF <- 0.000001
DIFF.TO.NO.DIFF <- 0.001

run.hmm <- function(object, counts, pvalues, lengths, values) {
    print("Computing observations")
    transitions <- t(matrix(data = c(1 - NO.DIFF.TO.DIFF, NO.DIFF.TO.DIFF, DIFF.TO.NO.DIFF, 1 - DIFF.TO.NO.DIFF), ncol = 2))
    transitions <- -log10(transitions / rowSums(transitions))
    emissions   <- t(matrix(data = c(1, 9, 9, 1), ncol = 2))
    emissions   <- -log10(emissions / rowSums(emissions))
    starts      <- c(DIFF.TO.NO.DIFF / NO.DIFF.TO.DIFF, NO.DIFF.TO.DIFF / DIFF.TO.NO.DIFF)
    starts      <- -log10(starts / sum(starts))
    ranges      <- rcpp_viterbi(object@chromosome.sizes, transitions, emissions, starts, counts, pvalues, lengths, values)
    return(Granges(ranges))
}

run.all.hmm <- function(object) {
    lengths      <- lapply(object@chromosomes, function(chromosome) lapply(lapply(object@coverages, `[[`, chromosome), slot, "lengths"))
    values       <- lapply(object@chromosomes, function(chromosome) lapply(lapply(object@coverages, `[[`, chromosome), slot, "values"))
    counts       <- build.data.hmm(object, lengths, values)
    pvalues      <- compute.pvalues(object, counts)
    intervals    <- run.hmm(object, counts, pvalues, lengths, values)
    print(intervals)
    return(intervals)
}
