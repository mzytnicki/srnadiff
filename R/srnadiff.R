run.all <- function(object) {
    set.annotation             <- run.all.annotation(object)
    set.naive                  <- run.all.naive(object)
    set.hmm                    <- run.all.hmm(object)
    all.sets                   <- c(set.annotation, set.naive, set.hmm)
    all.counts                 <- count.features(object, all.sets)
    all.pvalues                <- run.deseq2(object, all.counts)
    diff.sets                  <- all.sets[names(all.sets) %in% rownames(all.pvalues)]
    padj                       <- all.pvalues$padj
    overlaps                   <- findOverlaps(diff.sets, diff.sets)
    undominated                <- rep(TRUE, length(overlaps))
    undominated[overlaps@from] <- (padj[overlaps@from] <= padj[overlaps@to])
    undominated[overlaps@to]   <- undominated[overlaps@to] || (padj[overlaps@from] <= padj[overlaps@to])
    diff.sets                  <- diff.sets[undominated]
    print(diff.sets)
}
