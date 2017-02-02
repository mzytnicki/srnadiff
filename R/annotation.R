run.all.annotation <- function(object) {
    database         <- makeTxDbFromGFF(object@gtf.file.name, format="gtf", circ_seqs=character())
    annotation       <- exonsBy(database, by = "gene")
    counts           <- count.features(object, annotation)
    diff.data        <- run.deseq2(object, counts)
    print(diff.data)
    diff.gene.names  <- rownames(diff.data)
    diff.genes       <- genes(database, filter=list(gene_id=diff.gene.names))
    diff.genes
}
