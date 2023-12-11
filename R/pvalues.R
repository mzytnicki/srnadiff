##- Compute p-values and log2-fold changes of the selected counts ------------#
##----------------------------------------------------------------------------#
computePvalues <- function(object, counts, nThreads=1) {
    if (object@diffMethod == "deseq2") {
        return(useDESeq2(object, counts, nThreads))
    }
    if (object@diffMethod == "edger") {
        return(useEdgeR(object, counts, nThreads))
    }
    stop("Cannot understand differential analysis method ", object@diffMethod,
         ".")
}


##- Use DESeq2 to compute p-values and log2-fold changes ---------------------#
##----------------------------------------------------------------------------#
useDESeq2 <- function(object, counts, nThreads = 1) {
    colData <- sampleInfo(object)
    colData$Condition <- factor(colData$Condition)
    rse <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                colData = colData)
    dds <- DESeqDataSet(rse, design = ~Condition)
    sizeFactors(dds) <- normFactors(object)
    dds <- suppressMessages(DESeq(dds, parallel = (nThreads > 1)))
    if (! all(names(dds) == rownames(counts))) {
        stop("Error!  Region names have changed.")
    }
    padj <- results(dds)$padj
    padj <- ifelse(is.na(padj), 1.0, padj)
    log2FC <- results(dds)$log2FoldChange
    return(list(padj = padj, log2FC = log2FC))
}


##- Use edgeR to compute p-values and log2-fold changes ----------------------#
##----------------------------------------------------------------------------#
useEdgeR <- function(object, counts, nThreads = 1) {
    group <- factor(sampleInfo(object)$Condition)
    y <- DGEList(counts = counts, group = group)
    y <- calcNormFactors(y)
    y$samples$norm.factors <- 1 / normFactors(object)
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    padj <- p.adjust(qlf$table$PValue, method="BH")
    log2FC <- qlf$table$logFC
    return(list(padj = padj, log2FC = log2FC))
}
