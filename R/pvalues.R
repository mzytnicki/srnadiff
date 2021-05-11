##- Compute p-values of the selected counts ----------------------------------#
##----------------------------------------------------------------------------#
computePvalues <- function(object, counts, nThreads=1) {
    if (object@diffMethod == "deseq2") {
        return(useDESeq2(object, counts, nThreads))
    }
    if (object@diffMethod == "edger") {
        return(useEdgeR(object, counts, nThreads))
    }
    if (object@diffMethod == "bayseq") {
        return(useBaySeq(object, counts, nThreads))
    }
    stop("Cannot understand differential analysis method ", object@diffMethod,
         ".")
}


##- Use DESeq2 to compute p-values -------------------------------------------#
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
    return(padj)
}


##- Use edgeR to compute p-values --------------------------------------------#
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
    return(padj)
}


##- Use edgeR to compute p-values --------------------------------------------#
##----------------------------------------------------------------------------#
useBaySeq <- function(object, counts, nThreads = 1) {
    if (! is.factor(sampleInfo(object)$Condition)) {
        stop("Error!  Conditions should be factors.")
    }
    conditions <- factor(sampleInfo(object)$Condition)
    CD <- new("countData",
              data = counts,
              replicates = conditions,
              groups = list(NDE = rep.int(1, length(conditions)), 
                            DE = as.numeric(conditions)))
    libsizes(CD) <- getLibsizes(CD)
    CD <- suppressMessages(getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = NULL, verbose = FALSE))
    CD <- suppressMessages(getLikelihoods(CD, cl = NULL, bootStraps = 3, verbose = FALSE))
    res <- topCounts(CD, group = "DE", number = nrow(counts))
    return(res$FDR.DE)
}

