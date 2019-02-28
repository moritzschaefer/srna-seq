log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

contrast <- list(snakemake@params[["contrast"]][["positive"]], snakemake@params[["contrast"]][["negative"]])

listValues <- c(1/length(contrast[[1]]), -1/length(contrast[[2]]))

res <- results(dds, contrast=contrast, listValues=listValues, parallel=parallel)
# shrink fold changes for lowly expressed genes
# res <- lfcShrink(dds, res=res) # TODO shrinkage is not possible because betaPrior is TRUE
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])
