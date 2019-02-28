log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

mirna_counts <- read.table(snakemake@input[["mirna_counts"]], header=TRUE, row.names="Geneid", sep="\t")
mirna_counts <- mirna_counts[-c(1:5)]


normalization_counts <- rbind(
  read.table(snakemake@input[["trna_counts"]], header=TRUE, row.names="Geneid", sep="\t"),
  read.table(snakemake@input[["snrna_counts"]], header=TRUE, row.names="Geneid", sep="\t"),
  read.table(snakemake@input[["snorna_counts"]], header=TRUE, row.names="Geneid", sep="\t"),
  read.table(snakemake@input[["mt_trna_counts"]], header=TRUE, row.names="Geneid", sep="\t")
)

normalization_counts <- normalization_counts[-c(1:5)]

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", sep=",")
# cts <- cts[-c(1:5)]

dds.normal <- DESeqDataSetFromMatrix(countData=mirna_counts,
                              colData=coldata,
                              design=~ 0+type)

dds.special <- DESeqDataSetFromMatrix(countData=normalization_counts,
                              colData=coldata,
                              design=~ 0+type)

dds.normal <- estimateSizeFactors(dds.normal)
dds.special <- estimateSizeFactors(dds.special)

write.table(counts(dds.normal, normalized=TRUE), file=snakemake@output[["normal_normalized_counts"]])

sizeFactors(dds.normal) <- sizeFactors(dds.special)

write.table(counts(dds.normal, normalized=TRUE), file=snakemake@output[["special_normalized_counts"]])

# no need to remove uninformative columns -> commented
# dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds.normal, parallel=parallel)

saveRDS(dds, file=snakemake@output[["all"]])
