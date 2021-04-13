log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
coldata <- read.table(snakemake@params[["units"]], header=TRUE, row.names="sample", sep=",")
names(coldata)[names(coldata) == "unit"] <- "type"
coldata <- coldata[1]

target_counts <- read.delim(snakemake@input[["target_counts"]], row.names="Geneid")
# tDRs and mirbase_mmu21 have different number of index columns, so we have to cut the first columns
if (ncol(target_counts) > nrow(coldata)) {  # if there are more columns than samples
  index_cols = ncol(target_counts) - nrow(coldata)
  target_counts <- target_counts[-c(1:index_cols)]  # cut the index_cols
}

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

dds.normal <- DESeqDataSetFromMatrix(countData=target_counts,
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
