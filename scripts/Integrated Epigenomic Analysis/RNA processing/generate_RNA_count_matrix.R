rm(list=ls())

library(GenomicFeatures) # version 1.38.0
library(tximport) # version 1.14.2
library(DESeq2) # version 1.26.0
# R version 3.6.3

# Create TxDb object from GENCODE annotation file (PMACS) #

# TxDb objects stores transcript metadata (5' and 3' UTRs, protein-coding sequences, exons for mRNA transcripts)
txdb.vM24 <- makeTxDbFromGFF(file="/home/jingyaq/Minn/references/GRCm38/gencode.vM24.annotation.gtf", format="gtf", dataSource="GENCODE vM24", organism="Mus musculus")
saveDb(txdb.vM24, file="/home/jingyaq/Minn/references/GRCm38/txdb.GENCODE_vM24")

##########################################
### Import salmon counts with tximport ###
##########################################

# Do this on local computer with R >= 3.6, tximport >= 1.14 (quantification with salmon version 1.2.0)

dir.path <- "~/Dropbox/Minn/ifnar_epigenome/data/RNA/"

sampleIDs <- list.dirs(paste0(dir.path, "salmon"), full.names=F, recursive=F)
conditions <- c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))

coldata <- data.frame(sample=sampleIDs, condition=conditions)
rownames(coldata) <- sampleIDs

# Associate transcripts with gene IDs for gene-level summarization
txdb <- loadDb("~/Dropbox/Minn/references/GRCm38/txdb.GENCODE_vM24")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys=k, columns=c("GENEID"), keytype="TXNAME")

# Import transcript-level estimates
files <- paste0(dir.path, coldata$run, "salmon/", coldata$sample, "/quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene, varReduce=T)

####################################
### rlog transform and normalize ###
####################################

# Make DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)

# Filter out lowly expressed genes
keep <- rowSums(counts(dds)) >= 3
dds <- dds[keep, ]

# Run DESeq2 to estimate size factors
dds <- DESeq(dds)
rownames(dds) <- sapply(strsplit(rownames(dds), split="\\."), function(x) x[[1]]) # Remove gene version information

# rlog transform
rld <- rlog(dds, blind=FALSE)

# Write out rlog counts
edat <- assay(rld)
write.table(edat, file=paste0(dir.path, "mat/RNA_counts_rlog.txt"), sep="\t", quote=F, col.names=T, row.names=T)

# Save raw counts
write.table(assay(dds), file=paste0(dir.path, "mat/RNA_counts_raw.txt"), sep="\t", quote=F, col.names=T, row.names=T)

# # Save normalized data (divide raw counts by size factors)
# write.table(counts(dds, normalized=T), file=paste0(dir.path, "mat/RNA_counts_normalized.txt"), sep="\t", quote=F, col.names=T, row.names=T)

# Save dds object
saveRDS(dds, file=paste0(dir.path, "DE/RNA_counts_dds.rds"))

res.B16_v_R499 <- results(dds, contrast=c("condition", "R499_cas", "B16_cas"))
res.B16_v_B16_SKO <- results(dds, contrast=c("condition", "B16_SKO", "B16_cas"))
res.R499_v_R499_SKO <- results(dds, contrast=c("condition", "R499_SKO", "R499_cas"))
res.B16_SKO_v_R499_SKO <- results(dds, contrast=c("condition", "R499_SKO", "B16_SKO"))

DEGs <- list(res.B16_v_R499, res.B16_v_B16_SKO, res.R499_v_R499_SKO, res.B16_SKO_v_R499_SKO)
names(DEGs) <- c("B16_v_R499", "B16_v_B16_SKO", "R499_v_R499_SKO", "B16_SKO_v_R499_SKO")
saveRDS(DEGs, file=paste0(dir.path, "DE/DEG_list.rds"))




