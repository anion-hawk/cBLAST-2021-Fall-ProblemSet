# load libraries
if (!require("pacman")) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
pacman::p_load(pacman, rio, DESeq2, apeglm, ggplot2, pheatmap)


DiffExp <- function(countdata, colinfo){
  dds <-DESeqDataSetFromMatrix(countdata, colinfo, ~condition)
  
  # remove lowly expressed genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  return(DESeq(dds))
}

# load count matrix
data <- read.csv("Raw_counts.csv", header = T, row.names = 1)
info <- read.table("colData.txt", header = T, sep = '\t')

# Differential Expression Objects
ddsDE_VL <- DiffExp(data[, 1:6], info[1:6,])
ddsDE_VR <- DiffExp(data[, 7:12], info[7:12,])
ddsDE_BL <- DiffExp(data[, 13:18], info[13:18,])
ddsDE_BR <- DiffExp(data[, 19:24], info[19:24,])


# export normalized read counts
normCounts_VL <- counts(ddsDE_VL, normalized = T)
normCounts_VR <- counts(ddsDE_VR, normalized = T)
normCounts_BL <- counts(ddsDE_BL, normalized = T)
normCounts_BR <- counts(ddsDE_BR, normalized = T)

normCounts_V <- merge(normCounts_VL, normCounts_VR, by=0)
normCounts_B <- merge(normCounts_BL, normCounts_BR, by=0)
normCounts <- merge(normCounts_V[2:13], normCounts_B[2:13], by=0)
normCounts <- normCounts[2:25]
row.names(normCounts) <- normCounts_V$Row.names

write.csv(normCounts, "Normalized_counts.csv")


# results
res_VL <- results(ddsDE_VL, alpha = 0.01)
res_VR <- results(ddsDE_VR, alpha = 0.01)
res_BL <- results(ddsDE_BL, alpha = 0.01)
res_BR <- results(ddsDE_BR, alpha = 0.01)

resOrdered_VL <- res_VL[order(res_VL$padj),]
resOrdered_VR <- res_VR[order(res_VR$padj),]
resOrdered_BL <- res_BL[order(res_BL$padj),]
resOrdered_BR <- res_BR[order(res_BR$padj),]

write.csv(resOrdered_VL, "DESeq_results_VL.csv")
write.csv(resOrdered_VR, "DESeq_results_VR.csv")
write.csv(resOrdered_BL, "DESeq_results_BL.csv")
write.csv(resOrdered_BR, "DESeq_results_BR.csv")

summary(res_VL)
summary(res_VR)
summary(res_BL)
summary(res_BR)


# extract differentially expressed genes
deSeqRes_VL <- read.csv("DESeq_results_VL.csv", header = T, row.names = 1)
deSeqRes_VR <- read.csv("DESeq_results_VR.csv", header = T, row.names = 1)
deSeqRes_BL <- read.csv("DESeq_results_BL.csv", header = T, row.names = 1)
deSeqRes_BR <- read.csv("DESeq_results_BR.csv", header = T, row.names = 1)
deGenes_VL <- subset(deSeqRes_VL, padj <= 0.01)
deGenes_VR <- subset(deSeqRes_VR, padj <= 0.01)
deGenes_BL <- subset(deSeqRes_BL, padj <= 0.01)
deGenes_BR <- subset(deSeqRes_BR, padj <= 0.01)

write.csv(deGenes_VL, "DEgenes_VL.csv")
write.csv(deGenes_VR, "DEgenes_VR.csv")
write.csv(deGenes_BL, "DEgenes_BL.csv")
write.csv(deGenes_BR, "DEgenes_BR.csv")


# heatmap of normalized counts
normCounts <- read.csv("Normalized_counts.csv", header = T, row.names = 1)
pheatmap(log2(normCounts+1), scale = 'row', treeheight_row = 0, treeheight_col = 0)


# clear the environment
rm(list = ls()) 
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base
dev.off()  # But only if there IS a plot

