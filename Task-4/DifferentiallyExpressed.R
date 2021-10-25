# load libraries
if (!require("pacman")) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
pacman::p_load(pacman, rio, DESeq2, apeglm, ggplot2, pheatmap)


# load count matrix
data <- read.csv("Task-4/Raw_counts.csv", header = T, row.names = 1)
info <- read.table("Task-4/colData.txt", header = T, sep = '\t')

# DESeq
dds <-DESeqDataSetFromMatrix(data, info, ~species + organ + condition)

# remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep,]

ddsDE <- DESeq(dds)

# export normalized read counts
normCounts <- counts(ddsDE, normalized = T)
write.csv(normCounts, "Task-4/Normalized_counts.csv")

# results
res <- results(ddsDE, alpha = 0.01)
resOrdered <- res[order(res$padj),]

write.csv(resOrdered, "Task-4/DESeq_results.csv")
summary(res)

normalCounts <- read.csv("Task-4/Normalized_counts.csv", header = T, row.names = 1)
deSeqRes <- read.csv("Task-4/DESeq_results.csv", header = T, row.names = 1)
deGenes <- subset(deSeqRes, padj <= 0.01)
write.csv(deGenes, "Task-4/DEgenes.csv")

# heat map to show fold change
allSig <- merge(normalCounts, deGenes, by = 0)
sigCounts <- allSig[2:25] 
row.names(sigCounts) <- allSig$Row.names
pheatmap(log2(sigCounts + 1), scale = 'row', treeheight_row = 0, treeheight_col = 0)


# clear the environment
rm(list = ls()) 
p_unload(all)  # Remove all add-ons
detach("package:datasets", unload = TRUE)  # For base
dev.off()  # But only if there IS a plot
