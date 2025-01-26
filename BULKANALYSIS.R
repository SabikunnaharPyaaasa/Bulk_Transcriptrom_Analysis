install.packages("ggplot2")
install.packages("pheatmap")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
# Install org.Hs.eg.db via Bioconductor
BiocManager::install("org.Hs.eg.db")
BiocManager::install("sva")

# 1. Load Required Libraries
library(edgeR)          # For normalization and DEG analysis
library(ggplot2)        # For data visualization
library(pheatmap)       # For heatmaps
library(clusterProfiler) # For GO enrichment
library(org.Hs.eg.db)   # Human gene annotations (adjust for other organisms)
library(limma)          # For linear model analysis
library(sva)            # For batch correction
update.packages(ask = TRUE)

# 2. Load Data
# Import raw counts and metadata

setwd("C:/Users/pyaas/OneDrive/Desktop/R")

data <- read.table("data.txt", header = TRUE, row.names = 1)

annotation <- read.table("annotation.txt", header = TRUE)

# Check structure of the data
head(data)
head(annotation)

# Check Data
print("Counts Dimensions:")
print(dim(counts))
print("Metadata Dimensions:")
print(dim(metadata))

# Quality Control (QC)
# Boxplot for raw counts
boxplot(log2(data + 1), las = 2, main = "Raw Count Distribution")

# PCA for raw counts
pca_raw <- prcomp(t(log2(data + 1)))
plot(pca_raw$x[,1:2], col=as.factor(annotation$Condition), pch=16)

# Normalization using TMM
dge <- DGEList(counts = data, group = annotation$Condition)
dge <- calcNormFactors(dge)
norm_counts <- cpm(dge, log=TRUE)

# Boxplot of normalized data
boxplot(norm_counts, las = 2, main = "Normalized Count Distribution")

# PCA for normalized data
pca_norm <- prcomp(t(norm_counts))
plot(pca_norm$x[,1:2], col=as.factor(annotation$Condition), pch=16)

# Differential Expression Analysis
design <- model.matrix(~annotation$Condition)
fit <- lmFit(norm_counts, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="fdr", number=Inf)
upregulated <- subset(results, logFC > 1 & adj.P.Val < 0.05)
downregulated <- subset(results, logFC < -1 & adj.P.Val < 0.05)

# Visualization (Heatmap of significant genes)
library(pheatmap)
selected_genes <- rownames(results)[results$adj.P.Val < 0.05]
pheatmap(norm_counts[selected_genes, ], cluster_rows=TRUE, cluster_cols=TRUE)


# GO Enrichment Analysis
upregulated_genes <- rownames(upregulated)
downregulated_genes <- rownames(downregulated)
go_up <- enrichGO(gene = upregulated_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(go_up, showCategory=20)

go_down <- enrichGO(gene = downregulated_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(go_down, showCategory=20)


# Bonus 1: Linear Model Analysis
fit_lm <- lmFit(norm_counts, design)
fit_lm <- eBayes(fit_lm)
results_lm <- topTable(fit_lm, coef=2, adjust="fdr", number=Inf)
pheatmap(norm_counts[rownames(results_lm)[results_lm$adj.P.Val < 0.05], ], cluster_rows=TRUE, cluster_cols=TRUE)

# Bonus 2: Batch Correction
if (!"Batch" %in% colnames(annotation)) {
  annotation$Batch <- factor(rep(1, nrow(annotation)))
  message("Batch column was missing and has been added with a single batch level.")
}

if (length(unique(annotation$Batch[!is.na(annotation$Batch)])) <= 1) {
  combat_data <- norm_counts
  message("Batch correction skipped as there is only one batch or batch data is missing.")
} else {
  annotation$Batch <- as.factor(annotation$Batch)
  mod <- model.matrix(~annotation$Condition)
  combat_data <- ComBat(dat=as.matrix(norm_counts), batch=annotation$Batch, mod=mod)
}

# Re-evaluate differential expression after batch correction
fit_batch <- lmFit(combat_data, design)
fit_batch <- eBayes(fit_batch)
results_batch <- topTable(fit_batch, coef=2, adjust="fdr", number=Inf)
