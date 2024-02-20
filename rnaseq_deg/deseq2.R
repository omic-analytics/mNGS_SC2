#### RNA-SEQ mNGS Tohoku-SC2 data ####
#### SMRojas ####
#### DESeq2 Method ####



# Libraries -------------------------------------------------------------

library(DESeq2)
library(readxl)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ashr)

# Imports -----------------------------------------------------------------

countdata <- read.table("counts.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
countdata <- as.matrix(countdata)
head(countdata)

metadata <- "metadata.xlsx"
metadata <- read_excel(metadata)

sampleNames <- colnames(countdata)
sub_metadata <- metadata[metadata$lab_id %in% sampleNames, ]
sub_metadata <- sub_metadata[match(sampleNames, sub_metadata$lab_id), ]
table(sub_metadata$result)



# DESeq2 ------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sub_metadata,
                              design = ~ result)
# counts(dds)+1
#### FILTERING LOW DEPTH SAMPLES

sort(colSums(counts(dds)))
keep_samples <- colSums(counts(dds)) >= 9000
dds <- dds[, keep_samples]
table(dds$result)

#### FILTERING GENES WITH LOW COUNTS

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#### SPECIFY CONDITION & Run DESeq

dds$result <- factor(dds$result, levels = c("positive","negative"))
dds$result <- relevel(dds$result, ref = "negative")
dds <- DESeq(dds)
res <- results(dds)
summary(res)

resLFC <- lfcShrink(dds, contrast=c("result","positive","negative"), type="ashr")
resLFC
summary(resLFC)


# Plots -------------------------------------------------------------------

plotDispEsts(dds, main = "Dispersion Estimates")

expmatrix_DESeq <- DESeq2::vst(dds, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)
pca_data <- plotPCA(expmatrix_DESeq, intgroup = "result") + stat_ellipse()
pca_data

#sample label
pca_coords <- pca_data$data
pca_data +
geom_text(aes(x = PC1, y = PC2, label = colnames(expmatrix_DESeq)), 
            hjust = 1.2, size = 3)

# with background
topgenes_upregulated <- rownames(resLFC[resLFC$log2FoldChange > 0.5 & resLFC$pvalue < 0.00000001,])
GO_results <- enrichGO(gene = topgenes_upregulated, OrgDb = "org.Hs.eg.db", universe = rownames(resLFC), keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

# without background

GO_results <- enrichGO(gene = topgenes_upregulated, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

# with background
topgenes_downregulated <- rownames(resLFC[resLFC$log2FoldChange < -0.5 & resLFC$pvalue < 0.05,])
GO_results <- enrichGO(gene = topgenes_downregulated, OrgDb = "org.Hs.eg.db", universe = rownames(resLFC), keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

#without background

GO_results <- enrichGO(gene = topgenes_downregulated, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))


# other plots

res.df <- as.data.frame(resLFC)
res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = "ENSEMBL", column = "SYMBOL")


EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol)


DESeq2::plotMA(resLFC, 
               alpha = 0.05,
               main = "positive vs negative",
               xlab = "mean of normalized counts",
               ylab = "log fold change",
               ylim = c(-5,5))
