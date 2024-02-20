#### RNA-SEQ mNGS Tohoku-SC2 data ####
#### SMRojas ####
#### Wilcoxon rank-sum test ####


# Libraries ---------------------------------------------------------------


library(DESeq2)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ashr)
library(edgeR)


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



# Filtering ---------------------------------------------------------------

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


# Downstream ----------------------------------------------------

conditions <- factor(t(dds$result))

y <- DGEList(counts = counts(dds), group = conditions)

# y$counts <- y$counts + 0.5

#### Perform TMM normalization and convert to CPM (Counts Per Million)

y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)

#### Run the Wilcoxon rank-sum test for each gene

pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")

#### Calculate the fold-change for each gene

conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))

#### Output results based on the FDR threshold 0.05

outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
resLFC <- outRst[outRst$FDR<fdrThres,]
resLFC

# with background
topgenes_upregulated <- rownames(resLFC[resLFC$log2foldChange > 0.5,])
topgenes_upregulated
GO_results <- enrichGO(gene = topgenes_upregulated, OrgDb = "org.Hs.eg.db", universe = rownames(y), keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

# without background

GO_results <- enrichGO(gene = topgenes_upregulated, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

# with background
topgenes_downregulated <- rownames(resLFC[resLFC$log2foldChange < -0.5,])
GO_results <- enrichGO(gene = topgenes_downregulated, OrgDb = "org.Hs.eg.db", universe = rownames(y), keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

#without background

GO_results <- enrichGO(gene = topgenes_downregulated, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
KEGGmap <- plot(dotplot(GO_results, showCategory = 20))

# list genes upregulated
write.table(topgenes_upregulated, file = "genes_upregulated.tsv", sep = "\t", row.names = FALSE)

# other plots

res.df <- as.data.frame(resLFC)
res.df <- res.df[is.finite(res.df$log2foldChange), ]
res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = "ENSEMBL", column = "SYMBOL")

EnhancedVolcano(res.df, x = "log2foldChange", y = "pValues", lab = res.df$symbol)

