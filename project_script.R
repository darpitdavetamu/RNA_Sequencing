######
###### RNA-Seq data on human liver tissue, from three females and three males. 
######
######   http://www.ncbi.nlm.nih.gov/pubmed?term=20009012
######

####
#### Load the data. There are separate files for the read counts and the phenotype data. 
#### The phenotype data tells us that there were actually two technical replicates of 
#### each biological sample (each person). Those have been pooled to create a single row 
#### of read counts in the count table, for each sample. The first three samples are 
#### females, and the last three samples are males.
####
rm(list=ls())
Y <- read.csv("merged_count.csv", header = TRUE)
GRP <- read.csv("GroupInfor.csv", header = TRUE,sep="\t")

Y <- data.frame(Y[,-c(1,2)],row.names = Y[,2])
GRP <- GRP[, -1]
colnames(Y) <- substring(colnames(Y), 2)
Y <- Y[, order(colnames(Y))]
GRP <- GRP[order(rownames(GRP)), ]



## There are many ENSEMBL genes that had no counts in any of the samples. Remove them 
## from the analysis.
n_count <- drop(data.matrix(Y) %*% rep(1, 27))
Y <- Y[n_count > 0, ]

####
#### Analysis with edgeR.
####

library(edgeR)

## Create a DGElist object for holding the counts along with the library sizes (the total 
## number of counts) and the group membership information for each sample. 
edgeR_dge <- DGEList(Y, lib.size = colSums(Y), group = GRP[])
edgeR_dge

## Estimate the "common dispersion parameter." 
edgeR_disp <- estimateCommonDisp(edgeR_dge)

## Differential expression analysis. Extract the top 100 genes. There are 63 genes with 
## an estimated FDR < 0.1.
edgeR_test <- exactTest(edgeR_disp)
edgeR_out <- topTags(edgeR_test, n = nrow(Y), sort.by = "PValue")
edgeR_out_sig <- subset(edgeR_out$table, FDR < 0.1)

head(edgeR_out_sig)

## Look at the read counts for the significant genes. Interestingly, nearly every one of 
## them appears significant based on a single large count in one of the comparison groups.
nms_edgeR_sig <- rownames(edgeR_out_sig)
Y[nms_edgeR_sig, ]

## Histogram of p-values. There is a large spike at 1. In fact, there are nearly 1,600 
## genes with p-value = 1. Some, but not all, of those genes have low overall read counts.
hist(edgeR_out$table$PValue)
table(edgeR_out$table$PValue == 1)

head(edgeR_dge$counts[edgeR_out$table$PValue == 1, ], 100)

####
#### Analysis with DESeq2.
####

library(DESeq2)

## Create a DESeqDataSet object for holding the counts and phenotype data.    
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = Y, colData = GRP, design = ~ Group + Sex)

## Apply rlog transformation.
deseq2_rld <- rlog(deseq2_deseqds)

## Cluster the samples. There appear to be three clusters. The third female clusters by 
## itself, and the first male clusters by itself. Referring to Table S3 from the 
## supplementary information of this publication, there are some aspects of these two 
## samples that are distinct from the others.
deseq2_hc <- hclust(dist(t(assay(deseq2_rld))))

par(mfrow = c(1, 2))
plot(deseq2_hc)
plot(deseq2_hc$height)

## PCA plot and analysis. Using DESeq2's plotPCA, it is clear that two samples are 
## different from the others. Running the PCA ourselves and looking at the loadings, we 
## can see, as expected, that these are the two samples mentioned above.
plotPCA(deseq2_rld, ntop = nrow(assay(deseq2_rld)), intgroup = "Age")

deseq2_pca <- prcomp(t(assay(deseq2_rld)))
deseq2_pca$x

## SVD analysis. There appear to be at least two interesting eigengenes. Respectively, 
## these seem to highlight the two samples from above.
deseq2_svd <- svd(t(scale(t(assay(deseq2_rld)), center = TRUE, scale = FALSE)))
round(deseq2_svd$d ^ 2 / sum(deseq2_svd$d ^ 2), 2)

par(mfrow = c(1, 1))
plot(deseq2_svd$v[, 1], type = "l", col = "blue", lwd = 2, ylim = c(-0.8, 0.8), 
     xlab = "Sample Index", ylab = "Eigengene Coefficient")
lines(deseq2_svd$v[, 2], col = "red", lwd = 2)
legend(1, 0.6, legend = c("36%", "27%"), lwd = rep(2, 2), col = c("blue", "red"), 
       bty = "n")

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)
deseq2_de_out <- results(deseq2_de)
summary(deseq2_de_out)

mcols(deseq2_de_out, use.names = TRUE)

## There are only 4 genes with an estimated FDR < 0.1.
deseq2_de_out_sig <- subset(deseq2_de_out, padj < 0.1)
nms_deseq2_sig <- rownames(deseq2_de_out_sig)

## Plots of (normalized) read counts for the top 4 genes.
par(mfrow = c(2, 2))
for(i in 1:4)
  plotCounts(deseq2_de, gene = rownames(deseq2_de_out_sig)[i], intgroup = "gender")

## Histogram of p-values. Not of an expected shape.
par(mfrow = c(1, 1))
hist(deseq2_de_out$pvalue)

## 83 of the p-values are listed as NA. According to the documentation, p-values of NA 
## are assigned for genes with low read counts or for genes that are deemed to be 
## outliers (via use of Cook's distance). Some of these genes are clearly low counts, but 
## others are not clearly outliers.
Y[is.na(deseq2_de_out$pvalue), ]

####
#### Analysis with limma.
####

## Design matrix.
limma_X <- cbind(1, rep(0:1, each = 3))

## Voom transformation and differential expression analysis.
limma_voom <- voom(Y, design = limma_X, plot = FALSE)
limma_fit <- lmFit(limma_voom, design = limma_X)
limma_eb <- eBayes(limma_fit)

## For some reason, no estimated FDRs fall below 0.95.
limma_FDRs <- topTable(limma_eb, n = nrow(Y), coef = 2)
quantile(limma_FDRs$adj.P.Val)

## Histogram of p-values. The picture suggests no differential expression.
hist(limma_FDRs$P.Value)

####
#### Summary of three methods.
#### 

## Comparison of p-value distributions.
par(mfrow = c(1, 3))
hist(edgeR_out$table$PValue, xlab = "P Values", main = "edgeR")
hist(deseq2_de_out$pvalue, xlab = "P Values", main = "DESeq2")
hist(limma_FDRs$P.Value, xlab = "P Values", main = "limma")

## edgeR selected 63 genes as significant, based on an estimated FDR threshold of 0.1. 
## DESeq2 only selected 4 genes. And limma didn't select any. Three of the four genes 
## selected by DESeq2 were also selected by edgeR. Of the genes selected as significant 
## by edgeR, 2 / 3 were assigned p-values of NA by DESeq2, presumably because they were 
## deemed to be outlier genes. Most of the other DESeq2 FDRs were very large for the 
## genes selected by edgeR.
nms_deseq2_sig %in% nms_edgeR_sig

deseq2_de_out[nms_edgeR_sig, ]$padj


