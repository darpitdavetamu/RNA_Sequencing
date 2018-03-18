rm(list=ls())
###################################
#### Analysis with DESeq2.#########
###################################



library(DESeq2)

## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = DTA, colData = GRP, design = ~ Age+Group+Sex+Age*Group+Sex*Age+RIN)

## Apply rlog transformation.
deseq2_rld <- rlog(deseq2_deseqds)

## Differential expression analysis. 
deseq2_de <- DESeq(deseq2_deseqds)
deseq2_de_out <- results(deseq2_de)
summary(deseq2_de_out)

mcols(deseq2_de_out, use.names = TRUE)

## There are only 4 genes with an estimated FDR < 0.1.
deseq2_de_out_sig_005 <- subset(deseq2_de_out, padj < 0.05)
(nms_deseq2_sig_005 <- rownames(deseq2_de_out_sig_005))

deseq2_de_out_sig_01 <- subset(deseq2_de_out, padj < 0.1)
(nms_deseq2_sig_01 <- rownames(deseq2_de_out_sig_01))


## Histogram of p-values. Not of an expected shape.
#par(mfrow = c(1, 1))
#hist(deseq2_de_out$pvalue)

