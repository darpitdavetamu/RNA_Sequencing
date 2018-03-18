rm(list=ls())
####
#### Analysis with limma.
####

## Design matrix.
#limma_dge <- DGEList(DTA, lib.size = colSums(DTA),group=GRP[,2])
limma_X<-model.matrix(~Age+Group+Sex+Age*Group,data=GRP)

#GRP<-GRP[order(GRP$Group),]
#limma_X <- cbind(1, rep(0:1, each = 15))


## Voom transformation and differential expression analysis.
limma_voom <- voom(DTA, design = limma_X)
limma_fit <- lmFit(DTA, design = limma_X)
limma_eb <- eBayes(limma_fit)

## For some reason, no estimated FDRs fall below 0.95.
limma_FDRs <- topTable(limma_eb, n = nrow(DTA), coef = c(3,5))
#quantile(limma_FDRs$adj.P.Val)

## Histogram of p-values. The picture suggests no differential expression.
#hist(limma_FDRs$P.Value)


# Display rownames with FDR (condition)
rownames(limma_FDRs[limma_FDRs$adj.P.Value<0.1,])

