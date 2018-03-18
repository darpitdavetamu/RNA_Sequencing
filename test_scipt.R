####
#### Data cleaning
####

DTA <- read.csv('merged_count.csv')
GRP <- read.csv("GroupInfor.csv", header = TRUE,sep="\t")
RIN <- read.csv("miRNA_RIN.csv", header = TRUE,sep=",")


#Remove any characters in Sample column # Done for matching rows
RIN$Sample<-gsub("[[:punct:]]", "_", RIN$Sample)


# Converting all lower case to upper case
RIN<-data.frame(lapply(RIN, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))


#Merge data frames to extract RIN score
GRP<-merge(x = GRP, y = RIN, by.x='Sample', by.y = 'Sample')


DTA <- data.frame(DTA[,-c(1,2)], row.names=DTA[,2])
GRP <- data.frame(GRP[, -1], row.names = GRP[, 1])
RIN <- data.frame(RIN[, -1], row.names = RIN[, 1])

colnames(DTA) <- substring(colnames(DTA), 2)

DTA <- DTA[, order(colnames(DTA))]
GRP <- GRP[order(rownames(GRP)), ]
RIN<-RIN[order(rownames(RIN)),]

## There are many genes that had no counts in any of the samples. Remove them 
## from the analysis.
n_count_merge <- drop(data.matrix(DTA) %*% rep(1, 27))
DTA <- DTA[n_count_merge > 0, ]



###################################
#### Analysis with DESeq2.#########
###################################



library(DESeq2)

## Create a DESeqDataSet object for holding the counts and phenotype data.
deseq2_deseqds <- DESeqDataSetFromMatrix(countData = DTA, colData = GRP, design = ~ Group+Sex+Age)

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
par(mfrow = c(1, 1))
hist(deseq2_de_out$pvalue)

###############
#### EdgeR ####
###############

library(edgeR)

## Create a DGElist object for holding the counts along with the library sizes (the total 
## number of counts) and the group membership information for each sample. 
edgeR_dge <- DGEList(DTA, lib.size = colSums(DTA), group = GRP[, 2])

mydesign1<-model.matrix(~Age+Group+Sex+Age*Group, data=GRP)
D1 <- estimateGLMCommonDisp(edgeR_dge, mydesign1) 

# D is the DGEList of the normalized counts
#D1 <- estimateGLMTrendedDisp(D1, mydesign1)

fit1 <- glmFit(D1, mydesign1)
results<-glmLRT(fit1,coef = 2)
edgeR_out <- topTags(results, n = nrow(DTA), sort.by = "PValue")
edgeR_out_sig <- subset(edgeR_out$table, FDR < 0.1)

nms_edgeR_sig <- rownames(edgeR_out_sig)
DTA[nms_edgeR_sig, ]


design.add<-model.matrix(~Group+Sex+Age,data=GRP)
edgeR_disp <- estimateCommonDisp(edgeR_dge)

glm.stressadd<-glmFit(edgeR_disp, design.add)
lrt.stressadd<-glmLRT(glm.stressadd)


edgeR_out_n <- topTags(results, n = nrow(DTA), sort.by = "PValue")
edgeR_out_sig_n <- subset(edgeR_out_n$table, FDR < 0.05)


nms_edgeR_sig_n <- rownames(edgeR_out_sig_n)
DTA[nms_edgeR_sig_n, ]
