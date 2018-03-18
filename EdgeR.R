rm(list=ls())
###############
#### EdgeR ####
###############

library(edgeR)

## Create a DGElist object for holding the counts along with the library sizes (the total 
## number of counts) and the group membership information for each sample. 
edgeR_dge <- DGEList(DTAl4.9, lib.size = colSums(DTAl4.9), group = GRPl4.9[, 1])

mydesign1<-model.matrix(~Age+Group+Sex+Age*Group, data=GRPl4.9)
D1 <- estimateGLMCommonDisp(edgeR_dge, mydesign1) 

fit1 <- glmFit(D1, mydesign1)
results<-glmLRT(fit1,coef = 5)
edgeR_out <- topTags(results, n = nrow(DTAl4.9), sort.by = "PValue")

edgeR_out_sig <- subset(edgeR_out$table, FDR < 0.1)
edgeR_out_n <- topTags(results, n = nrow(DTAl4.9), sort.by = "PValue")
edgeR_out_sig_n <- subset(edgeR_out_n$table, FDR < 0.05)

nms_edgeR_sig <- rownames(edgeR_out_sig)
rownames(DTA[nms_edgeR_sig, ])

nms_edgeR_sig_n <- rownames(edgeR_out_sig_n)
rownames(DTA[nms_edgeR_sig_n, ])



#design.add<-model.matrix(~Group+Sex+Age+Age*Group,data=GRP)
#edgeR_disp <- estimateCommonDisp(edgeR_dge)

#glm.stressadd<-glmFit(edgeR_disp, design.add)
#lrt.stressadd<-glmLRT(glm.stressadd)


#edgeR_out_n <- topTags(results, n = nrow(DTA), sort.by = "PValue")
#edgeR_out_sig_n <- subset(edgeR_out_n$table, FDR < 0.05)


#nms_edgeR_sig_n <- rownames(edgeR_out_sig_n)
#DTA[nms_edgeR_sig_n, ]
#######
#mydesign1<-model.matrix(~Age+Group+Sex+Age*Group, data=GRP)
#D1 <- estimateGLMCommonDisp(edgeR_dge, mydesign1) 
# D is the DGEList of the normalized counts
#D1 <- estimateGLMTrendedDisp(D1, mydesign1)
#fit1 <- glmFit(D1, mydesign1)
#results<-glmLRT(fit1,coef = 5)
#edgeR_out <- topTags(results, n = nrow(DTA), sort.by = "PValue")
#edgeR_out_sig <- subset(edgeR_out$table, FDR < 0.1)
#edgeR_out_n <- topTags(results, n = nrow(DTA), sort.by = "PValue")
#edgeR_out_sig_n <- subset(edgeR_out_n$table, FDR < 0.05)

#nms_edgeR_sig <- rownames(edgeR_out_sig)
#DTA[nms_edgeR_sig, ]

#nms_edgeR_sig_n <- rownames(edgeR_out_sig_n)
#DTA[nms_edgeR_sig_n, ]

#design.add<-model.matrix(~Group+Sex+Age,data=GRP)
#edgeR_disp <- estimateCommonDisp(edgeR_dge)
#glm.stressadd<-glmFit(edgeR_disp, design.add)
#design.add
#lrt.stressadd<-glmLRT(glm.stressadd)
#topTags(lrt.stressadd)
#edgeR_out_n <- topTags(results, n = nrow(DTA), sort.by = "PValue")
#edgeR_out_sig_n <- subset(edgeR_out_n$table, FDR < 0.05)


#nms_edgeR_sig_n <- rownames(edgeR_out_sig_n)
#DTA[nms_edgeR_sig_n, ]
