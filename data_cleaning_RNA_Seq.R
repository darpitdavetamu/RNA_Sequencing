rm(list=ls())
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
# Subsetting GRP data frame by RIN score for analysis on the subset of data
GRP6<-subset(GRP,RIN>6)
GRPg4.9<-subset(GRP,RIN>=4.9)
GRPl4.9<-subset(GRP,RIN<4.9)

GRP6$Sample<-droplevels(GRP6$Sample)
GRPg4.9$Sample<-droplevels(GRPg4.9$Sample)
GRPl4.9$Sample<-droplevels(GRPl4.9$Sample)


# Removing extra columns from countmatrix
DTA <- data.frame(DTA[,-c(1,2)], row.names=DTA[,2])


# Subset countmatrix to perform the necessary analysis

# Extract samples which have RIN score less than (condition)
col6<-GRP6$Sample
colg4.9<-GRPg4.9$Sample
coll4.9<-GRPl4.9$Sample


# Subset countmatrix where column name matches that which have a RIN score (condition)
DTA6<-subset(DTA,select = col6)
DTAg4.9<-subset(DTA,select = colg4.9)
DTAl4.9<-subset(DTA,select = coll4.9)


GRP <- data.frame(GRP[, -1], row.names = GRP[, 1])
GRP6 <- data.frame(GRP6[, -1], row.names = GRP6[, 1])
GRPg4.9 <- data.frame(GRPg4.9[, -1], row.names = GRPg4.9[, 1])
GRPl4.9 <- data.frame(GRPl4.9[, -1], row.names = GRPl4.9[, 1])



RIN <- data.frame(RIN[, -1], row.names = RIN[, 1])

colnames(DTA) <- substring(colnames(DTA), 2)
colnames(DTA6) <- substring(colnames(DTA6), 2)
colnames(DTAg4.9) <- substring(colnames(DTAg4.9), 2)
colnames(DTAl4.9) <- substring(colnames(DTAl4.9), 2)



DTA <- DTA[ ,order(colnames(DTA))]
DTA6 <- DTA[ ,order(colnames(DTA6))]
DTAg4.9 <- DTA[ ,order(colnames(DTAg4.9))]
DTAl4.9 <- DTA[ ,order(colnames(DTAl4.9))]


GRP <- GRP[order(rownames(GRP)), ]
GRP6 <- GRP[order(rownames(GRP6)), ]
GRPg4.9 <- GRP[order(rownames(GRPg4.9)), ]
GRPl4.9 <- GRP[order(rownames(GRPl4.9)), ]


## There are many genes that had no counts in any of the samples. Remove them 
## from the analysis.
n_count_merge <- drop(data.matrix(DTA) %*% rep(1, 27))
n_count_merge6 <- drop(data.matrix(DTA6) %*% rep(1, 16))
n_count_mergeg4.9 <- drop(data.matrix(DTAg4.9) %*% rep(1, 20))
n_count_mergel4.9 <- drop(data.matrix(DTAl4.9) %*% rep(1, 7))


DTA <- DTA[n_count_merge > 0, ]
DTA6 <- DTA6[n_count_merge > 0, ]
DTAg4.9 <- DTAg4.9[n_count_merge > 0, ]
DTAl4.9 <- DTAl4.9[n_count_merge > 0, ]



