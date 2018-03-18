df1<-read.csv("1_1PP.txt",sep="\t",header=TRUE)
df2<-read.csv("4MH.txt",sep="\t",header=TRUE)

new<-rbind(df1[,1],df2[,1])

Reduce(union, list(x, y, z))


x <- scan("1_1PP.txt", sep="\t",header=TRUE)
