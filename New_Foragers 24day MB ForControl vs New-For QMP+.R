library("DESeq2")
setwd("E:/24.03.2020/Dokumente/Experiment/2019.04.09 Honeybee Cage experiment/Results/Result_geneid/New-Foragers 24day old/R analysis/New_Foragers 24day MB")
sampleNames <- c("NeFA1MB_QMP","NeFA2MB_QMP","NeFC1MB_QMP","NeFC2MB_QMP","NeFD2MB_QMP","NeFD3MB_QMP",
                 "CoFA1MB","CoFA4MB","CoFC1MB","CoFC2MB","CoFD1MB","CoFD2MB")
mydata <- read.table("New_Foragers 24day MB ForControl vs New-For QMP+.txt", header = TRUE)
names(mydata)[2:13] <- sampleNames
countMatrix <- as.matrix(mydata[2:13])
rownames(countMatrix) <-mydata$Geneid
table2 <- data.frame(name = c("QMP1","QMP2","QMP3","QMP4","QMP5","QMP6",
                              "ColFor1","ColFor2","ColFor3","ColFor4","ColFor5","ColFor6"),
                     condition = c("QMP","QMP","QMP","QMP","QMP","QMP",
                                   "ColFor","ColFor","ColFor","ColFor","ColFor","ColFor"),
                     colony = c("A","A","C","C","D","D",
                                "A","A","C","C","D","D"))
rownames(table2) <- sampleNames
#########Filter:The number of counts beyond 9 at least 5 samples
countMatrix_filter1<-countMatrix[rowSums(countMatrix>9)>=5,]
head(countMatrix_filter1)
dds <- DESeqDataSetFromMatrix(countMatrix_filter1, colData=table2, design= ~ colony+condition)
dds<-DESeq(dds,test="LRT", reduced = ~colony)

res <- results(dds)
write.table(res,"New_Foragers 24day MB ForControl vs New-For QMP+ result.csv", sep = ",", row.names = TRUE)
head(res)
summary(res)
rld <- rlog(dds)
#plotPCA(rld, intgroup=c("condition","colony"), ntop=length(rownames(matrix_cook)))

library(ggplot2)
data1 <- plotPCA(rld,intgroup=c("condition","colony"),ntop=length(rownames(countMatrix_filter1)),returnData=TRUE)
percentVar <- round(100*attr(data1, "percentVar"))
p<-ggplot(data1,aes(PC1,PC2,color=condition, shape=colony))+
  scale_color_manual(values = c("#7fbf7b","#91bfdb"))+
  geom_point(size=8)+
  xlab(paste0("PC1:",percentVar[1],"%variance"))+
  ylab(paste0("PC2:",percentVar[2],"%variance"))
p
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p1<-p+geom_point()+theme
p1
theme2<-theme(axis.text.x=element_text(size=26,face = "plain",hjust = 0.5),
              axis.text.y=element_text(size=26,face = "plain",hjust = 0.5))
p2<-p1+theme2
p2


pn<-p2+geom_text(aes(label=sampleNames),vjust=1.5,colour='red')
pn


library("pheatmap")
select<-
  order(rowMeans(counts(dds,normalized=TRUE)),decreasing = TRUE)[1:1000]
nt<-normTransform(dds)
log2.norm.counts<-assay(nt)[select,]
df<-as.data.frame(colData(dds)[,c("name","condition")])
pdf('heatmap1000.pdf',width = 7,height = 8)
pheatmap(log2.norm.counts,cluster_rows = TRUE,show_rownames = FALSE,cluster_cols = TRUE,annotation_col = df)








