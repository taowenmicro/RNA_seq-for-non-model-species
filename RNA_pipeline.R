#设置工作目录，所有文件放到此目录下
setwd("E:/Shared_Folder/RNA_pipline/ref_M")

path = getwd()
path = paste(path,"/result",sep = "")

dir.create(path)

mytheme <- theme_bw()+
  #theme_classic()+
  # scale_color_manual(values = mi, guide = guide_legend(title = NULL))+
  # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
  theme(
    
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    
    plot.title = element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
    axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
    axis.text = element_text(size = 20,face = "bold"),
    axis.text.x = element_text(colour = "black",size = 24),
    axis.text.y = element_text(colour = "black",size = 24),
    legend.text = element_text(size = 20,face = "bold"))+
  #theme(legend.position = c(0.1,0.2))+
  
  theme(strip.text.x = element_text(size=15, angle=0),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="blue", fill="#CCCCFF"))


library("DESeq2")  
#-----------------------------------------导入数据---------------------------------------------------
countdata <- read.delim("CountMatrix.txt",row.names = 1)  
head(countdata, 10)  
#构造分组文件
SampleType<- sapply(strsplit(colnames(countdata), "_"), `[`, 1)
SampleType
coldata =data.frame(ID = colnames(countdata),SampleType)
head(coldata)

Desep_group <- as.character(levels(coldata$SampleType))
Desep_group

aaa = combn(Desep_group,2)




## step1 封装成DESeqDataSetFromMatrix对象
dds <- DESeqDataSetFromMatrix(countData = countdata,colData = coldata,design = ~ SampleType)  

##--查看DESeqDataSet对象  
dim(dds)  
assay(dds)#count表格 
assayNames(dds)  
colSums(assay(dds))#每个样本read数目统计  
rowRanges(dds)  
colData(dds) #mapping文件 
nrow(dds) #查看基因数量

#过滤没有reads比对上的基因，所有reads数为零  
dds <- dds[rowSums(counts(dds)) > 1,]  
nrow(dds)  

#--------------------数据整理------------------------------------------------------------
# 将数据通过rolg方法与vst方法转换，这样可以用于后面计算距离矩阵  
## ----rlog方法-------
rld <- rlog(dds, blind = FALSE)  
head(assay(rld), 3)  
## ----vst方法------
vsd <- vst(dds, blind = FALSE)  
head(assay(vsd), 3)  



######---------------------------------------基于矩阵差异-整体样本差异--------------------------
#利用转换后的结果计算样品之间距离关系  
#方法1--欧氏距离
sampleDists <- dist(t(assay(rld)))  
sampleDists  

library("pheatmap")  
library("RColorBrewer")  
library("ggplot2")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5-----
sampleDistMatrix <- as.matrix( sampleDists )  
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )  
# colnames(sampleDistMatrix) <- NULL  
p = pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,  
         clustering_distance_cols = sampleDists)  

p

filename = paste(path,"/a1_sample_heatmap.pdf",sep = "")
ggsave(filename, p, width = 6, height =6 )

##----------------------------------------------------------PCA-----------------------------------------
plotPCA(rld, intgroup = c("SampleType"))  
pcaData <- plotPCA(rld, intgroup = c( "SampleType"), returnData = TRUE)  
pcaData  
percentVar <- round(100 * attr(pcaData, "percentVar"))  
library(ggplot2)  
p2 = ggplot(pcaData, aes(x = PC1, y = PC2, color = SampleType)) +  
  geom_point(size =3) +  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  # stat_ellipse( linetype = 2,level = 0.65,aes(group  =SampleType, colour =  SampleType))+
  coord_fixed()

p2 = p2 + mytheme
p2

filename = paste(path,"/a1_sample_PCA.pdf",sep = "")
ggsave(filename, p2, width = 8, height =6 )

##----------------------------------------DEsep2差异基因分析-----------------------------------------for 循环---------
# 差异表达计算  
library(DESeq2)
dep <- DESeq(dds)  
#导出标准化表格
count_nor = as.data.frame(counts(dep,normalize=TRUE))
head(count_nor)

#查看差异分析分组及其矩阵
resultsNames(dep)
unclass(dep)


# res <- results(dep)  
# res  
i= 1
for (i in 1:dim(aaa)[2]) {
  
  Desep_group = aaa[,i]

# Desep_group <- c("Y3" , "G1" )
res <-  results(dep, contrast=c("SampleType",Desep_group ),alpha=0.05)
# summary一下，看一下结果的概要信息
summary(res)
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
WT <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

dim(WT)
res$level = as.factor(ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "enriched",ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "depleted","nosig")))
# dim(res)
# head(res)

#保存差异分析全部结果文件
group = paste(Desep_group[1],Desep_group[2],sep = "-")


filename = paste(path,"/","a1_DEsep2-all-",group,".csv",sep = "")
write.csv(x = res,file = filename)  


library(tidyverse)
library("dplyr")
#筛选出p值小于0.05的基因  
res.05 <- results(dep, alpha = 0.05)  
table(res.05$padj < 0.05)  
head(res.05)
res.05$ID = row.names(res.05)
res.05  = as.data.frame(res.05 )
res_filter<- filter(res.05, padj < 0.05)
dim(res_filter)


head(countdata, 10) 
countdata$ID = row.names(countdata)
sub_gene = filter(countdata, row.names(countdata) %in% res_filter$ID)
head(sub_gene)


# filename = paste(path,"/a1_DEsep2_05.csv",sep = "")
# write.csv(x = sub_gene,file = filename,row.names = T)  

# write.csv(x =sub_gene,file = "res.05.csv",row.names = T)  

#统计p值小于0.05差异表达基因数目  
sum(res$pvalue < 0.05, na.rm=TRUE)  
sum(!is.na(res$pvalue))  
sum(res$padj < 0.1, na.rm=TRUE)  

#筛选出差异表达明显的基因Significant，设定标准为p值小于0.01，至于使用0.05还是0.01，具体问题具体分析  
resSig <- subset(res, padj < 0.05)  
#按log2差异倍数排序，先升序，设置decreasing = TRUE降序  
head(resSig[ order(resSig$log2FoldChange), ])  
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ]) 

filename = paste(path,"/a1_DEsep2_05_order",group,".csv",sep = "")
write.csv(resSig, file = filename)

#-------------------------------------------------------展示差异-----------
m <- res
head(m)    
m <- na.omit(m)    
m <- transform(m,padj=-1*log10(m$padj))    
down <- m[m$log2FoldChange<=-1,]     
up <- m[m$log2FoldChange>=1,]    
no <- m[m$log2FoldChange>-1 & m$log2FoldChange <1,]     
filename3 = paste(path,"/a2_DEsep2_Volcano_plot-",group,".pdf",sep = "")
pdf(file=filename3,width = 8, height = 6)
plot(no$log2FoldChange,no$padj,xlim = c(-10,10),ylim=c(0,100),col="blue",pch=16,cex=0.8,main = "Gene Expression",xlab = "log2FoldChange",ylab="-log10(Qvalue)")    
points(up$log2FoldChange,up$padj,col="red",pch=16,cex=0.8)    
points(down$log2FoldChange,down$padj,col="green",pch=16,cex=0.8) 
dev.off()
}
#-------------------------------------------------------------------------------------------------for 循环结束---------


##--------------------富集分析------------------------
#加载各种包，如果加载失败就自行安装，注意除了前两个，都使用BiocManager::install()进行安装。
library(dplyr)
library(tidyr)
library(DOSE)
library(BiocManager)
library(GO.db)
# library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)

SampleType<- sapply(strsplit(colnames(countdata), "_"), `[`, 1)
SampleType
coldata =data.frame(ID = colnames(countdata),SampleType)
head(coldata)

Desep_group <- as.character(levels(coldata$SampleType))
Desep_group
aaa = combn(Desep_group,2)
aaa


path = getwd()
path = paste(path,"/result",sep = "")

dir.create(path)

path_GO = paste(path,"/result_GO",sep = "")
dir.create(path_GO)


data <- read.table(
  "melon_v4_GO_anno.txt",
  header = F,
  sep = "\t")
go2gene <- data[, c(2, 1)]
head(go2gene)
user_data = buildGOmap(go2gene)
bg = user_data[, c(2, 1)]

length(unique(bg$GO))
#----------------构建人造数据库----------------------------------------
#-每个基因可能对应有多个GO  ID，因为有多个等级的GO通路等级-------------
library(GO.db)
# write.csv(goname_all,"./goname_GO.csv")
goname_all <- AnnotationDbi::select(x=GO.db, keys = as.character(unique(bg$GO)),  keytype = "GOID",columns = c("TERM","ONTOLOGY","DEFINITION") )
dim(goname_all)


i = 1
for (i in 1:dim(aaa)[2]) {
Desep_group = aaa[,i]
group = paste(Desep_group[1],Desep_group[2],sep = "-")
group
filename = paste(path,"/a1_DEsep2_05_order",group,".csv",sep = "")

gene <- read.csv( filename ,header = T,row.names = 1,stringsAsFactors = F)
dim( gene)
geneList = sort(row.names(gene), decreasing = TRUE)
geneList = as.factor(geneList)
geneList 
# write.csv(geneList,"./genglist.csv",row.names = F,col.names = T)

goname_CC = filter(goname_all, ONTOLOGY == "CC")
dim(goname_CC)
# write.csv(goname_CC,"./goname_GO_CC.csv")
bg1 = bg[,c(2,1)]
head(bg1)
bg_CC = filter(bg1, GO  %in% goname_CC$GOID)
head(bg_CC)

x <- enricher(
  geneList,
  TERM2GENE =bg1,
  TERM2NAME = goname_all,
  pvalueCutoff = 0.2,
  minGSSize = 10, 
  maxGSSize = 500,
  qvalueCutoff  = 0.9,
  pAdjustMethod = "BH"
  # pvalueCutoff  = 0.2,
  
  )
  
x

##一共富集了几个通路，有多少个显著
aa = x@result[1:2]
head(aa)

length(aa$ID)

length(unique(aa$ID))

aaaa <- AnnotationDbi::select(x=GO.db, keys = as.character(aa$ID),  keytype = "GOID",columns = c("TERM","ONTOLOGY","DEFINITION") )
head(aaaa)

x@result$Description = aaaa$TERM
x@result$ONTOLOGY = aaaa$ONTOLOGY

p = barplot(x, title=group, showCategory = 62) 
list(p)
###出图
plot_GO$Count
plot_GO = subset(x@result,qvalue < 0.9  & p.adjust < 0.2)
dim(plot_GO )
colnames(plot_GO)
str(plot_GO)
##排序
plot_GO<- arrange(plot_GO, ONTOLOGY, Count)
head(plot_GO)
plot_GO$Description = factor(plot_GO$Description,levels = plot_GO$Description)

library(ggplot2)
p = ggplot(plot_GO,aes( x = Count, y = Description,fill =ONTOLOGY,size = Count )) +
  geom_point(pch = 21) 
p = p + mytheme
p

FileName2 <- paste(path_GO,"/GO_all_point",group,".pdf", sep = "")

ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)

p = ggplot(plot_GO,aes(y = Count, x = Description,fill =ONTOLOGY )) +
  geom_bar(stat = "identity") + coord_flip()
p
p = p + mytheme
p

FileName2 <- paste(path_GO,"/GO_all_bar",group,".pdf", sep = "")

ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)



# dotplot(x,title=group)#点图，按富集的数从大到小的
filename = paste(path_GO,"/a1_DEsep2_05_order",group,"_GO.csv",sep = "")
write.csv(x@result,filename)
wt_cs = x@result
x@result = filter(wt_cs, ONTOLOGY == "CC")
x@ontology = "CC"
#DNG
# cnetplot(ego, foldChange=geneList)
library("topGO")
library("Rgraphviz")

# groupGO()
# ?groupGO

filename3 = paste(path_GO,"/GO_DNG_plot-CC-",group,".pdf",sep = "")
pdf(file=filename3,width = 8, height = 6)
plotGOgraph(x)
dev.off()

x@result = filter(wt_cs, ONTOLOGY == "BP")
x@ontology = "BP"
#DNG
filename3 = paste(path_GO,"/GO_DNG_plot-BP-",group,".pdf",sep = "")
pdf(file=filename3,width = 8, height = 6)
plotGOgraph(x)
dev.off()
colnames(wt_cs)
x@result = filter(wt_cs, ONTOLOGY == "MF")
x@ontology = "MF"
#DNG
filename3 = paste(path_GO,"/GO_DNG_plot-MF-",group,".pdf",sep = "")
pdf(file=filename3,width = 8, height = 6)
plotGOgraph(x)
dev.off()

}


##------------------------------------------------KEGG富集------------------------------------
##甜瓜数据库:cmo


SampleType<- sapply(strsplit(colnames(countdata), "_"), `[`, 1)
SampleType
coldata =data.frame(ID = colnames(countdata),SampleType)
head(coldata)

Desep_group <- as.character(levels(coldata$SampleType))
Desep_group
aaa = combn(Desep_group,2)
aaa


path = getwd()
path = paste(path,"/result",sep = "")

dir.create(path)

path_KEGG = paste(path,"/result_KEGG",sep = "")
dir.create(path_KEGG)


wt_gene = read.delim("./kegg_background",header = F,row.names = 1)
head(wt_gene)
aaab
i = 1

SampleType<- sapply(strsplit(colnames(countdata), "_"), `[`, 1)
SampleType
coldata =data.frame(ID = colnames(countdata),SampleType)
head(coldata)

Desep_group <- as.character(levels(coldata$SampleType))
Desep_group

aaab = combn(Desep_group,2)
aaab




for (i in 1:dim(aaab)[2]) {
  Desep_group = aaab[,i]
  group = paste(Desep_group[1],Desep_group[2],sep = "-")
  group
  filename = paste(path,"/a1_DEsep2_05_order",group,".csv",sep = "")
  



dif = read.csv(filename,row.names = 1)  
head(dif)

index = merge(dif,wt_gene,by = "row.names",all = F)
dim(index)
#去除重复，因为会有多个基因注释到同一个KO上
index = distinct(index, V2, .keep_all = TRUE)  
head(index)
dim(index) 



kk <- enrichKEGG(gene = index$V2,
                 organism='ko', keyType='kegg',
                 pvalueCutoff = 0.2)
head(kk,6)

colnames(kk@result)

plot_GO = subset(kk@result,qvalue < 0.9  & p.adjust < 0.05)
dim(plot_GO)
colnames(plot_GO)

##排序
plot_GO<- arrange(plot_GO,  Count)
head(plot_GO)
plot_GO$Description = factor(plot_GO$Description,levels = plot_GO$Description)

library(ggplot2)
p = ggplot(plot_GO,aes( x = Count, y = Description,size = Count )) +
  geom_point(pch = 21) 
p = p + mytheme
p

FileName2 <- paste(path_KEGG,"/KEGG_all_point",group,".pdf", sep = "")

ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)

p = ggplot(plot_GO,aes(y = Count, x = Description )) +
  geom_bar(stat = "identity") + coord_flip()
p
p = p + mytheme
p

FileName2 <- paste(path_KEGG,"/KEGG_all_bar",group,".pdf", sep = "")

ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)



# install("pathview")
library("pathview")
# dotplot(kk,title="Enrichment KEGG_dot")
head(index)

aaa = index[,c(3,9)]
head(aaa)
str(aaa)
# aaa = distinct(aaa, V2, .keep_all = TRUE)  
dim(aaa)

aaa<- arrange(aaa, V2)
head(aaa)


# row.names(index) = gsub("K","",index$V2)
# index$Row.names = NULL
# aa = index[,3]
# length(index$V2)
# unique(index$V2)
# row.names(aaa) = gsub("K","",aaa$V2)
row.names(aaa) = aaa$V2
aaa$V2 = NULL

aaa = as.matrix(aaa)
head(aaa)
aaa[,1]




#http://www.bioinfo-scrounger.com/archives/639 代谢通路整合转录组和代谢组

path_ko_plot = paste(path_KEGG,"/Ko",group,sep = "")
dir.create(path_ko_plot)
setwd(path_ko_plot )
tmp = sapply(gsub("ko","",plot_GO$ID[1:10]), function(pid) pathview(gene.data=aaa[,1],
                                                                    pathway.id=pid,
                                                                    species="ko",
                                                                    kegg.native = F,
                                                                    kegg.dir = path_ko_plot))


setwd(path)
setwd("../")
}


##------------------------------GSEA------------------------------------------------------------

#-------------------------------------------------GSVA-------------------------------------
head(count_nor)
wt_gene = read.delim("./kegg_background",header = F,row.names = 1)
head(wt_gene)

GS_tb = merge(count_nor,wt_gene,by = "row.names",all = F)
dim(GS_tb)
#去除重复，因为会有多个基因注释到同一个KO上
library(tidyverse)
GS_tb = distinct(GS_tb, V2, .keep_all = TRUE)  
head(GS_tb)
dim(GS_tb) 
row.names(GS_tb) = GS_tb$V2
GS_tb$V2 = NULL
GS_tb$Row.names = NULL
GS_tb = as.matrix(GS_tb)

library(clusterProfiler)
ko_kegg <- clusterProfiler::download_KEGG('ko')

head(ko_kegg )
names(ko_kegg)
PATH2ID <- ko_kegg$KEGGPATHID2EXTID
PATH2NAME <- ko_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')

head(PATH_ID_NAME )
##查看有多少个ko号
length(unique(PATH_ID_NAME$KEGGID))
#转化成list
group <- split(PATH_ID_NAME$KO,PATH_ID_NAME$KEGGID)

#去除重复
KO_tax = distinct(PATH_ID_NAME, KEGGID, .keep_all = TRUE)  
head(KO_tax)
names(group ) = KO_tax$DESCRPTION


path = getwd()
path = paste(path,"/result",sep = "")
path_GSVA = paste(path,"/result_GSVA",sep = "")
dir.create(path_GSVA)
path_GSVA
#Perform GSVA   
topMatrixGSVA <- gsva(GS_tb, group, min.sz=10, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)
dim(topMatrixGSVA)
filename = paste(path_GSVA,"/a1_table_GSVA",".csv",sep = "")
write.csv(topMatrixGSVA , file = filename)
# design <- model.matrix( ~ factor(coldata$SampleType, levels=c("G1", "G2")))
# colnames(design) <- c("Intercept", "G1VsG2")
# design
design =model.matrix(~ 0 + coldata$SampleType)
colnames(design)=levels(coldata$SampleType)
fit <- lmFit(topMatrixGSVA, design)
fit <- eBayes(fit)

Desep_group <- as.character(levels(coldata$SampleType))
Desep_group

aaa = combn(Desep_group,2)

filename3 = paste(path_GSVA,"/GSVE_",".pdf",sep = "")
pdf(file=filename3,width = 18, height = 15)
heatmap.2(heat, col=myCol, breaks=myBreaks, main="Title", key=TRUE, keysize=1.0, key.title="", key.xlab="Enrichment Z-score", scale="none", density.info="none", reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), trace="none", cexRow=1.0, cexCol=1.0, distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"), margin=c(10,25))
dev.off()

i = 1
for (i in 1:dim(aaa)[2]) {
  sds = paste(aaa[1,i],aaa[2,i],sep = "-")
  contrast.matrix<- makeContrasts(sds, levels=design)
  fit2 <-contrasts.fit(fit, contrast.matrix)
  fit2 <-eBayes(fit2)
  #"fdr"(equivalent to "BH");lfc(log2-fold-change )
  results<-decideTests(fit2,method="global", adjust.method="BH", p.value=0.05, lfc=0)
  summary(results)
  x<-topTable(fit2,coef=1, number=10000, adjust.method="BH", sort.by="p",resort.by=NULL)
  sigPathways <- x
  dim(sigPathways)
  WT <-subset(sigPathways,adj.P.Val < 0.05 )
  dim(WT)
  
  
  filename = paste(path_GSVA,"/a1_limma_GSVA_",sds,".csv",sep = "")
  write.csv(sigPathways , file = filename)
  
  # sigPathways <- topTable(fit2, coef="DiseaseVsControl", number=Inf, p.value=0.05, adjust="BH")
  # sigPathways <- sigPathways[abs(sigPathways$logFC)>1,]
  res <- results
  # head(res)
  #Create an object that can easily be written to disc
  wObject <- data.frame(rownames(sigPathways), sigPathways)
  colnames(wObject) <- c("Pathway","Log2FoldChange","MeanExpression","tStat","Pvalue","AdjustedPval","Bvalue")
  
  #Filter the GSVA object to only include significant pathways
  topMatrixGSVA <- topMatrixGSVA[rownames(sigPathways),]
  dim(topMatrixGSVA )
  
  
  
  
  #Set colour for heatmap
  require(RColorBrewer)
  myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
  myBreaks <- seq(-1.5, 1.5, length.out=101)
  heat <- t(scale(t(topMatrixGSVA)))
  head(heat)
  heat = as.data.frame(heat)
  head(WT)
  heat$id = row.names(heat)
  library(tidyverse)
  sub_heat<- filter(heat,id %in% row.names(WT))
  
  sub_map<- filter(coldata,SampleType %in% aaa[,i])
  
  
  
  vars<- c(as.character(sub_map$ID),"id")
  sub_heat = dplyr::select(sub_heat,one_of(vars))
  head(sub_heat) 
  row.names(sub_heat) = sub_heat$id
  sub_heat$id = NULL
  sub_heat = as.matrix(sub_heat)
  
  library(gplots)
  filename3 = paste(path_GSVA,"/GSVE_",sds,".pdf",sep = "")
  pdf(file=filename3,width = 18, height = 15)
  heatmap.2(sub_heat, col=myCol, breaks=myBreaks, main="Title", key=TRUE, keysize=1.0, key.title="", key.xlab="Enrichment Z-score", scale="none", density.info="none", reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), trace="none", cexRow=1.0, cexCol=1.0, distfun=function(x) dist(x, method="euclidean"), hclustfun=function(x) hclust(x, method="ward.D2"), margin=c(10,25))
  dev.off()
  
}

###-----------------------------------GSEA--------------------------------------------------------------
path = getwd()
path = paste(path,"/result",sep = "")
path_GSEA = paste(path,"/result_GSEA",sep = "")
dir.create(path_GSEA)



for (i in 1:dim(aaa)[2]) {
  
  
  library(clusterProfiler)
  # for (i in 1:dim(aaa)[2]) {
  Desep_group = aaa[,i]
  group = paste(Desep_group[1],Desep_group[2],sep = "-")
  group
  filename = paste(path,"/a1_DEsep2_05_order",group,".csv",sep = "")
  
  gene <- read.csv( filename ,header = T,row.names = 1,stringsAsFactors = F)
  head( gene)
  gene$GSEA = row.names(gene)
  gene_GSEA<- arrange(gene,  desc(log2FoldChange))
  head(gene_GSEA)
  geneList_GSE = factor(gene_GSEA$GSEA,levels = gene_GSEA$GSEA)
  
  # write.csv(geneList,"./genglist.csv",row.names = F,col.names = T)
  
  goname_CC = filter(goname_all, ONTOLOGY == "CC")
  dim(goname_CC)
  # write.csv(goname_CC,"./goname_GO_CC.csv")
  bg1 = bg[,c(2,1)]
  head(bg1)
  bg_CC = filter(bg1, GO  %in% goname_CC$GOID)
  head(bg_CC)
  
  geneList =gene_GSEA[,2]
  ## feature 2: named vector
  names(geneList) = as.character(gene_GSEA[,8])
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  
  head(goname_CC)
  geneList
  # GSEA富集分析
  x <- GSEA(
    geneList,
    TERM2GENE =bg1,
    TERM2NAME = goname_CC,
    pvalueCutoff = 0.2)
  
  if(dim(x@result)[1] == 0){
    filename = paste(path_GSEA,"/a1_table_GSEA",group,".csv",sep = "")
    write.csv(x@result , file = filename)
  }else{
    
    filename = paste(path_GSEA,"/a1_table_GSEA",group,".csv",sep = "")
    write.csv(x@result , file = filename)
    
    
    ##----下面可视化通络，通路上调和下调
    
    qw = rep("up",length(x@result$ID))
    for (i in 1:length(x@result$ID)) {
      if (x@result$NES[i] > 0) {
        qw[i] = "up"
      }else{
        qw[i] = "down"
      }
    }
    x@result$.sign = qw
    
    aa = x@result[1:2]
    head(aa)
    
    length(aa$ID)
    
    length(unique(aa$ID))
    library(GO.db)
    aaaa <- AnnotationDbi::select(x=GO.db, keys = as.character(aa$ID),  keytype = "GOID",columns = c("TERM","ONTOLOGY","DEFINITION") )
    head(aaaa)
    
    x@result$Description = aaaa$TERM
    x@result$ONTOLOGY = aaaa$ONTOLOGY
    
    
    p = dotplot(x,showCategory=30,split = ".sign") + facet_grid(~.sign)
    
    p
    FileName2 <- paste(path_GSEA,"/GSEA_all_point",group,".pdf", sep = "")
    
    ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)
    
    
    p = ridgeplot(x, 30)
    
    FileName2 <- paste(path_GSEA,"/GSEA_all_bridge",group,".pdf", sep = "")
    ggsave(FileName2, p, width = 18, height =24, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)
    
    
    
    # gseaplot(x, 1) #输出第一个结果
    
    library("enrichplot")
    # gseaplot2(x, 1)
    p = gseaplot2(x, 1:3, pvalue_table = TRUE)
    
    p
    FileName2 <- paste(path_GSEA,"/GSEA_all_gseaplot2",group,".pdf", sep = "")
    ggsave(FileName2, p, width = 8, height =6, device = cairo_pdf, family = "Times New Roman" ,limitsize = FALSE)
    
  }
}



