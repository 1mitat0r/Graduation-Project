rm(list=ls())

#获取一些包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install(version = "3.10")
#BiocManager::install("cgdsr")
BiocManager::install("TCGAbiolinks")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("DESeq2")
BiocManager::install("WGCNA")
BiocManager::install(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))



if (!requireNamespace("survminer", quietly = TRUE))
  install.packages("survminer")


if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")


library(BiocManager)
library("TCGAbiolinks")
library("survival")
library("survminer")
library("tidyverse")
library("ggplot2")
#library("ggpubr")
#library(org.Hs.eg.db)
library(AnnotationDbi)
library("DESeq2")
library("pheatmap")
#library("cgdsr")
library(WGCNA)

## 1.差异表达分析，16个阳性16个阴性
query_mRNA <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts"
                  )
### 获得癌症样本编号
samplesDown <- getResults(query_mRNA,cols=c("cases")) 

#TP代表PRIMARY SOLID TUMOR；
#NT代表Solid Tissue Normal

####从samplesDown中筛选出TP样本的barcodes
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

####从samplesDown中筛选出NT样本的barcodes
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
####TP和NT太多了,选择其中的一些
dataSmNT <- dataSmNT[1:16]#选择前16个TP
dataSmTP <- dataSmTP[1:16]#选择前16个NT
### 最后设置barcodes参数
queryDown_mRNA <-GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          experimental.strategy = "RNA-Seq",
                          workflow.type = "STAR - Counts",
                          barcode = c(dataSmTP, dataSmNT)
                          )
###下载
GDCdownload(queryDown_mRNA,
            method = "api",
            directory = "GDCdata1",
            files.per.chunk = 6
            )  
###创建SummarizedExperiment对象
SE_mRNA <-  GDCprepare(query = queryDown_mRNA)

###将SE对象简化为表格
dat_mRNA <- data.frame(assay(SE_mRNA))

###将ensembID转变为geneSymbol
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
geneID_list=select(org.Hs.eg.db,keys=k,columns = c("SYMBOL"), keytype="ENSEMBL")#获得EnsemblID与geneSymbol的对照表
ENSEMBL <- gsub("\\..*", "",  rownames(dat_mRNA))#将EnsemblID的版本号去掉
dat_mRNA <- cbind(dat_mRNA,ENSEMBL)
dat_mRNA <- merge(dat_mRNA,geneID_list)#横向合并两个数据框

###筛选出铁死亡相关基因

library("dplyr")

ferrGenes <- read.table("ferrGene_list.txt",col.names = "SYMBOL")



ferr_dat_mRNA <- dat_mRNA %>%
  filter(`SYMBOL` %in% ferrGenes$SYMBOL) %>%
  select(-c("ENSEMBL")) %>%
  select(`SYMBOL`,everything())
FerrCoLNames <- c("SYMBOL",rep("TP",time=16),rep("NT",time=16))
colnames(ferr_dat_mRNA) <- FerrCoLNames#修改列名
rownames(ferr_dat_mRNA) <- ferr_dat_mRNA$SYMBOL#修改行名
ferr_dat_mRNA <- ferr_dat_mRNA[,-1]#删除symbol列

###筛选完之只剩下224个基因,做差异表达分析

write.csv(ferr_dat_mRNA,file="ferr_dat_mRNA.csv") 
#### 1.差异表达矩阵
ferr_dat_mRNA <- read.csv("ferr_dat_mRNA.csv")

#### 2.样本分组矩阵

group_list <- c(rep("TP",time=16),rep("NT",time=16))
group_list <- factor(group_list,levels = c("TP","NT"))#分组信息因子
colData <- data.frame(row.names=colnames(ferr_dat_mRNA), 
                      group_list=group_list
                     )#得到样本分组矩阵

#### 3.构建dds
dds <- DESeqDataSetFromMatrix(countData = ferr_dat_mRNA,
                              colData = colData,
                              design = ~ group_list)

#### 4.做normalization
dds2 <- DESeq(dds)#标准化
suppressMessages(dds2)

#### 5.结果
res <-  results(dds2)
res <- res[order(res$padj),]#按照p值排序
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <-row.names(diff_gene_deseq2)#差异表达的基因，共有62个
resdata <-  merge(as.data.frame(res),
                  as.data.frame(counts(dds2,normalize=TRUE)),
                  by="row.names",sort=FALSE)
#### 6.输出为表格
write.csv(res,file="F:/R/Rdata/diff_res/diff_res.csv")#总的结果，排过序
write.table(diff_gene_deseq2,file="diff_res/diff_genes.txt",sep="\n",
            row.names = FALSE,col.names = FALSE,quote=FALSE)


#### 7.绘制热图
diff_ferr_mRNA <- ferr_dat_mRNA %>% 
  filter(row.names(ferr_dat_mRNA) %in% diff_gene_deseq2)#筛选出差异表达的基因

##### 列注释
annotation_col <- data.frame(
  type = rep(c("Cancer", "Noraml"), each = 16)
)
rownames(annotation_col) <- colnames(diff_ferr_mRNA)

##### 绘制
heatmap <- pheatmap(
                    diff_ferr_mRNA ,
                    scale = "row",# 归一化
                    color = rainbow(7),
                    annotation_col = annotation_col, #列注释
                    fontsize_row = 8,#行标签大小
                    fontsize_col = 8,#列标签大小
                    angle_col = 45#倾斜度
                   )

#### 8.绘制火山图
df <- as.data.frame(res)

df$group = ifelse(df$log2FoldChange>=1&df$padj<=0.05,"Up",
                  ifelse(df$log2FoldChange<=-1&df$padj<=0.05,"Down","Not sig"))

ggplot(df,aes(x=log2FoldChange,y = -log10(padj)))+
      geom_point(aes(color=group))+
      scale_color_manual(values = c("red","grey","blue"),limits = c('Up','Not sig',"Down"))+  
      theme_bw(base_size = 20)+
      theme(plot.title = element_text(size=30,hjust = 0.5))+
      coord_cartesian(xlim = c(-5,5),ylim = c(0,40)
      )


##2. 根据基因表达差异做生存分析


###阳性样本也就是癌症患者编号
samplesTP <- TCGAquery_SampleTypes(samplesDown, typesample = c("TP"))#1106个

samplesTP <- samplesTP[1:200]#太多了只要前200个
###所有阳性样本的基因表达矩阵
queryDown_all_mRNA <-GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          experimental.strategy = "RNA-Seq",
                          workflow.type = "STAR - Counts",
                          barcode = samplesTP
)

GDCdownload(queryDown_all_mRNA,
            method = "api",
            directory = "GDCdata_all_mRNA",
            files.per.chunk = 6
)  
####创建SummarizedExperiment对象
SE_all_mRNA <-  GDCprepare(query = queryDown_all_mRNA,
                           directory="GDCdata_all_mRNA"
                          )

####将SE对象简化为表格
dat_all_mRNA <- data.frame(assay(SE_all_mRNA))

####将ensembID转变为geneSymbol
#k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
#geneID_list=select(org.Hs.eg.db,keys=k,columns = c("SYMBOL"), keytype="ENSEMBL")#获得EnsemblID与geneSymbol的对照表
ENSEMBL <- gsub("\\..*", "",  rownames(dat_all_mRNA))#将EnsemblID的版本号去掉
dat_all_mRNA <- cbind(dat_all_mRNA,ENSEMBL)
dat_all_mRNA <- merge(dat_all_mRNA,geneID_list)#横向合并两个数据框

####筛选出铁死亡相关基因

library("dplyr")

ferrGenes <- read.table("ferrGene_list.txt",col.names = "SYMBOL")

ferr_dat_all_mRNA <- dat_all_mRNA %>%
  filter(`SYMBOL` %in% ferrGenes$SYMBOL) %>%
  select(-c("ENSEMBL")) %>%
  select(`SYMBOL`,everything())
#FerrCoLNames <- c("SYMBOL",rep("TP",time=16),rep("NT",time=16))
#colnames(ferr_dat_mRNA) <- FerrCoLNames#修改列名
rownames(ferr_dat_all_mRNA) <- ferr_dat_all_mRNA$SYMBOL#修改行名
ferr_dat_all_mRNA <- ferr_dat_all_mRNA[,-1]#删除symbol列

names(ferr_dat_all_mRNA) <- sapply(strsplit(names(ferr_dat_all_mRNA),'[.]'),
                                   function(x) paste0(x[1:3],collapse="-")
                                  )#修改样本名称，原是TCGA-DD-AAD5-01A-11R-A41C-07，
                                   #但是clin.BRCA数据中的样本名称是前12位，
                              

write.csv(ferr_dat_all_mRNA,"ferr_dat_all_mRNA.csv")#200个阳性样本，铁死亡相关基因的表达数据矩阵

ferr_dat_all_mRNA <- read.csv("ferr_dat_all_mRNA.csv")
rownames(ferr_dat_all_mRNA) <- ferr_dat_all_mRNA$X
ferr_dat_all_mRNA <- ferr_dat_all_mRNA[,-1]


genelist <- rownames(ferr_dat_all_mRNA)#获得要做差异表达的基因列表
### 临床数据处理

####所有的临床数据
clin_BRCA <- GDCquery_clinic("TCGA-BRCA", "clinical",save.csv = FALSE )
clin_BRCA <- clin_BRCA %>%
  #只要病人ID，结局，死亡时间，删失时间
  select(`submitter_id`,`vital_status`,`days_to_death`,`days_to_last_follow_up`)%>%
  filter(!is.na(`vital_status` ))#去除na值
write.csv(clin_BRCA,"clin_BRCA.csv")

clin_BRCA <- read.csv("clin_BRCA.csv")
#将status表示患者结局，1表示删失，2表示死亡
df2 <- clin_BRCA
df2 <- df2[,-1]
df2[df2$vital_status=='Dead',]$vital_status <- '2'
df2[df2$vital_status=='Alive',]$vital_status <- '1'
df2$vital_status <- as.numeric(df2$vital_status)

#time列为患者死亡或者删失时间
df2$time <- df2$days_to_death
df2$time[which(is.na(df2$time))] <- df2$days_to_last_follow_up[which(is.na(df2$time))]

write.csv(df2,"clin_BRCA_sim.csv")
genelist <- rownames(ferr_dat_all_mRNA)#获得用于生存分析的gene名称


####循环给每个基因做生存分析
#####dat用于存p值很小的基因
dat <- data.frame(
  symbol = '',
  pval=''
)
for(i in 1:length(genelist)){
  gene_symbol=genelist[i]
  dat_gene_exp <- ferr_dat_all_mRNA[gene_symbol,]#取出单个基因在阳性样本中的表达量
  dat_gene_exp<-as.data.frame(cbind(names(dat_gene_exp),t(dat_gene_exp)))#转置，行为样本编号，列为表达量
  colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)
  dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])#将表达量变成数字格式
  dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
  dat_clin$exp <- ''
  dat_clin[dat_clin[,6] >= mean(dat_clin[,6]),]$exp <- "H"#高表达基因为H
  dat_clin[dat_clin[,6] < mean(dat_clin[,6]),]$exp <- "L"#地表达基因为L
  #return(dat_clin)
  # 根据表达建模
  fit <- survfit(Surv(time, vital_status)~exp, data=dat_clin)
  #显示P value
  pval <- surv_pvalue(fit)$pval
  #return(pval)
  print(pval)
  dat <- rbind(dat,c(gene_symbol,pval))
}

####结果为：
dat

####选出p值很小的基因
surv_result <- dat
surv_result$pval <- as.numeric(surv_result$pval)#字符串转化为数字
surv_result <- surv_result[-1,]#去掉第一行
surv_result <- surv_result[surv_result$pval<0.05,]#筛选出p值小于0.05的基因

write.csv(surv_result,"surv_result.csv")

surv_result <- read.csv("surv_result.csv")

surv_result$symbol


driver_genes <- read.csv("driver_Hi.csv")
suppressor_genes <- read.csv("suppressor_Hi.csv")
marker_genes <- read.csv("marker_Hi.csv")

driver_genes <- driver_genes$symbol
suppressor_genes <- suppressor_genes$symbol
marker_genes <- marker_genes$symbol

dat <- data.frame(
  symbol = '',
  Type=''
)
for(i in 1:8){
  if(surv_result$symbol[i] %in% driver_genes) dat <- rbind(dat,c(surv_result$symbol[i],"driver_gene"))
  if(surv_result$symbol[i] %in% suppressor_genes) dat <- rbind(dat,c(surv_result$symbol[i],"suppressor_gene"))
  if(surv_result$symbol[i] %in% marker_genes) dat <- rbind(dat,c(surv_result$symbol[i],"marker_gene"))
}

write.csv(dat,"surv_res_genes")
surv_result$symbol[1] %in% driver_gene
#####针对p值小于0.05的基因绘制生存曲线

ferr_dat_all_mRNA <- read.csv("ferr_dat_all_mRNA.csv")
row.names(ferr_dat_all_mRNA) <- ferr_dat_all_mRNA$X#定义行名
ferr_dat_all_mRNA <- ferr_dat_all_mRNA[,-1]#删除第一列的genesymbol，因为已经放在行名中
names(ferr_dat_all_mRNA) <- sapply(strsplit(names(ferr_dat_all_mRNA),'[.]'),
                                   function(x) paste0(x[1:3],collapse="-")
                                  )#把列名的'.'，变成'-'

pval_gene <- read.csv("surv_result.csv")


for(i in 1:length(pval_gene$symbol)){
  gene_symbol=pval_gene$symbol[i]
  dat_gene_exp <- ferr_dat_all_mRNA[gene_symbol,]#取出单个基因在阳性样本中的表达量
  dat_gene_exp<-as.data.frame(cbind(names(dat_gene_exp),t(dat_gene_exp)))#转置，行为样本编号，列为表达量
  colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)
  dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])#将表达量变成数字格式
  dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
  dat_clin$exp <- ''
  dat_clin[dat_clin[,6] >= mean(dat_clin[,6]),]$exp <- "H"#高表达基因为H
  dat_clin[dat_clin[,6] < mean(dat_clin[,6]),]$exp <- "L"#地表达基因为L
  #return(dat_clin)
  # 根据表达建模
  fit <- survfit(Surv(time, vital_status)~exp, data=dat_clin)
  #显示P value
  File=paste("surv_res/mRNA/","surv_pict_",gene_symbol,".pdf",sep = '')
  pdf(file = File,width = 11.69,height = 8.27,
      onefile = F)
  ggsurvplot(fit,
             pval = TRUE, #P值
             surv.median.line = "hv",#中位生存时间
             #conf.int = TRUE,#置信区间
             #add.all = TRUE,#总患者生存曲线
             palette = "lancet",# 自定义调色板
             xlab = "Follow up time(d)", # 指定x轴标签
             legend = c(0.8,0.75), # 指定图例位置
             legend.title = gene_symbol, # 设置图例标题
             break.x.by = 1000 #设置x轴刻度间距
             )
  dev.off()
  
}


#wgcna


group_list <- c(rep("TP",time=0),rep("NT",time=200))
group_list <- factor(group_list,levels = c("TP","NT"))#分组信息因子
colData <- data.frame(row.names=colnames(ferr_dat_all_mRNA), 
                      group_list=group_list
)#得到样本分组矩阵

#### 3.构建dds
dds <- DESeqDataSetFromMatrix(countData = ferr_dat_all_mRNA,
                              colData = colData,
                              design = ~ 1)
dds_var <- log2(x+1)(dds)
#对表达矩阵进行标准化
ferr_dat_all_mRNA_nor <- assay(dds_var)

write.csv(ferr_dat_all_mRNA_nor,"ferr_dat_all_mRNA_nor.csv")
Expr <- t(ferr_dat_all_mRNA_nor)

#检查离群样本
sampleTree = hclust(dist(Expr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",xlim=10000) 
#去除离群值
Expr <- Expr[rownames(Expr)!="TCGA.A2.A3XV",]
Expr <- Expr[rownames(Expr)!="TCGA.A7.A5ZV",]

#选择软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(Expr, powerVector = powers,networkType = "signed hybrid", verbose = 5)
power = sft$powerEstimate

#绘制软阈值参考图
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")



net = blockwiseModules(Expr, power = power, maxBlockSize = 1000,
                       TOMType = "signed",minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson", maxPOutliers = 1,
                       loadTOMs=TRUE,saveTOMFileBase = "data.tom",
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#模块之间相关性
library(stringr)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

#可视化基因网络
TOM = TOMsimilarityFromExpr(Expr, power=power)
probes = colnames(Expr)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("edges.txt", sep=""),
                               nodeFile = paste("nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.005,
                               nodeNames = probes, nodeAttr = moduleColors)