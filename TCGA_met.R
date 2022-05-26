BiocManager::install("sesameData")
BiocManager::install("sesame")
#library(BiocManager)
library("TCGAbiolinks")
library("survival")
library("survminer")
library("tidyverse")
#library("ggplot2")
#library("ggpubr")
#library(org.Hs.eg.db)
#library(AnnotationDbi)
#library("DESeq2")
#library("pheatmap")
#library("sesameData")
#library("sesame")


query_met <- GDCquery(project = "TCGA-BRCA",
                      data.category = "DNA Methylation",
                      data.type = "Methylation Beta Value",
                      experimental.strategy = "Methylation Array",
                      platform = "Illumina Human Methylation 450"
                     )

samplesDown_met <- getResults(query_met,cols=c("cases"))#共有895个样本

samplesDown_met_TP <- TCGAquery_SampleTypes(barcode = samplesDown_met,
                                            typesample = "TP")#阳性样本共有793个
samplesDown_met_TP <- samplesDown_met_TP[1:100]#太多了下载太慢，只取前100个病例

queryDown_met <- GDCquery(project = "TCGA-BRCA",
                          data.category = "DNA Methylation",
                          data.type = "Methylation Beta Value",
                          experimental.strategy = "Methylation Array",
                          platform = "Illumina Human Methylation 450",
                          barcode = samplesDown_met_TP
                          )

GDCdownload(queryDown_met,
            method = "api",
            directory = "GDCdata_met",
            files.per.chunk = 6
           ) 

SE_met <- GDCprepare(query = queryDown_met,directory = "GDCdata_met",
                     save = TRUE, save.filename ="SE_met.rda")

dat_met <- data.frame(assay(SE_met))

write.csv(dat_met,"dat_met.csv")

#dat_met <- read.csv("dat_met.csv")
####得到甲基化矩阵，行为甲基化芯片探针，列为病人样本
dat_met <- read.csv("dat_met.csv")
colnames(dat_met)[1] <- "IlmnID"
####对探针编号做注释
dat_ref <- read.csv(file = "GPL13534_HumanMethylation450_15017482_v.1.1.csv")

dat_met_ref <- dat_ref%>%
  filter(`UCSC_RefGene_Name`!='')%>%
  select(`IlmnID`,`UCSC_RefGene_Name`)

dat_met_ref$UCSC_RefGene_Name <- sapply(strsplit(dat_met_ref$UCSC_RefGene_Name,';'),
                                        function(x) paste(x[1:1])
                                        ) #多个symbol只选前一个
dat_met <- merge(dat_met,dat_met_ref,by = "IlmnID")#通过探针序列号合并两个矩阵，每个探针都对应一个基因id

ferrGenes <- read.table("ferrGene_list.txt")
ferrGenes <- ferrGenes[,1]#铁死亡相关基因
dat_met[is.na(dat_met)] <- 0#去掉na值
ferr_dat_met <- dat_met %>%
  #基因id列提前，删除探针id
  select(SYMBOL=`UCSC_RefGene_Name`,everything())%>%
  select(-c("IlmnID"))%>%
  #筛选出铁死亡相关基因
  filter(`SYMBOL` %in% ferrGenes)%>%
  #将相同的基因名合并并且取甲基化值的中值
  group_by(`SYMBOL`)%>%
  summarise_each(funs(mean))

ferr_dat_met <- data.frame(ferr_dat_met)
rownames(ferr_dat_met) <- ferr_dat_met$SYMBOL#修改行名
ferr_dat_met <- ferr_dat_met[,-1]#删除symbol列

names(ferr_dat_met) <- sapply(strsplit(names(ferr_dat_met),'[.]'),
                                   function(x) paste0(x[1:3],collapse="-")
                              )#修改样本名称，原是TCGA-DD-AAD5-01A-11R-A41C-07，
                              #但是clin.BRCA数据中的样本名称是前12位，
write.csv(ferr_dat_met,"ferr_dat_met.csv")

### 生存分析
#读取文件则要进行下面操作
ferr_dat_met <- read.csv("ferr_dat_met.csv")
rownames(ferr_dat_met) <- ferr_dat_met$X
ferr_dat_met <- ferr_dat_met[,-1]
names(ferr_dat_met) <- sapply(strsplit(names(ferr_dat_met),'[.]'),
                              function(x) paste0(x[1:3],collapse="-")
                              )#修改样本名称，原是TCGA-DD-AAD5-01A-11R-A41C-07，
                                #但是clin.BRCA数据中的样本名称是前12位，
  
df2 <- read.csv("clin_BRCA_sim.csv")
df2 <- df2[,-1]#删除第一列
genelist <- rownames(ferr_dat_met)#获得用于生存分析的gene名称
####循环给每个基因做生存分析
#####dat用于存p值很小的基因
dat <- data.frame(
  symbol = '',
  pval=''
)
for(i in 1:length(genelist)){
  gene_symbol=genelist[i]
  dat_gene_exp <- ferr_dat_met[gene_symbol,]#取出单个基因在阳性样本中的甲基化值
  dat_gene_exp<-as.data.frame(cbind(names(dat_gene_exp),t(dat_gene_exp)))#转置，行为样本编号，列为甲基化值
  colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)#修改列名与临床数据中的列名相同
  dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])#将甲基化值变成数字格式
  dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
  dat_clin$exp <- ''
  A <- try(dat_clin[dat_clin[,6] >= 0.6,]$exp <- "A",silent = T)
  N <- try(dat_clin[dat_clin[,6] < 0.2,]$exp <- "N",silent = T)
  S <- try(dat_clin[(dat_clin[,6] >=0.2) & (dat_clin[,6] <0.6),]$exp <- "S",silent = T)
  #return(dat_clin)
  # 根据表达建模
  fit <- survfit(Surv(time, vital_status)~exp, data=dat_clin)
  #ggsurvplot(fit)
  #显示P value
  Error <- tryCatch(pval <- surv_pvalue(fit)$pval,silent = T)
  #return(pval)
  #print(pval)
  if(!is.na(Error))
  dat <- rbind(dat,c(gene_symbol,pval))
}
####结果为：
dat

####选出p值很小的基因
surv_result <- dat
surv_result$pval <- as.numeric(surv_result$pval)#字符串转化为数字
surv_result <- surv_result[-1,]#去掉第一行
surv_result <- surv_result[surv_result$pval<0.05,]#筛选出p值小于0.05的基因

write.csv(surv_result,"surv_met_result.csv")


#画图


gene_symbol=surv_result$symbol[14]
dat_gene_exp <- ferr_dat_met[gene_symbol,]#取出单个基因在阳性样本中的甲基化值
dat_gene_exp<-as.data.frame(cbind(names(dat_gene_exp),t(dat_gene_exp)))#转置，行为样本编号，列为甲基化值
colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)#修改列名与临床数据中的列名相同
dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])#将甲基化值变成数字格式
dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
dat_clin$exp <- ''
A <- try(dat_clin[dat_clin[,6] >= 0.6,]$exp <- "A",silent = T)
N <- try(dat_clin[dat_clin[,6] < 0.2,]$exp <- "N",silent = T)
S <- try(dat_clin[(dat_clin[,6] >=0.2) & (dat_clin[,6] <0.6),]$exp <- "S",silent = T)
#return(dat_clin)
# 根据表达建模
fit <- survfit(Surv(time, vital_status)~exp, data=dat_clin)
File=paste("surv_res/met/","surv_pict_",gene_symbol,".pdf",sep = '')
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


surv_result <- read.csv("surv_met_result.csv")

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
for(i in 1:14){
  if(surv_result$symbol[i] %in% driver_genes) dat <- rbind(dat,c(surv_result$symbol[i],"driver_gene"))
  if(surv_result$symbol[i] %in% suppressor_genes) dat <- rbind(dat,c(surv_result$symbol[i],"suppressor_gene"))
  if(surv_result$symbol[i] %in% marker_genes) dat <- rbind(dat,c(surv_result$symbol[i],"marker_gene"))
}

write.csv(dat,"surv_res_met_genes")
