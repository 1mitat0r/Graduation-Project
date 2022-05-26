library(BiocManager)
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


query_mut <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Simple Nucleotide Variation",
                      access = "Open"
)

samplesDown_mut <- getResults(query_mut,cols=c("cases"))#共有988个样本

samplesDown_mut_TP <- TCGAquery_SampleTypes(barcode = samplesDown_mut,
                                            typesample = "TP")#阳性样本共有982个
#samplesDown_mut_TP <- samplesDown_mut_TP[1:400]#太多了下载太慢，只取前400个病例

queryDown_mut <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Simple Nucleotide Variation",
                          data.type = "Masked Somatic Mutation",
                          access = "Open",
                          barcode = samplesDown_mut_TP
                         )

GDCdownload(queryDown_mut,
            method = "api",
            directory = "GDCdata_mut",
            files.per.chunk = 6
) 

SE_mut <-  GDCprepare(query = queryDown_mut,directory = "GDCdata_mut",
                      save = TRUE, save.filename ="SE_mut.rda")

dat_mut <- data.frame(SE_mut)

write.csv(dat_mut,"dat_mut.csv")

dat_mut <- read.csv("dat_mut.csv")

##处理变异数据
dat_mut_1 <- dat_mut %>% 
  select(SYMBOL = `Hugo_Symbol`,submitter_id = `Tumor_Sample_Barcode`)%>%#只要基因名和样本编号
  mutate(mix = paste(SYMBOL,submitter_id,sep = '-'),Index="1")%>%
  filter(!duplicated(mix))%>%#除去相同样本中的重复基因
  select(SYMBOL,submitter_id,Index)

#每一行为样本编号，每一列为基因，1为变异，na为没有变异  
dat_mut_2 <- spread(dat_mut_1,SYMBOL,Index)#将长数据变为宽数据，数据扁平化


write.csv(dat_mut_2,"dat_mut_2.csv")

dat_mut_2 <- read.csv("dat_mut_2.csv")
#铁死亡相关基因
ferrGenes <- read.table("ferrGene_list.txt")
ferrGenes <- ferrGenes[,1]#铁死亡相关基因

x <- c(names(dat_mut_2)[2:length(names(dat_mut_2))] %in% ferrGenes)
submitter_id <- dat_mut_2$submitter_id
dat_mut_2 <- dat_mut_2[,c(-1,-2)]
ferr_dat_mut <- dat_mut_2[, names(dat_mut_2) %in% ferrGenes ]#筛选出铁死亡相关基因




ferr_dat_mut$submitter_id<- sapply(strsplit(submitter_id,'-'),
                                   function(x) paste0(x[1:3],collapse="-")
                                   )#修改样本名称，原是TCGA-DD-AAD5-01A-11R-A41C-07，
                                    #但是clin.BRCA数据中的样本名称是前12位，

ferr_dat_mut<- filter(ferr_dat_mut,!duplicated(ferr_dat_mut$submitter_id))#去除重复
rownames(ferr_dat_mut) <- ferr_dat_mut$submitter_id#重新命名列名
ferr_dat_mut <- subset(ferr_dat_mut,select = -c(submitter_id))#去除不需要的列

write.csv(ferr_dat_mut,"ferr_dat_mut.csv")
ferr_dat_mut <- read.csv("ferr_dat_mut.csv")
rownames(ferr_dat_mut) <- ferr_dat_mut$X
ferr_dat_mut <- ferr_dat_mut[,-1]
#获得临床数据
df2 <- read.csv("clin_BRCA_sim.csv")
df2 <- df2[,-1]#删除第一列
genelist <- names(ferr_dat_mut)#获得用于生存分析的gene名称

dat <- data.frame(
  symbol = '',
  pval=''
)
for(i in 1:length(genelist)){
  gene_symbol=genelist[i]
  dat_gene_exp <- ferr_dat_mut[,gene_symbol]#取出单个基因在阳性样本是否变异的数据
  dat_gene_exp<-as.data.frame(cbind(rownames(ferr_dat_mut),dat_gene_exp))#行为样本编号，列为是否变异
  colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)#修改列名与临床数据中的列名相同
  dat_gene_exp[is.na(dat_gene_exp[,2]),2] <- "0"
  dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])
  dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
  dat_clin$index <- ''
  try(dat_clin[dat_clin[,6] == 0,]$index <- "nor",silent = T)
  try(dat_clin[dat_clin[,6] == 1,]$index <- "mut",silent = T)
  #return(dat_clin)
  # 根据表达建模
  fit <- survfit(Surv(time, vital_status)~index, data=dat_clin)
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

write.csv(surv_result,"surv_mut_result.csv")
surv_result <- read.csv("surv_mut_result.csv")
genelist <- surv_result$symbol

#画图：
gene_symbol=genelist[14]
dat_gene_exp <- ferr_dat_mut[,gene_symbol]#取出单个基因在阳性样本是否变异的数据
dat_gene_exp<-as.data.frame(cbind(rownames(ferr_dat_mut),dat_gene_exp))#行为样本编号，列为是否变异
colnames(dat_gene_exp) <- c("submitter_id",gene_symbol)#修改列名与临床数据中的列名相同
dat_gene_exp[is.na(dat_gene_exp[,2]),2] <- "0"
dat_gene_exp[,gene_symbol]<- as.numeric(dat_gene_exp[,gene_symbol])
dat_clin <- merge(df2,dat_gene_exp,by="submitter_id")
dat_clin$index <- ''
try(dat_clin[dat_clin[,6] == 0,]$index <- "nor",silent = T)
try(dat_clin[dat_clin[,6] == 1,]$index <- "mut",silent = T)
#return(dat_clin)
# 根据变异建模
fit <- survfit(Surv(time, vital_status)~index, data=dat_clin)
File=paste("surv_res/mut/","surv_pict_",gene_symbol,".pdf",sep = '')
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


#多个基因的生存分析

library("tidyverse")
driver_genes <- read.csv("driver_Hi.csv")
suppressor_genes <- read.csv("suppressor_Hi.csv")
marker_genes <- read.csv("marker_Hi.csv")

driver_gene <- genelist[which(genelist %in% driver_genes$symbol)]
suppressor_gene <- genelist[which(genelist %in% suppressor_genes$symbol)]
marker_gene <- genelist[which(genelist %in% marker_genes$symbol)]

driver_dat_mut <- ferr_dat_mut[,names(ferr_dat_mut) %in% driver_gene]
suppressor_dat_mut <- ferr_dat_mut[,names(ferr_dat_mut) %in% suppressor_gene]
marker_dat_mut <- ferr_dat_mut[,names(ferr_dat_mut) %in% marker_gene]

driver_submitter_id <- rownames(driver_dat_mut)
suppressor_submitter_id <- rownames(driver_dat_mut)
marker_submitter_id <- rownames(driver_dat_mut)

driver_dat_mut[is.na(driver_dat_mut)] <- 0
suppressor_dat_mut[is.na(suppressor_dat_mut)] <- 0
marker_dat_mut[is.na(marker_dat_mut)] <- 0


driver_dat_mut <- lapply(driver_dat_mut,as.numeric)
suppressor_dat_mut <- lapply(suppressor_dat_mut,as.numeric)
marker_dat_mut <- lapply(marker_dat_mut,as.numeric)


driver_dat_mut <- as.tibble(driver_dat_mut)
suppressor_dat_mut <- as.tibble(suppressor_dat_mut)
marker_dat_mut <- as.tibble(marker_dat_mut)

driver_dat_mut <- driver_dat_mut %>%
  mutate(rowsum = rowSums(.[1:10]))
suppressor_dat_mut <- suppressor_dat_mut %>%
  mutate(rowsum = rowSums(.[1:length(suppressor_dat_mut)]))
marker_dat_mut <- marker_dat_mut %>%
  mutate(rowsum = rowSums(.[1:length(marker_dat_mut)]))

driver_dat_mut$submitter_id <- driver_submitter_id
suppressor_dat_mut$submitter_id <- suppressor_submitter_id
marker_dat_mut$submitter_id <- marker_submitter_id

dat_clin <- merge(df2,marker_dat_mut,by="submitter_id")
dat_clin$index=''
dat_clin[dat_clin[,"rowsum"] > 0,]$index <- "mut"
dat_clin[dat_clin[,"rowsum"] ==0,]$index <- "nor"

#return(dat_clin)
# 根据表达建模
fit <- survfit(Surv(time, vital_status)~index, data=dat_clin)
ggsurvplot(fit)
#显示P value
pval <- surv_pvalue(fit)$pval
#return(pval)
print(pval)

File=paste("surv_res/met/","surv_pict_","marker",".pdf",sep = '')
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

surv_result <- read.csv("surv_mut_result.csv")

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
for(i in 1:24){
  if(surv_result$symbol[i] %in% driver_genes) dat <- rbind(dat,c(surv_result$symbol[i],"driver_gene"))
  if(surv_result$symbol[i] %in% suppressor_genes) dat <- rbind(dat,c(surv_result$symbol[i],"suppressor_gene"))
  if(surv_result$symbol[i] %in% marker_genes) dat <- rbind(dat,c(surv_result$symbol[i],"marker_gene"))
}

write.csv(dat,"surv_res_met_genes")

for(i in 1:23){
  if(surv_result$symbol[i] %in% suppressor_genes) print(i)
}