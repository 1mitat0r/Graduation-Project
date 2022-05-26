#下载包
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
library("tidyverse")

#读取下载好的三个类别的基因数据
driver <- read.csv("driver.csv")
suppressor <- read.csv("suppressor.csv")
marker <- read.csv("marker.csv")

#使用dplyr筛选出高质量基因

## 1.driver基因
driver.Hi <- driver %>%
  ###物种限制为人类
  filter(grepl('Human',`testin`)) %>%
  
  ###可信度限制为Validated
  filter(`confidence` == 'Validated') %>%
  
  ##去除没有uniprot编号的基因，因为该基因可能对应的不是蛋白质
  filter(`uniprotid` != '_NA_') %>%
  
  ###只需要基因名称，蛋白质全称和uniprot序列号
  select(`symbol`,`name`,`testin`,`proteinfullname`,`uniprotac`,`uniprotid`)

## 2.suppressor基因
suppressor.Hi <- suppressor %>%
  ###物种限制为人类
  filter(grepl('Human',`testin`)) %>%
  
  ###可信度限制为Validated
  filter(`confidence` == 'Validated') %>%
  
  ##去除没有uniprot编号的基因，因为该基因可能对应的不是蛋白质
  filter(`uniprotid` != '_NA_') %>%
  
  ###只需要基因名称，蛋白质全称和uniprot序列号
  select(`symbol`,`name`,`testin`,`proteinfullname`,`uniprotac`,`uniprotid`)

## 3.marker基因
marker.Hi <- marker %>%
  ##物种限制为人类
  filter(grepl('Human',`testin`)) %>%
  
  ##不限制可信度，因为marker基因中Validated的基因太少
  ##filter(`confidence` == 'Validated') %>%
  
  ##去除没有uniprot编号的基因，因为该基因可能对应的不是蛋白质
  filter(`uniprotid` != '_NA_') %>%
  
  ##只需要基因名称，蛋白质全称和uniprot序列号
  select(`symbol`,`name`,`testin`,`proteinfullname`,`uniprotac`,`uniprotid`)

#将筛选过的基因导出为csv格式
write.csv(driver.Hi,"driver_Hi.csv")
write.csv(suppressor.Hi,"suppressor_Hi.csv")
write.csv(marker.Hi,"marker_Hi.csv")

#只需要基因名称用于string数据库构建网络
write.table(driver.Hi$uniprotac,"driver_ID.txt",sep="\n",row.names = FALSE,
            col.names = FALSE,quote=FALSE)
write.table(suppressor.Hi$uniprotac,"suppressor_ID.txt",sep="\n",row.names = FALSE,
            col.names = FALSE,quote=FALSE)
write.table(marker.Hi$uniprotac,"marker_ID.txt",sep="\n",row.names = FALSE,
            col.names = FALSE,quote=FALSE)

#将所有基因名称合成一个表
ferrGene_list <- c(driver.Hi$symbol,suppressor.Hi$symbol,marker.Hi$symbol)
write.table(ferrGene_list,"ferrGene_list.txt",sep="\n",row.names = FALSE,
            col.names = FALSE,quote=FALSE)
ferrGene_list <- read.table("ferrGene_list.txt")
