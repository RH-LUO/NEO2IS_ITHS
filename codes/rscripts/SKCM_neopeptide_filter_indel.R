#!/usr/local/bin/Rscript
myPaths <- .libPaths()  
new <- c('/public/home/lorihan/miniconda3/lib/R/library','/usr/local/lib64/R/library','/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')#c('/public/home/lorihan/miniconda3/lib/R/library')#
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
#####
library(maftools)## 
library(data.table)
library(stringr)
load("/public/home/lorihan/lrh/SKCM/SKCM_202203_result.Rdata")
setwd("/public/home/lorihan/lrh/SKCM/SNV/mutect2_results/filter.vcfs/vep.vcfs/")
path_to_maf = "vep_merge.maf"
maf.passed = data.table::fread(file = path_to_maf, stringsAsFactors = FALSE)
head(maf.passed,5)
dim(maf.passed)##77493   133
unique(maf.passed$Tumor_Sample_Barcode)
path <- "/public/home/lorihan/lrh/SKCM/SNV/mutect2_results/filter.vcfs/vep.vcfs"
maf_files <- list.files(path)
maf_files <- maf_files[grep("_vep.maf",maf_files)]
filePath <- sapply(maf_files, function(x){ 
  paste(path,x,sep='/')})  
maf_data <- lapply(filePath, function(x){
  if (file.info(x)$size!=0) {
    fread(x,header = T,sep = "\t",fill = F)}})  
n <- 1:length(maf_data)
names <- 1:nrow(maf.passed)
for (i in 1:length(maf_data)) {
  n[i] <- nrow(maf_data[[i]])
}
#for (i in 1:length(maf_data)) {
# maf_data[[i]]$Tumor_Sample_Barcode <- rep(names(maf_data)[i],nrow(maf_data[[i]]))}
names <- rep(substr(names(maf_data),1,16),n)
maf.passed$Tumor_Sample_Barcode <- names
#TCGA_SKCM_VCF$Matched[match(names,TCGA_SKCM_VCF$Sample)]
maf.passed$Matched_Norm_Sample_Barcode <- TCGA_SKCM_VCF$Matched[match(names,TCGA_SKCM_VCF$Sample)]


table(is.na(maf.passed$t_alt_count))#FALSE
table(is.na(maf.passed$t_depth))#FALSE
#maf.passed[!is.na(maf.passed$t_alt_count),]$Tumor_Sample_Barcode
maf.passed$Score <- round(maf.passed$t_alt_count/maf.passed$t_depth,3)
table(maf.passed$Score>0.02)#3124 
table(maf.passed$FILTER)

#maf.annovar <- read.table("/public/home/lorihan/lrh/NSCLC/AnnoFilter/tmp.hg38_multianno.txt",header = T,sep = "\t")
maf.annovar <- read.table("tmp.hg38_multianno.txt",header = T,sep = "\t")
dim(maf.annovar)
setdiff(maf.annovar$Start,maf.passed$Start_Position)
colnames(maf.annovar)
maf.annovar <- maf.annovar[-1,]
dim(maf.annovar)
setdiff(maf.annovar$Start,maf.passed$Start_Position)
maf.annovar <- cbind(maf.passed,maf.annovar[,-(1:5)])
colnames(maf.annovar)
unique(maf.annovar$Tumor_Sample_Barcode)#5

grep("FILTER",colnames(maf.annovar))
table(maf.annovar$FILTER!="PASS")
table(maf.annovar$FILTER)#3024
unique(maf.annovar$Tumor_Sample_Barcode)#5

table(na.omit(maf.annovar$X1000g2015aug_all)>0.05)#3-->1
filter.1000g <- which(maf.annovar$X1000g2015aug_all > 0.05)
maf <- maf.annovar[-filter.1000g,]
unique(maf$Tumor_Sample_Barcode)#5

table(na.omit(maf.annovar$ExAC_ALL)>0.05)#8 
table(na.omit(maf$ExAC_ALL)>0.05)#7
filter.ExAC <- which(maf$ExAC_ALL > 0.05)
#maf <- maf[-filter.ExAC,]
unique(maf$Tumor_Sample_Barcode)#5

table(na.omit(maf.annovar$gnomAD_genome_ALL)>0.05)#8
table(na.omit(maf$gnomAD_genome_ALL)>0.05)#1
filter.gnomeAD <- which(maf$gnomAD_genome_ALL > 0.05)
#maf <- maf[-filter.gnomeAD,]
unique(maf$Tumor_Sample_Barcode)#5

filter.ExAC <- which(maf.annovar$ExAC_ALL > 0.05)
filter.gnomeAD <- which(maf.annovar$gnomAD_genome_ALL > 0.05)
filter.maf <- unique(c(filter.1000g,filter.ExAC,filter.gnomeAD))
filter.maf <- unique(c(filter.1000g,filter.ExAC))
length(unique(filter.maf))#-->1
table(!is.na(maf.annovar$genomicSuperDups))#T:615
table(maf.annovar$rmsk)
table(!is.na(maf.annovar$rmsk))#1518
dim(maf.annovar)#3258  152

maf.annovar <- maf.annovar[-filter.maf,]
dim(maf.annovar)#3014  152 -->3257
dim(maf)#3013  152
maf.annovar$Score <- round(maf.annovar$t_alt_count/maf.annovar$t_depth,3)
#maf$Score <- round(maf$t_alt_count/maf$t_depth,3)
table(maf.annovar$Score>0.02)
table(maf$Score>0.02)

table(maf.annovar$t_depth<5)
table(maf.annovar$t_alt_count<5)
min(maf.annovar$t_depth)
maf.annovar<-maf.annovar[maf.annovar$Score>=0.02,]#0.4%
maf <- data.frame(maf.annovar)
maf<-maf[maf$Score>=0.02,]#0.4%
unique(maf$Tumor_Sample_Barcode)#5

dim(maf.annovar)#77288   152
table(maf.annovar$t_depth>2)

#######
maf <- data.frame(maf.annovar)
unique(maf$Tumor_Sample_Barcode)#"TCGA-44-2661-01A"
table(!is.na(maf$genomicSuperDups))#612
table(!is.na(maf$rmsk))#18423

table(maf$Variant_Type)
table(maf$Variant_Type !="SNP")
table(maf[!is.na(maf$rmsk),]$Variant_Type)
table(maf$Score>=0.004)
table(maf[!is.na(maf$genomicSuperDups),]$Score>0.04)
table(maf[!is.na(maf$rmsk) & maf$Variant_Type !="SNP",]$Score<0.1)
table(maf[!is.na(maf$genomicSuperDups) & maf$Variant_Type !="SNP",]$Score<0.1)
#F:1144  T:1809 
####
filter.rmsk <- which(!is.na(maf$rmsk) & maf$Variant_Type !="SNP")
filter.dups <- which(!is.na(maf$genomicSuperDups) & maf$Variant_Type !="SNP")

filter.repeat <- union(filter.rmsk,filter.dups)

maf.repeat <- maf[filter.repeat,]
max(maf.repeat$Score)
dim(maf)
dim(maf[-filter.repeat,])
maf <- maf[-filter.repeat,]
unique(maf$Tumor_Sample_Barcode)#5
dim(maf)#30575   152
unique(maf$Tumor_Sample_Barcode)
#maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode,1,15)
maf.SKCM <- maf

library("maftools")## 
library(data.table)
d1<-read.maf(maf=maf)#,clinicalData=clin.all)
#unique(d1@data$Hugo_Symbol)
unique(d1@data$Tumor_Sample_Barcode)
#2944
d2.SKCM <- d1@data
d2.SKCM$Exon_Number <- sapply(str_split(d2.SKCM$Exon_Number,"[/]"),"[",1)
dim(d2.SKCM)
unique(d2.SKCM$Hugo_Symbol)[1:20]#38
min(d2.SKCM$Score)
table(d2.SKCM$Tumor_Sample_Barcode)
unique(d2.SKCM$Tumor_Sample_Barcode)

#######
table(d2.SKCM$FILTER)#20316
table(d2.SKCM$Variant_Classification)### 
#mv test ../../
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/fastaFiles/test.fa" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
data <- lapply(filePath, function(x){
  fread(x, sep = "",fill = T,header = F)})  ## ， list

####
data.new <- data
names(data.new)
name.new <- sapply(str_split(names(data.new),"[.]"),"[",1)
table(name.new %in% d2.SKCM$Tumor_Sample_Barcode)
name.new[!name.new %in% d2.SKCM$Tumor_Sample_Barcode]
data.new <- data.new[name.new %in% d2.SKCM$Tumor_Sample_Barcode]
length(data.new)#330
#######
table(d2.SKCM$Tumor_Sample_Barcode %in% name.new)
d2.test <- d2.SKCM[d2.SKCM$Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", 
                                                     "In_Frame_Del", "Missense_Mutation,", 
                                                     "Nonsense_Mutation") & d2.SKCM$Variant_Type %in% c("INS","DEL"),]

d2.test$Tumor_Sample_Barcode[!d2.test$Tumor_Sample_Barcode %in% name.new]
#0
library('org.Hs.eg.db')
library(clusterProfiler)
library(AnnotationDbi)
load("/public/home/lorihan/lrh/NSCLC/gene_info.RData")
#library(enrichplot)
## ， ， for 。
for(i in 1:length(data.new)){
  data.new[[i]] <- gsub(";"," ",data.new[[i]])
  test <-str_split(data.new[[i]],">")
  #length(test[[1]])#89
  test <- test[[1]][-grep("WILDTYPE",test[[1]])]
  test <- test[-match("",test)]
  test.p <- sapply(str_split(test,"[ ]"),'[',4)
  ######INDEL
  test.p.2 <- sapply(str_split(test.p,"[.]"),'[',2)
  test.p.2
  test.p.2 <- sapply(str_split(test.p.2,"[fs*]"),'[',1)
  test.p.2
  test.p.2 <- sapply(str_split(test.p.2,"[_]"),'[',1)
  test.p.2
  d2.test <- d2.SKCM[d2.SKCM$Tumor_Sample_Barcode %in% sapply(str_split(names(data.new)[[i]],"[.]"),"[",1),]
  #####
  d2.test <- d2.test[d2.test$Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", 
                                                           "In_Frame_Del", "Missense_Mutation,", 
                                                           "Nonsense_Mutation") & d2.test$Variant_Type %in% c("INS","DEL"),]
  d2.HGVSp_Short <- d2.test$HGVSp_Short
  test.logi.1 <- as.logical(test.p %in% d2.test$HGVSp_Short)
  test.logi.2 <- test.logi.1
  table(test.logi.2)
  test.logi.1[which(!test.logi.1)]
  ####indel
  for (j in which(!test.logi.2)) {
    test.logi.2[j] <- length(unique(grepl(test.p.2[j],as.vector(d2.HGVSp_Short))))==2
  }
  test.new <- test[test.logi.2]
  test.new[1:5]
  length(test.new)
  
  test.g <- sapply(str_split(test.new,"[ ]"),'[',2)
  test.p.3 <- sapply(str_split(test.new,"[ ]"),'[',4)
  test.p.3 <- sapply(str_split(test.p.3,"[.]"),'[',2)
  test.p.3
  test.p.3 <- sapply(str_split(test.p.3,"[fs*]"),'[',1)
  test.p.3
  test.p.3 <- sapply(str_split(test.p.3,"[_]"),'[',1)
  test.p.3
  ###indel
  d2.p <- d2.HGVSp_Short
  d2.p <- d2.p[1:length(test.new)]
  if(length(test.new)==0) {
    d2.p <- d2.p
  } else {
    for (h in 1:length(test.new)) {
      d2.p[h] <- d2.HGVSp_Short[grepl(test.p.3[h],as.vector(d2.HGVSp_Short))]
    }  }
  
  d2.p
  test.p.3
  #######
  #gene.df <- bitr(test.g, fromType = "REFSEQ",
  #               toType = "SYMBOL",
  #              OrgDb = "org.Hs.eg.db",drop=TRUE)
  # test.g.tr <-  gene.df$SYMBOL[match(test.g,gene.df$REFSEQ)]
  test.g.tr <-  gene_info$symbols[match(test.g,gene_info$geneRefSeqs)]
  d2.test[match(d2.p[is.na(test.g.tr)],d2.test$HGVSp_Short),1:50]
  d2.p[is.na(test.g.tr)]
  test.p.3[is.na(test.g.tr)]
  test.g.tr[is.na(test.g.tr)] <- d2.test$Hugo_Symbol[match(d2.p[is.na(test.g.tr)],d2.test$HGVSp_Short)]
  
  test.split <- str_split(test.new,"[ ]")
  test.new <- test.split
  if(length(test.new)==0) {
    test.new <- test.new
  } else {
    for(k in 1:length(test.split)){
      test.new[[k]][1] <- test.g.tr[k]
      test.new[[k]][2] <- test.new[[k]][4]
      test.new[[k]] <- paste(">",paste(test.new[[k]][-4],collapse = " "),sep = "")
    }    }  
  
  length(test)
  test.new <- as.character(test.new)
  data.new[[i]] <- test.new
}

outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/test.fasta"#path## 
out_fileName <- sapply(names(data),function(x){
  paste(x, "sta", sep='')}) ##csv 
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ## 
out_filePath
## 
for(i in 1:length(data.new)){
  write.table(data.new[[i]], file=out_filePath[i],col.names = F, quote = F,row.names = F)
}
