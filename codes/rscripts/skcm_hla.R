#!/usr/local/bin/Rscript
setwd("/public/home/lorihan/lrh/SKCM")
.libPaths()  
myPaths <- .libPaths()  
#new <- '/public/apps/R-4.1.0/library'
new <- '/public/home/lorihan/miniconda3/lib/R/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
library(stringr)
load("/public/home/lorihan/lrh/SKCM/SKCM_202203_result.Rdata")
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% clin.all$Tumor_Sample_Barcode)
HLA.type[1:5,1:5]
dim(HLA.type)#469
length(unique(TCGA_SKCM_VCF$Patients))#469
table(TCGA_SKCM_VCF$Sample %in% HLA.type$patientBarcode)#469
output <- 1:nrow(HLA.type)
for(i in 1:nrow(HLA.type)){
  output[i] <- paste(HLA.type[i,-1],collapse = ",")
  output[i] <- gsub("NA,","",output[i])
}

outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/hla" ## 
out_fileName <- sapply(HLA.type[,1],function(x){
  paste(x, ".hla", sep='')}) ##csv 
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ## 
## 
for(i in 1:nrow(HLA.type)){
  write.table(output[[i]], file=out_filePath[i], quote = F,row.names = F,col.names = F)
}