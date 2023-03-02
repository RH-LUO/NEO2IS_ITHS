#!/usr/local/bin/Rscript
myPaths <- .libPaths()  
#new <- c('/usr/local/lib64/R/library','/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')#c('/public/home/lorihan/miniconda3/lib/R/library')#
new <- '/public/home/lorihan/miniconda3/lib/R/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
library(seqinr)
#mv test ../../
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/test.fasta" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
faPath <- filePath[grep("snp.fa",filePath)]
#faPath <- faPath[-grep("snp.fasta",faPath)]

path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/pos" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath

bedPath <- filePath[grep("snp.bed",filePath)]

snp <- lapply(faPath, function(x){
  read.table(x,sep = ",",fill = T,header = F)})  ## ， list
snp.fasta <- lapply(faPath, function(x){
  read.fasta(x)})  ## ， list
snp.bed <- lapply(bedPath, function(x){
  read.table(x,fill = T,header = F)})
#snp.bed.1 <- lapply(snp.bed, function(x){x[!(is.na(x$V4)),]})

length(snp.fasta[[1]])
n1 <- snp.fasta[[1]]
#count(n1,1)
#annotation <- 1:length(snp.fasta[[1]])
n <- snp.fasta
remove(annotation)
for(j in 1:length(snp)){
  annotation <- getAnnot(snp.fasta[[j]])
  for(i in 1:length(snp.fasta[[j]])){
    n[[j]][[i]] <- toupper(paste(n[[j]][[i]][snp.bed[[j]][i,3]:snp.bed[[j]][i,4]],collapse = ""))
  }
  names(n[[j]])<-annotation
}

for(j in 1:length(snp)){
  for(i in 1:length(snp.fasta[[j]])){
    snp[[j]][(2*i),] <- as.character(n[[j]][i])
  }
}
outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/test.fasta"#path## 
out_fileName <- sapply(names(snp),function(x){
  #paste(x, "sta", sep='')
  gsub("fa","fasta",x)
}) ##csv 
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath,x,sep='/')}) ## 
out_filePath
## 
for(i in 1:length(snp)){
  write.table(snp[[i]], file=out_filePath[i],col.names = F, quote = F,row.names = F)
}