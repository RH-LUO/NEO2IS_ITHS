#!/usr/local/bin/Rscript
myPaths <- .libPaths()  
#new <- c('/usr/local/lib64/R/library','/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')#c('/public/home/lorihan/miniconda3/lib/R/library')#
new <- '/public/home/lorihan/miniconda3/lib/R/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
library(seqinr)
library(stringr)
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/test.fasta" ## 

fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
faPath <- filePath[grep("indel.fa",filePath)]

path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/pos" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
bedPath <- filePath[grep("indel.bed",filePath)]
#####
indel <- lapply(faPath, function(x){
  read.table(x,sep = ",",fill = T,header = F)})  ## ， list
indel.fasta <- lapply(faPath, function(x){
  read.fasta(x)})  ## ， list
indel.bed <- lapply(bedPath, function(x){
  if (file.info(x)$size!=0) {
    read.table(x,sep = " ",fill = T,header = F)}})  

length(indel.fasta[[1]])
n1 <- indel.fasta[[1]]
#count(n1,1)
n <- indel.fasta
remove(annotation)
#annotation <- getAnnot(indel.fasta)
for(j in 1:length(indel)){
  annotation <- getAnnot(indel.fasta[[j]])
  indel.bed[[j]]$V6 <- nchar(indel.bed[[j]]$V5)
  indel.bed[[j]]$V4 <- ifelse(grepl("end",indel.bed[[j]]$V4), 
                              indel.bed[[j]]$V2+indel.bed[[j]]$V6-1,indel.bed[[j]]$V4)
  for(i in 1:length(indel.fasta[[j]])){
    n[[j]][[i]] <- toupper(paste(n[[j]][[i]][indel.bed[[j]][i,3]:indel.bed[[j]][i,4]],collapse = ""))
  }
  names(n[[j]])<-annotation
}

for(j in 1:length(indel)){
  for(i in 1:length(indel.fasta[[j]])){
    indel[[j]][(2*i),] <- as.character(n[[j]][i])
    indel[[j]][2*i-1,] <- sapply(str_split(indel[[j]][(2*i-1),],"[(]"),"[",1)
  }
}

outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/test.fasta"#path## 
out_fileName <- sapply(names(indel),function(x){
  #  gsub("fa","fsa",x)
  gsub("fa","fasta",x)
}) ##csv 
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ## 
out_filePath
## 
for(i in 1:length(indel)){
  write.table(indel[[i]], file=out_filePath[i],col.names = F, quote = F,row.names = F)
}
