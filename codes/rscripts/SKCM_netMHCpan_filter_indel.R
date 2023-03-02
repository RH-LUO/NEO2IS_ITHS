#!/usr/local/bin/Rscript
library(stringr)
###netMHCpan: >slurm-out
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/netMHCpan" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
sbPath <- filePath[grep("indel.sb",filePath)]
sbPath
indel.pep <- lapply(sbPath, function(x){
  read.csv(x,sep = "",header=T)})  ## ， list
indel.samples <- substr(names(indel.pep),1,16)
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/pos" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
bedPath <- filePath[grep("indel.bed",filePath)]
indel.bed <- lapply(bedPath, function(x){
  read.table(x,sep = " ",header = F)})
length(indel.pep)
length(indel.bed)
table(substr(names(indel.bed),1,16) %in% indel.samples)
indel.bed <- indel.bed[substr(names(indel.bed),1,16) %in% indel.samples]
indel.pep.new <- indel.pep
for(j in 1:length(indel.pep)){
  if (nrow(indel.pep[[j]])!=0) {
    col <- colnames(indel.pep[[j]])
    indel.bed[[j]]$V6 <- nchar(indel.bed[[j]]$V5)
    indel.bed[[j]]$V4 <- ifelse(grepl("end",indel.bed[[j]]$V4), 
                                indel.bed[[j]]$V2+indel.bed[[j]]$V6-1,indel.bed[[j]]$V4)
    identity=str_split(indel.pep[[j]]$Identity,"_p",simplify = T)[,1]
    identity <- ifelse(substr(identity,nchar(identity),nchar(identity))=="_",
                       str_split(identity,"_",simplify = T)[,1],identity)
    indel.pep[[j]]$site <- indel.bed[[j]]$V2[match(identity,indel.bed[[j]]$V1)]
    indel.pep[[j]]$site_1 <- indel.bed[[j]]$V3[match(identity,indel.bed[[j]]$V1)]
    indel.pep[[j]]$site_2 <- indel.bed[[j]]$V4[match(identity,indel.bed[[j]]$V1)]
    
    passed <- indel.pep[[j]]$site_2 > indel.pep[[j]]$site
    filter.1 <- indel.pep[[j]][passed,]$Pos+nchar(indel.pep[[j]][passed,]$Peptide)-1>=indel.pep[[j]][passed,]$site-indel.pep[[j]][passed,]$site_1+1# >=11
    table(filter.1)
    indel.pep.new[[j]] <- indel.pep[[j]][passed,][filter.1,]
    indel.pep[[j]] <- rbind(indel.pep[[j]][!passed,],indel.pep.new[[j]])
  }
}
#j=443 more TCGA-XC-AA0X-01A.indel.pep
outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.indel_results/test/netMHCpan"#path## 
out_fileName <- sapply(names(indel.pep),function(x){
  gsub("sb","pep",x)
}) ##csv 

out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ## 
out_filePath
## 
for(i in 1:length(indel.pep)){
  write.table(indel.pep[[i]][,col], file=out_filePath[i],sep = "\t",quote = F,row.names = F)
}
