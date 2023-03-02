#!/usr/local/bin/Rscript
library(stringr)
###netMHCpan: >nohup.out
path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/netMHCpan" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
sbPath <- filePath[grep("snp.sb",filePath)]
sbPath
snp.pep <- lapply(sbPath, function(x){
  read.csv(x,sep = "",header=T)})  ## ， list

path <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/pos" ## 
fileNames <- dir(path)  ## 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ## 
filePath
bedPath <- filePath[grep("snp.bed",filePath)]
snp.bed <- lapply(bedPath, function(x){
  read.table(x)})
for(j in 1:length(snp.pep)){
  if (nrow(snp.pep[[j]])!=0) {
    col <- colnames(snp.pep[[j]])
    identity <- str_split(snp.pep[[j]]$Identity,"_p",simplify = T)[,1]
    identity <- ifelse(substr(identity,nchar(identity),nchar(identity))=="_",
                       str_split(identity,"_",simplify = T)[,1],identity)
    #unique(sapply(str_split(snp.pep$Identity,"[_]"),"[",1))
    #match(sapply(str_split(snp.pep$Identity,"[_]"),"[",1),snp.bed$V1)
    #table(sapply(str_split(snp.pep[[j]]$Identity,"[_]"),"[",1) %in% snp.bed[[j]]$V1)
    #Error: sapply(str_split(snp.pep$Identity,"[_]"),"[",1) causing NA in snp.pep[[293]]
    # table(str_split(snp.pep[[j]]$Identity,"_p_",simplify = T)[,1] %in% snp.bed[[j]]$V1)
    table(is.na(str_split(snp.pep[[j]]$Identity,"_p_",simplify = T)[,1]))
    snp.pep[[j]][which(is.na(str_split(snp.pep[[j]]$Identity,"_p_",simplify = T)[,1])),]
    
    table(is.na(match(identity,snp.bed[[j]]$V1)))
    snp.bed[[j]][is.na(match(identity,snp.bed[[j]]$V1)),]
    table(identity %in% snp.bed[[j]]$V1)
    snp.pep[[j]][!identity %in% snp.bed[[j]]$V1,]
    identity[!identity %in% snp.bed[[j]]$V1]
    
    snp.pep[[j]]$site <- snp.bed[[j]]$V2[match(identity,snp.bed[[j]]$V1)]
    snp.pep[[j]]$site_1 <- snp.bed[[j]]$V3[match(identity,snp.bed[[j]]$V1)]
    snp.pep[[j]]$site_2 <- snp.bed[[j]]$V4[match(identity,snp.bed[[j]]$V1)]
    
    snp.pep[[j]][as.numeric(snp.pep[[j]]$Pos)+nchar(snp.pep[[j]]$Peptide)-1 < snp.pep[[j]]$site-snp.pep[[j]]$site_1+1,] 
    snp.pep[[j]][snp.pep[[j]]$site_1+snp.pep[[j]]$Pos-1 >snp.pep[[j]]$site,] 
    filter.1 <- snp.pep[[j]]$Pos+nchar(snp.pep[[j]]$Peptide)-1>=snp.pep[[j]]$site-snp.pep[[j]]$site_1+1
    table(filter.1)
    which(is.na(filter.1))
    #snp.pep[[j]][which(is.na(filter.1)),]
    snp.pep[[j]] <-snp.pep[[j]][filter.1,]
    filter.2 <- snp.pep[[j]]$Pos<=snp.pep[[j]]$site-snp.pep[[j]]$site_1+1
    table(filter.2)
    snp.pep[[j]] <-snp.pep[[j]][filter.2,]
  }
}
outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.snp_results/test/netMHCpan"#path## 
out_fileName <- sapply(names(snp.pep),function(x){
  #paste(x, "sta", sep='')
  gsub("sb","pep",x)
}) ##csv 

out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ## 
out_filePath
## 
for(i in 1:length(snp.pep)){
  write.table(snp.pep[[i]][,col], file=out_filePath[i],sep = "\t",quote = F,row.names = F)
}
