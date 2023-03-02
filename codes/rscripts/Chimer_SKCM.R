# title: "SKCM fusions"
# author: "Ruihan Luo"
# date: "Aug 10th,2022"
myPaths <- .libPaths()  
new <- '/public/apps/R-4.1.0/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
library(GenomicRanges)
library(stringr)
chimer_SKCM <- read.csv("/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM.txt",
                        fill = T,header = T,sep = "\t")
table(chimer_SKCM$Source)
length(unique(chimer_SKCM$Fusion_pair))#2248
length(unique(chimer_SKCM$Barcode.ID))#485
length(grep("-11A",chimer_SKCM$Barcode.ID))#0
#chimer_SKCM <- chimer_SKCM[-grep("-11A",chimer_SKCM$Barcode.ID),]
chimer_SKCM.star <- chimer_SKCM[chimer_SKCM$Source=="STARFusion",]
chimer_SKCM.FAWG <- chimer_SKCM[chimer_SKCM$Source=="TCGA FAWG",]
chimer_SKCM.TumorFusion <- chimer_SKCM[chimer_SKCM$Source=="PRADA(TumorFusion)",]
chimer_SKCM.Tophat <- chimer_SKCM[chimer_SKCM$Source=="Tophat-Fusion",]
chimer_SKCM.FusionScan <- chimer_SKCM[chimer_SKCM$Source=="FusionScan",]


chimer_SKCM.FAWG <- chimer_SKCM.FAWG[!chimer_SKCM.FAWG$Fusion_pair %in% c(
  chimer_SKCM.star$Fusion_pair,chimer_SKCM.TumorFusion$Fusion_pair,
  chimer_SKCM.Tophat$Fusion_pair,chimer_SKCM.FusionScan$Fusion_pair),]

chimer_SKCM.TumorFusion <- chimer_SKCM.TumorFusion[!chimer_SKCM.TumorFusion$Fusion_pair %in% c(
  chimer_SKCM.star$Fusion_pair,chimer_SKCM.FAWG$Fusion_pair,
  chimer_SKCM.Tophat$Fusion_pair,chimer_SKCM.FusionScan$Fusion_pair),]

chimer_SKCM.Tophat <- chimer_SKCM.Tophat[!chimer_SKCM.Tophat$Fusion_pair %in% c(
  chimer_SKCM.star$Fusion_pair,chimer_SKCM.TumorFusion$Fusion_pair,
  chimer_SKCM.FAWG$Fusion_pair,chimer_SKCM.FusionScan$Fusion_pair),]


chimer_SKCM.FusionScan <- chimer_SKCM.FusionScan[!chimer_SKCM.FusionScan$Fusion_pair %in% c(
  chimer_SKCM.star$Fusion_pair,chimer_SKCM.TumorFusion$Fusion_pair,
  chimer_SKCM.Tophat$Fusion_pair,chimer_SKCM.FAWG$Fusion_pair),]

length(unique(chimer_SKCM.star$Fusion_pair))#458
length(grep("-11A",chimer_SKCM.star$Barcode.ID))

chimer_SKCM.all <- rbind(chimer_SKCM.star,chimer_SKCM.TumorFusion,chimer_SKCM.FAWG
                         ,chimer_SKCM.Tophat,chimer_SKCM.FusionScan)#4036
length(unique(substr(chimer_SKCM.all$Barcode.ID,1,12)))#419
#chimer_SKCM <- rbind(chimer_SKCM.star,chimer_SKCM.TumorFusion#,chimer_SKCM.FAWG,
#                    ,chimer_SKCM.Tophat,chimer_SKCM.FusionScan)
length(unique(chimer_SKCM.all$Fusion_pair))#2248
dim(chimer_SKCM.all)#2650   26
table(chimer_SKCM.all$Source)
table(chimer_SKCM.all$Genome_build_ver)
chimer_SKCM <- chimer_SKCM.all[chimer_SKCM.all$Genome_build_ver=="hg38",]
chimer_SKCM.all <- chimer_SKCM.all[chimer_SKCM.all$Genome_build_ver!="hg38",]
length(unique(chimer_SKCM$Barcode.ID))
length(unique(substr(chimer_SKCM$Barcode.ID,1,12)))#192

length(unique(chimer_SKCM.all$Barcode.ID))
patients <- substr(chimer_SKCM.all$Barcode.ID,1,12)
length(unique(patients))#410

#chain <- import.chain("bosTau6.hg19.all.chain")
#BiocManager::install("rtracklayer")
hg19 <- str_split(chimer_SKCM.all$Gene5.Junction,"[(]",simplify = T)[,1]
hg19.left <- str_split(hg19, '[:]',simplify = T)[,2]
hg19.left <- as.numeric(str_split(hg19.left, '[ ]',simplify = T)[,2])
hg19.left.chr <- str_split(hg19, '[:]',simplify = T)[,1]
hg19.left.chr <- str_split(hg19, '[ ]',simplify = T)[,1]
hg19.left.strand<- str_split(chimer_SKCM.all$Gene5.Junction, '[)]',simplify = T)[,1]
hg19.left.strand <- str_split(hg19.left.strand, '[(]',simplify = T)[,2]
table(hg19.left.strand)
hg19.left.strand[hg19.left.strand=="'+'"] <- "+"
#hg19.left.strand[hg19.left.strand=="'-'"] <- "-"
hg19.left.strand <- as.character(hg19.left.strand)
Hg19.left <- GRanges(seqnames=Rle(hg19.left.chr),ranges=IRanges(hg19.left,(hg19.left+1)),strand = hg19.left.strand)
write.table(Hg19.left,"/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_left.bed",sep="\t",quote = F,row.names = F,col.names = F)
#####
hg19 <- str_split(chimer_SKCM.all$Gene3.Junction,"[(]",simplify = T)[,1]
hg19.right <- str_split(hg19, '[:]',simplify = T)[,2]
hg19.right <- as.numeric(str_split(hg19.right, '[ ]',simplify = T)[,2])
hg19.right.chr <- str_split(hg19, '[:]',simplify = T)[,1]
hg19.right.chr <- str_split(hg19, '[ ]',simplify = T)[,1]
hg19.right.strand<- str_split(chimer_SKCM.all$Gene3.Junction, '[)]',simplify = T)[,1]
hg19.right.strand <- str_split(hg19.right.strand, '[(]',simplify = T)[,2]
table(hg19.right.strand)
hg19.right.strand[hg19.right.strand=="'+'"] <- "+"
#hg19.right.strand[hg19.right.strand=="'-'"] <- "-"
hg19.right.strand <- as.character(hg19.right.strand)
Hg19.right <- GRanges(seqnames=Rle(hg19.right.chr),ranges=IRanges(hg19.right,(hg19.right+1)),strand = hg19.right.strand)
write.table(Hg19.right,"/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_right.bed",sep = "\t",quote = F,row.names = F,col.names = F)
setwd("/public/home/lorihan/lrh/SKCM/Fusion")
system("/public/home/lorihan/bin/liftOver chimer_SKCM_left.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz chimer_SKCM_left_hg38.bed chimer_SKCM_left_hg38_unMapped.bed")
system("/public/home/lorihan/bin/liftOver chimer_SKCM_right.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz chimer_SKCM_right_hg38.bed chimer_SKCM_right_hg38_unMapped.bed")
hg38.right_unmapped <- read.table("/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_right_hg38_unMapped.bed")
##1
unique(hg38.right_unmapped$V2)
hg38.right_unmapped <- hg38.right_unmapped[match(unique(hg38.right_unmapped$V2),hg38.right_unmapped$V2),2]
tmp.unmapped.right <- lapply(hg38.right_unmapped,function(x){
  tmp.unmapped <- grep(x,chimer_SKCM.all$Gene3.Junction)
})
unlist(tmp.unmapped.right)
hg38.left_unmapped <- read.table("/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_left_hg38_unMapped.bed")
unique(hg38.left_unmapped$V2)
hg38.left_unmapped <- hg38.left_unmapped[match(unique(hg38.left_unmapped$V2),hg38.left_unmapped$V2),2]
tmp.unmapped.left <- lapply(hg38.left_unmapped,function(x){
  tmp.unmapped <- grep(x,chimer_SKCM.all$Gene5.Junction)
})
unlist(tmp.unmapped.left)
hg38.left_unmapped[3]
grep(hg38.left_unmapped[3],chimer_SKCM.all$Gene5.Junction)
hg19_unmapped <- unique(c(unlist(tmp.unmapped.left),unlist(tmp.unmapped.right)))
chimer_SKCM.hg19 <- chimer_SKCM.all[-hg19_unmapped,]
hg19 <- str_split(chimer_SKCM.hg19$Gene5.Junction,"[(]",simplify = T)[,1]
hg19.left <- str_split(hg19, '[:]',simplify = T)[,2]
hg19.left <- as.numeric(str_split(hg19.left, '[ ]',simplify = T)[,2])
hg19.left.chr <- str_split(hg19, '[:]',simplify = T)[,1]
hg19.left.chr <- str_split(hg19, '[ ]',simplify = T)[,1]
hg19.left.strand<- str_split(chimer_SKCM.hg19$Gene5.Junction, '[)]',simplify = T)[,1]
hg19.left.strand <- str_split(hg19.left.strand, '[(]',simplify = T)[,2]
table(hg19.left.strand)
hg19.left.strand[hg19.left.strand=="'+'"] <- "+"
#hg19.left.strand[hg19.left.strand=="'-'"] <- "-"
hg19.left.strand <- as.character(hg19.left.strand)
Hg19.left <- GRanges(seqnames=Rle(hg19.left.chr),ranges=IRanges(hg19.left,(hg19.left+1)),strand = hg19.left.strand)
write.table(Hg19.left,"/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_hg19_left.bed",sep="\t",quote = F,row.names = F,col.names = F)
#####
hg19 <- str_split(chimer_SKCM.hg19$Gene3.Junction,"[(]",simplify = T)[,1]
hg19.right <- str_split(hg19, '[:]',simplify = T)[,2]
hg19.right <- as.numeric(str_split(hg19.right, '[ ]',simplify = T)[,2])
hg19.right.chr <- str_split(hg19, '[:]',simplify = T)[,1]
hg19.right.chr <- str_split(hg19, '[ ]',simplify = T)[,1]
hg19.right.strand<- str_split(chimer_SKCM.hg19$Gene3.Junction, '[)]',simplify = T)[,1]
hg19.right.strand <- str_split(hg19.right.strand, '[(]',simplify = T)[,2]
table(hg19.right.strand)
hg19.right.strand[hg19.right.strand=="'+'"] <- "+"
#hg19.right.strand[hg19.right.strand=="'-'"] <- "-"
hg19.right.strand <- as.character(hg19.right.strand)
Hg19.right <- GRanges(seqnames=Rle(hg19.right.chr),ranges=IRanges(hg19.right,(hg19.right+1)),strand = hg19.right.strand)
write.table(Hg19.right,"/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_hg19_right.bed",sep = "\t",quote = F,row.names = F,col.names = F)

system("/public/home/lorihan/bin/liftOver chimer_SKCM_hg19_left.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz chimer_SKCM_left_hg38_new.bed chimer_SKCM_left_hg38_unMapped_new.bed")
system("/public/home/lorihan/bin/liftOver chimer_SKCM_hg19_right.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz chimer_SKCM_right_hg38_new.bed chimer_SKCM_right_hg38_unMapped_new.bed")
#####
hg38.left <- read.table("/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_left_hg38_new.bed")
#hg38.left <- GRanges(seqnames=Rle(hg38.left[,1]),ranges=IRanges(hg38.left[,2],hg38.left[,3]),strand = hg38.left[,5])
table(hg19.left.strand==hg38.left$V5)
hg38.right <- read.table("/public/home/lorihan/lrh/SKCM/Fusion/chimer_SKCM_right_hg38_new.bed")
#hg38.right <- GRanges(seqnames=Rle(hg38.right[,1]),ranges=IRanges(hg38.right[,2],hg38.right[,3]),strand = hg38.right[,5])
table(hg19.right.strand==hg38.right$V5)

chimer_SKCM.hg38 <- chimer_SKCM
chimer_SKCM <- rbind(chimer_SKCM.hg19,chimer_SKCM.hg38)
Left.g <- str_split(chimer_SKCM$Fusion_pair,"[-]",simplify = T)[,1]
LeftGene <- chimer_SKCM$H_gene
table(is.na(match(Left.g,LeftGene)))
LeftGene[is.na(match(Left.g,LeftGene))]
Left.g[is.na(match(LeftGene,Left.g))]

length(unique(LeftGene))#1814
l <- 1:length(chimer_SKCM$Fusion_pair)
RightGene <- strsplit(chimer_SKCM$Fusion_pair,"-")
grep("CTD-2589C9.4",chimer_SKCM$Fusion_pair)#i=2029

for (i in 1:length(l)) {
  l[i] <- length(strsplit(chimer_SKCM$Fusion_pair[i],"-")[[1]])
}
table(l)
for (i in 1:length(l)) {
  l[i] <- length(strsplit(chimer_SKCM$Fusion_pair[i],"-")[[1]])
  logi.left.g <- grepl(LeftGene[i],RightGene[[i]][1])
  RightGene[[i]] <- ifelse(l[i]>2,ifelse(logi.left.g,paste(RightGene[[i]][l[i]-1],RightGene[[i]][l[i]],sep = "-"),
                                         RightGene[[i]][l[i]]),RightGene[[i]][l[i]])
}

RightGene <- unlist(RightGene)
length(unique(RightGene))
#####
LeftBreakpoint=paste(hg38.left[,1],hg38.left[,2],sep = ":")
left.breakpoint.hg38 <- chimer_SKCM.hg38$Gene5.Junction
hg38.left.strand<- str_split(left.breakpoint.hg38, '[)]',simplify = T)[,1]
hg38.left.strand <- str_split(hg38.left.strand, '[(]',simplify = T)[,2]

left.breakpoint.hg38 <- gsub(" : ",":",left.breakpoint.hg38)
left.breakpoint.hg38 <- gsub("[(+)]","",left.breakpoint.hg38)
left.breakpoint.hg38 <- str_replace(left.breakpoint.hg38,"[-]","")#gsub("[(-)]","",left.breakpoint.hg38)
LeftBreakpoint <- c(LeftBreakpoint,left.breakpoint.hg38)

RightBreakpoint=paste(hg38.right[,1],hg38.right[,2],sep = ":")
right.breakpoint.hg38 <- chimer_SKCM.hg38$Gene3.Junction
hg38.right.strand<- str_split(right.breakpoint.hg38, '[)]',simplify = T)[,1]
hg38.right.strand <- str_split(hg38.right.strand, '[(]',simplify = T)[,2]

right.breakpoint.hg38 <- gsub(" : ",":",right.breakpoint.hg38)
right.breakpoint.hg38 <- gsub("[(+)]","",right.breakpoint.hg38)
right.breakpoint.hg38 <- str_replace(right.breakpoint.hg38,"[-]","")#gsub("[(-)]","",right.breakpoint.hg38)
RightBreakpoint <- c(RightBreakpoint,right.breakpoint.hg38)

patients <- substr(chimer_SKCM$Barcode.ID,1,12)
length(unique(patients))#485

fusions_SKCM <- data.frame(LeftGene=LeftGene,LeftBreakpoint=LeftBreakpoint,
                      LeftStrand=c(hg38.left[,5],hg38.left.strand),
                      RightGene=RightGene,RightBreakpoint=RightBreakpoint,
                      RightStrand=c(hg38.right[,5],hg38.right.strand),
                      Barcode = paste(patients,"-01A",sep = ""))#460unique patients chimer_SKCM$Barcode.ID)
chimer_skcm <- cbind(fusions_SKCM,chimer_SKCM)
length(unique(LeftGene))#1646
length(unique(RightGene))#1813
#######transfer######
load("/public/home/lorihan/lrh/NSCLC/gene_info.RData")
gene_SKCM.l <- bitr(unique(LeftGene), toType = "REFSEQ",
                    fromType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db",drop=TRUE)
#test.g.tr <- ifelse(!is.na(match(test.g,gene.df$REFSEQ)),
#                   gene.df$SYMBOL[match(test.g,gene.df$REFSEQ)],test.g)
#test.g.tr <-  gene_SKCM.l$SYMBOL[match(test.g,gene_SKCM.l$REFSEQ)]
table(unique(LeftGene) %in% gene_SKCM.l$SYMBOL)#1510
table(unique(LeftGene) %in% gene_info$symbols)#1516

gene_SKCM.r <- bitr(unique(RightGene), toType = "REFSEQ",
                    fromType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db",drop=TRUE)
#test.g.tr <- ifelse(!is.na(match(test.g,gene.df$REFSEQ)),
#                   gene.df$SYMBOL[match(test.g,gene.df$REFSEQ)],test.g)
#test.g.tr <-  gene_SKCM.l$SYMBOL[match(test.g,gene_SKCM.l$REFSEQ)]
table(unique(RightGene) %in% gene_SKCM.r$SYMBOL)#1630
table(unique(RightGene) %in% gene_info$symbols)#1666


####start to find CDS regions####
#######
b=read.table('/public/home/lorihan/Pyscript/g2t2cds.txt')
colnames(b)=c('gene','transcript','exon','chr','start','end','exon_length')
head(b)
########
library(org.Hs.eg.db)
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2ensembl=toTable(org.Hs.egENSEMBL)
length(unique(eg2ensembl$ensembl_id))
eg2name=toTable(org.Hs.egGENENAME)
eg2alias=toTable(org.Hs.egALIAS2EG)
###
eg2ncbi=toTable(org.Hs.egREFSEQ)
eg2alis_list=lapply(split(eg2alias,eg2alias$gene_id),function(x){paste0(x[,2],collapse = ";")})
GeneList=mappedLkeys(org.Hs.egSYMBOL)
if( GeneList[1] %in% eg2symbol$symbol ){
  symbols=GeneList
  geneIds=eg2symbol[match(symbols,eg2symbol$symbol),'gene_id']
}else{
  geneIds=GeneList
  symbols=eg2symbol[match(geneIds,eg2symbol$gene_id),'symbol']
}
geneNames=eg2name[match(geneIds,eg2name$gene_id),'gene_name']
geneEnsembl=eg2ensembl[match(geneIds,eg2ensembl$gene_id),'ensembl_id']
head(eg2ncbi)
geneRefSeq=eg2ncbi[match(geneIds,eg2ncbi$gene_id),'accession']
geneAlias=sapply(geneIds,function(x){ifelse(is.null(eg2alis_list[[x]]),"no_alias",eg2alis_list[[x]])})
createLink <- function(base,val) {
  sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>',base,val) 
  ##  target="_blank" 
}
gene_info=data.frame(   symbols=symbols,
                        #                      geneIds=createLink(paste0("http://www.ncbi.nlm.nih.gov/gene/",geneIds),geneIds),
                        #                     geneNames=geneNames,
                        geneEnsembls=geneEnsembl,
                        geneAlias=geneAlias,
                        geneRefSeqs=geneRefSeq,
                        stringsAsFactors = F
) 
id_mapped <- read.table("/public/home/lorihan/lrh/SKCM/gencode.v22.annotation.gene.probeMap",header = T)
table(unique(LeftGene) %in% id_mapped$gene)
table(unique(LeftGene) %in% gene_info$symbols)
#str_split(mapped_id$V1,"[.]",simplify = T)[,1]

#####LeftGene#####
#load("/public/home/lorihan/lrh/NSCLC/gene_info.RData")
length(unique(LeftGene))#1814 --> 2077
gene.df.left <- gene_info[gene_info$symbols %in% unique(LeftGene),]
#gene.df.left <- id_trans[id_trans$id %in% unique(LeftGene),]
#length(unique(LeftGene)[na.omit(match(gene.df.left$symbols,unique(LeftGene)))])
head(gene.df.left[1:2],10)
head(unique(LeftGene)[na.omit(match(gene.df.left$symbols,unique(LeftGene)))],10)
unique(LeftGene)[!unique(LeftGene) %in% gene_info$symbols]
gene.df.left$LeftGene <- unique(LeftGene)[na.omit(match(gene.df.left$symbols,unique(LeftGene)))]
#gene.df.left$LeftGene <- unique(LeftGene)[na.omit(match(gene.df.left$id,unique(LeftGene)))]

#match(gene.df.left$symbols,gene.df.left$LeftGene)
#grepl(unique(LeftGene)[1],gene_info$geneAlias)
######
nonspecific<-unique(LeftGene)[!unique(LeftGene) %in% gene_info$symbols]
##nonspecific<-unique(LeftGene)[!unique(LeftGene) %in% id_trans$id]
nonspecific
test.logi.2 <- rep("FALSE",length(nonspecific))
alias <- unlist(str_split(as.vector(gene_info$geneAlias),";"))
for (j in 1:length(nonspecific)) {
  #test.logi.2[j] <- length(unique(grepl(nonspecific[j],as.vector(gene_info$geneAlias))))==2
  test.logi.2[j] <- !is.na(match(nonspecific[j],alias))
}
nonspecific <- nonspecific[test.logi.2=="TRUE"]
test.n<-vector("list",length(nonspecific))
names(test.n) <- nonspecific
count <- 1:length(nonspecific)
for (j in 1:length(nonspecific)) {
  # test.n[[j]] <- match(nonspecific[j], as.vector(unlist(alias[[j]])))
  test.n[[j]] <- grep(nonspecific[j], as.vector(gene_info$geneAlias))
  count[j] <- length(test.n[[j]])
}
unlist(test.n)
sum(count)
nonspecific_list<-gene_info[unlist(test.n),]
nonspecific_list$LeftGene <- rep(names(test.n),count)

alias <- str_split(as.vector(nonspecific_list$geneAlias),";")
names(alias) <- nonspecific_list$LeftGene
test.logi.1 <- rep("FALSE",length(alias))
for (j in 1:length(alias)) {
  test.logi.1[[j]] <- !is.na(match(nonspecific_list$LeftGene[j], as.vector(unlist(alias[[j]]))))
  # count[j] <- length(test.n[[j]])
}
table(test.logi.1)
nonspecific_list <- nonspecific_list[test.logi.1=="TRUE",]
dim(nonspecific_list)#125
dim(gene.df.left)
gene.df.l <- rbind(gene.df.left,nonspecific_list)
#colnames(nonspecific_list) <- colnames(gene.df.left)
#gene.df.l <- rbind(gene.df.left[,-3],nonspecific_list[,-3])
#####
dim(gene.df.l)
table(gene.df.l$symbols!=gene.df.l$LeftGene)#121
#table(gene.df.l$id!=gene.df.l$LeftGene)#65
table(is.na(gene.df.l$geneEnsembls))#15
table(unique(LeftGene) %in% gene.df.l$LeftGene)#1630
table(unique(LeftGene) %in% gene_SKCM.l$SYMBOL)#1510
gene.df.left <- gene.df.l#[gene.df.l$geneEnsembls %in% b$gene,]

#####RightGene#####
gene.df.right <- gene_info[gene_info$symbols %in% unique(RightGene),]
unique(RightGene)[!unique(RightGene) %in% gene_info$symbols]#214
gene.df.right$RightGene <- unique(RightGene)[na.omit(match(gene.df.right$symbols,unique(RightGene)))]

nonspecific<-unique(RightGene)[!unique(RightGene) %in% gene_info$symbols]
nonspecific
test.logi.2 <- rep("FALSE",length(nonspecific))
alias <- unlist(str_split(as.vector(gene_info$geneAlias),";"))
for (j in 1:length(nonspecific)) {
  #test.logi.2[j] <- length(unique(grepl(nonspecific[j],as.vector(gene_info$geneAlias))))==2
  test.logi.2[j] <- !is.na(match(nonspecific[j],alias))
}
nonspecific <- nonspecific[test.logi.2=="TRUE"]
test.n<-vector("list",length(nonspecific))
names(test.n) <- nonspecific
count <- 1:length(nonspecific)
for (j in 1:length(nonspecific)) {
  # test.n[[j]] <- match(nonspecific[j], as.vector(unlist(alias[[j]])))
  test.n[[j]] <- grep(nonspecific[j], as.vector(gene_info$geneAlias))
  count[j] <- length(test.n[[j]])
}
unlist(test.n)
sum(count)
nonspecific_list<-gene_info[unlist(test.n),]
nonspecific_list$RightGene <- rep(names(test.n),count)

alias <- str_split(as.vector(nonspecific_list$geneAlias),";")
names(alias) <- nonspecific_list$RightGene
test.logi.1 <- rep("FALSE",length(alias))
for (j in 1:length(alias)) {
  test.logi.1[[j]] <- !is.na(match(nonspecific_list$RightGene[j], as.vector(unlist(alias[[j]]))))
  # count[j] <- length(test.n[[j]])
}
table(test.logi.1)
nonspecific_list <- nonspecific_list[test.logi.1=="TRUE",]
dim(nonspecific_list)#123
gene.df.r <- rbind(gene.df.right,nonspecific_list)
table(gene.df.r$symbols!=gene.df.r$RightGene)#123
table(is.na(gene.df.r$geneEnsembls))#32
#####
table(unique(RightGene) %in% gene.df.r$symbols)#1666
table(unique(RightGene) %in% gene.df.r$RightGene)#1786
table(unique(RightGene) %in% gene_SKCM.r$SYMBOL)#1630
gene.df.right <- gene.df.r#[gene.df.r$geneEnsembls %in% b$gene,]
table(is.na(gene.df.right$geneEnsembls))#1629
gene.df.right[1:5,]

#####merge#####
gene.df <- data.frame(LeftGene=LeftGene, RightGene=RightGene)
gene.df$LeftEnsembl <- gene.df.left$geneEnsembls[match(gene.df$LeftGene,gene.df.left$LeftGene)]
gene.df$RightEnsembl <- gene.df.right$geneEnsembls[match(gene.df$RightGene,gene.df.right$RightGene)]
gene.df[1:5,]
length(unique(gene.df$LeftGene));length(unique(gene.df$RightGene))
#rownames(gene.df) <- rownames(left.g)
table(unique(LeftGene) %in% gene_SKCM.l$SYMBOL)#1510
table(unique(LeftGene) %in% gene.df.l$LeftGene)#1630
gene_SKCM.l <- gene.df.l
table(gene_SKCM.l$symbols!=gene_SKCM.l$LeftGene)#134
gene_SKCM.l.missed <- gene_SKCM.l[gene_SKCM.l$symbols!=gene_SKCM.l$LeftGene,]
####
gene_SKCM.r <- gene.df.r
table(gene_SKCM.r$symbols!=gene_SKCM.r$RightGene)#187
gene_SKCM.r.missed <- gene_SKCM.r[gene_SKCM.r$symbols!=gene_SKCM.r$RightGene,]
####
length(unique(gene_SKCM.r$RightGene))#1786
length(unique(RightGene))# 1813
gene.fusion_SKCM <- rbind(gene_SKCM.r[,c("symbols","geneRefSeqs")],gene_SKCM.l[,c("symbols","geneRefSeqs")])
length(unique(gene.fusion_SKCM$symbols))#2939
length(unique(gene.fusion_SKCM$geneRefSeqs))#2898
######
#######
colnames(gene_SKCM.r)[5] <- "fusions"
colnames(gene_SKCM.l)[5] <- "fusions"
gene.fusion <- rbind(gene_SKCM.l,gene_SKCM.r)
length(unique(gene.fusion$fusions))#3025
length(unique(gene.fusion$symbols))#2939
length(unique(gene.fusion$geneRefSeqs))#2898
gene.fusion$geneRefSeqs[is.na(gene.fusion$geneRefSeqs)]=""
table(gene.fusion$geneRefSeqs=="")
gene.fusion <- gene.fusion[match(unique(gene.fusion$symbols),gene.fusion$symbols),]
#####
#refFlat <- read.table("/public/home/lorihan/lrh/LUAD/PairedData/Genefuse/refFlat.txt.gz",header = F)
refFlat <- read.table("/public/home/lorihan/lrh/LUAD/PairedData/Genefuse/refFlat.txt",header = F)
table(gene.fusion$symbols %in% refFlat$V1)#2889
table(gene.fusion$fusions %in% refFlat$V1)#2770
fusion.missed <- gene.fusion[!gene.fusion$fusions %in% refFlat$V1,]
fusion.missed <- gene.fusion[!gene.fusion$symbols %in% refFlat$V1,]
table(is.na(match(gene.fusion$symbols,refFlat$V1)))
gene.fusion <- gene.fusion[gene.fusion$symbols %in% refFlat$V1,]
gene.fusion$gene2ncbi <- refFlat$V2[match(gene.fusion$symbols,refFlat$V1)]
gene.fusion.unconsist <- gene.fusion[gene.fusion$geneRefSeqs!=gene.fusion$gene2ncbi,]
gene.fusion$gene2ncbi[is.na(gene.fusion$gene2ncbi)]=""
table(gene.fusion$gene2ncbi=="")
system("mkdir /public/home/lorihan/lrh/SKCM/Genefuse") 
write.table(gene.fusion[,c(1,6)],quote = F,row.names = F,col.names = F,sep = "\t",
            "/public/home/lorihan/lrh/SKCM/Genefuse/TCGA-SKCM.fusion.txt")

#####
r.fusions <- split(fusions_SKCM, fusions_SKCM$Barcode)
for (i in 1:length(names(r.fusions))) {
  for (j in 1:nrow(r.fusions[[i]])) {
    r.fusions[[i]]$Left[j] <- gene.df.left$geneEnsembls[match(r.fusions[[i]]$LeftGene[j],gene.df.left$LeftGene)]
    r.fusions[[i]]$Right[j] <- gene.df.right$geneEnsembls[match(r.fusions[[i]]$RightGene[j],gene.df.right$RightGene)]
    
  }
}
#######
#rm(list = ls())
options(stringsAsFactors = F)
f <- readLines('/public/home/lorihan/lrh/hg38/Homo_sapiens.GRCh38.cds.all.fa.gz')#release-100(2021-11-05)
p <- readLines('/public/home/lorihan/Pyscript/out.pep.fa')
head(p)
length(p)
length(f)
pep=paste0(p,collapse = '')
pep=strsplit(pep,'>')[[1]]
pep=pep[-1]
fa<-f
head(fa)
head(fa[grep(">ENST",fa)],3)
fa.changed<- str_split(fa[grep(">ENST",fa)],'[ cds]',simplify = T)[,1]
head(fa.changed)
fa[grep(">ENST",fa)] <- fa.changed
fa=paste0( fa ,collapse = '')
fa=strsplit(fa,'>')[[1]]
fa=fa[-1]
head(fa)
length(pep)
length(fa)

#tid=str_split(fa,'[.]',simplify = T)[,1]
tid=sapply(str_split(fa,'[.]'),"[",1)
length(tid)

#seq=str_split(fa,']',simplify = T)[,2]
seq=sapply(str_split(fa,'[.]'),"[",2)
#seq=str_split(fa,'.',simplify = T)[,2]
seq <- substr(seq,2,nchar(seq))
head(tid)
head(seq)
length(seq)
table(seq=="")
table(is.na(seq))
#tid.p=str_split(pep,'[.]',simplify = T)[,1]
tid.p=sapply(str_split(pep,'[.]'),"[",1)
length(tid.p)
#seq.p=str_split(pep,']',simplify = T)[,3]
seq.p=sapply(str_split(pep,'standard]'),"[",2)
table(seq.p=="")
table(is.na(seq.p))

head(tid.p)
head(seq.p)

match(tid.p,tid)
length(seq)
length(seq.p)
#######
library(stringr)
re.fusions <- r.fusions
for (k in 1:length(r.fusions)) {
  c=r.fusions[[k]][,-7]
  re.l=NULL
  re.r=NULL
  re=NULL
  tmp = apply(c,1,function(x){
    ###LeftBreak###
    #x=c[1,]
    g = x[7]
    pos=as.numeric(str_split(x[2],':')[[1]][2])
    strand=x[3]
    info=b[b[,1]==g,] ## one gene might has more than 1 transcripts.
    lapply(split(info,info[,2]), function(y){
      ## each transcript 
      #y<-split(info,info[,2])[[2]]
      apply(y,1,function(z){
        #      z<-y[1,]
        ## each exon.  
        if(strand=="+"){
          if(z[5] <= pos & z[6] >= pos ){
            #print(c(z,pos))
            new = sum(as.numeric(y[y[,3] < as.numeric(z[3]),7])) + pos - as.numeric(z[5])+1 
            print(c(z,pos,new))
            re.l <<- rbind(re.l,c(z,pos,new,strand))
          }
        }else{#strand=="-"
          if(z[5] <= pos & z[6] >= pos){
            new = sum(as.numeric(y[y[,3] > as.numeric(z[3]),7])) + pos - as.numeric(z[5])+1 
            print(c(z,pos,new))
            re.l <<- rbind(re.l,c(z,pos,new,strand))}
        }
      })## each exon  
    }
    )## each transcript 
    ###RightBreak#
    g = x[8]
    pos=as.numeric(str_split(x[5],':')[[1]][2])
    strand=x[6]
    info=b[b[,1]==g,] ## one gene might has more than 1 transcripts.
    lapply(split(info,info[,2]), function(y){
      ## each transcript 
      #y<-split(info,info[,2])[[2]]
      apply(y,1,function(z){
        # z<-y[1,]
        ## each exon.  
        if(strand=="+"){
          if(z[5] <= pos & z[6] >= pos ){
            #print(c(z,pos))
            new = sum(as.numeric(y[y[,3] < as.numeric(z[3]),7])) + pos - as.numeric(z[5])+1 
            print(c(z,pos,new))
            re.r <<- rbind(re.r,c(z,pos,new,strand))
          }
        }
        else{#strand=="-"
          if(z[5] <= pos & z[6] >= pos){
            new = sum(as.numeric(y[y[,3] > as.numeric(z[3]),7])) + pos - as.numeric(z[5])+1 
            print(c(z,pos,new))
            re.r <<- rbind(re.r,c(z,pos,new,strand))}
        }
      })## each exon  
    })  ## each transcript
  }) ## each position 
  if(nrow(as.data.frame(re.r))!=0){
    colnames(re.r)=c('rightgene','right_transcript','exon','chr','start','end','exon_length','pos','right_new_pos','right_strand')
  }#else{  re.r=NULL}
  if(nrow(as.data.frame(re.l))!=0){
    colnames(re.l)=c('leftgene','left_transcript','exon','chr','start','end','exon_length','pos','left_new_pos','left_strand')
  }#else{ re.l=NULL}
  
  re.right <- as.data.frame(re.r)
  re.left <- as.data.frame(re.l) 
  #cbind(re.left,re.right)
  #  if (nrow(re.left)!=0) {
  if (length(grep("left",colnames(re.left)))!=0) {
    unique(re.left$left_transcript)
    match(re.left$left_transcript,tid)
    table(is.na(match(re.left$left_transcript,tid)))
    table(re.left$left_transcript %in% tid)
    re.left=re.left[re.left$left_transcript %in% tid,]
    re.left$left_tseq=seq[match(re.left$left_transcript,tid)]
    re.left$left_pseq=seq.p[match(re.left$left_transcript,tid.p)]
    ######
    re.left$up=unlist(lapply(1:nrow(re.left),function(i){
      #i=1
      new_pos.1=as.numeric(re.left[i,"left_new_pos"])
      l=nchar(re.left[i,"left_tseq"])
      strand=re.left[i,"left_strand"]
      if(strand=="-"){
        return(substring(re.left[i,"left_pseq"],(l-new_pos.1)/3,(l-new_pos.1)/3+12))
      }
      else{#strand=="+"
        return(substring(re.left[i,"left_pseq"],new_pos.1/3-12,new_pos.1/3))
      }
    }))
    re.left$l_symbol <- c$LeftGene[match(re.left$leftgene,c$Left)]#gene.df$LeftGene[match(re.left$leftgene,gene.df$LeftEnsembl)]
    re.left$r_symbol <- c$RightGene[match(re.left$leftgene,c$Left)]
  } else{
    re.left <- NULL
  }
  #####
  if (length(grep("right",colnames(re.right)))!=0) {
    unique(re.right$right_transcript)
    match(re.right$right_transcript,tid)
    table(is.na(match(re.right$right_transcript,tid)))
    table(re.right$right_transcript %in% tid)
    re.right=re.right[re.right$right_transcript %in% tid,]
    re.right$right_tseq=seq[match(re.right$right_transcript,tid)]
    re.right$right_pseq=seq.p[match(re.right$right_transcript,tid.p)]
    re.right[1,]
    #####  
    re.right$down=unlist(lapply(1:nrow(re.right),function(i){
      #i=1
      new_pos.1=as.numeric(re.right[i,"right_new_pos"])
      l=nchar(re.right[i,"right_tseq"])
      strand=re.right[i,"right_strand"]
      if(strand=="-"){
        return(substring(re.right[i,"right_pseq"],(l-new_pos.1)/3-12,(l-new_pos.1)/3))
      }
      else{#strand=="+"
        return(substring(re.right[i,"right_pseq"],new_pos.1/3+1,new_pos.1/3+13))
      }
    }))
    re.right$r_symbol <- c$RightGene[match(re.right$rightgene,c$Right)]#gene.df$RightGene[match(re.right$rightgene,gene.df$RightEnsembl)]
    re.right$l_symbol <- c$LeftGene[match(re.right$rightgene,c$Right)]
  }else{
    re.right <- NULL
  }
  ####
  if (length(grep("left",colnames(re.left)))!=0) {
    #re.right$down[re.right$l_symbol %in% re.left$l_symbol]
    for (i in 1:nrow(re.left)) {
      # re.left$down[re.left$l_symbol %in% re.right$l_symbol] <- paste(unique(re.right$down[re.right$l_symbol %in% re.left$l_symbol]),collapse = ",")
      re.left$down[i] <- paste(unique(re.right$down[grep(re.left$l_symbol[i],re.right$l_symbol)]),collapse = ",") 
    }
  }
  #re.left$down <- re.right$down[match(re.left$l_symbol,re.right$l_symbol)]
  #re.right$up <- re$up[match(re.right$r_symbol,re.1$r_symbol)]
  if (length(grep("right",colnames(re.right)))!=0) {
    #re.left$up[re.left$r_symbol %in% re.right$r_symbol]
    #re.right$up[re.right$r_symbol %in% re.left$r_symbol] <- paste(unique(re.left$up[re.left$r_symbol %in% re.right$r_symbol]),collapse = ",")
    for (j in 1:nrow(re.right)) {
      re.right$up[j] <- paste(unique(re.left$up[grep(re.right$r_symbol[j],re.left$r_symbol)]),collapse = ",") 
    }
  }
  if(nrow(as.data.frame(re.left))!=0 | nrow(as.data.frame(re.right))!=0){
    re.f <- rbind(re.left[,c("l_symbol","r_symbol","up","down")],
                  re.right[,c("l_symbol","r_symbol","up","down")]) 
    re.f$up[is.na(re.f$up)]=""
    re.f$down[is.na(re.f$down)]=""
    
    re.f$candidate <- apply(re.f,1,function(x){
      candidate <- ifelse(grepl(",",x["up"]),
                          paste(paste(unlist(str_split(x["up"],",")),x["down"],sep = ""),collapse = ";"),   
                          paste(paste(x["up"],unlist(str_split(x["down"],",")),sep = ""),collapse = ";")    
      )
    })
    re.f$fusion_pair <- paste(re.f$l_symbol,re.f$r_symbol,sep = "-")
    unique(re.f$dusion_pair)
    unique(re.f$candidate)
    
    re <- re.f[match(unique(re.f$candidate),re.f$candidate),c("fusion_pair","up","down","candidate")]
    class(c)
    print(colnames(re))
    data.fusion <- NULL
    tmp <- apply(re,1,function(i){
      data.fusion <<- #ifelse(grepl(";",i["candidate"]),
        c(data.fusion,paste(i["fusion_pair"],
                            unlist(str_split(i["candidate"],";")),sep = ";"))#,
    })
    length(unique(data.fusion))    
    data.fusion <- unique(data.fusion)
    data.fusion <- paste(">",data.fusion,sep = "")
    data.fa<-NULL
    tmp <- sapply(data.fusion, function(x){
      pep<- str_split(x,";")[[1]][2]
      name <- str_split(x,";")[[1]][1]
      data.fa <<- c(data.fa,c(name,pep))
    })
    re.fusions[[k]] <- data.fa
  }else{
    re.fusions[[k]]="NULL"
  }
}
re.fusions[[174]]
save(chimer_SKCM,fusions_SKCM,chimer_skcm,r.fusions, re.fusions,
     file="/public/home/lorihan/lrh/SKCM/Neoantigen/new.fusion.Rdata")#/Input.fusion_tsvs/re.fusions.Rdata")
system("mkdir /public/home/lorihan/lrh/SKCM/Neoantigen/Input.fusion_tsvs")
outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/Input.fusion_tsvs"#path
out_fileName <- sapply(names(re.fusions),function(x){
  paste(x, "fasta", sep='.')
}) ##csv
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath,x,sep='/')}) ##
out_filePath
##
for(i in 1:length(re.fusions)){
  write.table(re.fusions[[i]], file=out_filePath[i], quote = F,row.names = F,col.names = F)
}
setwd("/public/home/lorihan/lrh/SKCM")
load("SKCM_202203_result.Rdata")
HLA.type <- read.table("/public/home/lorihan/lrh/SKCM/TCGA-SKCM-hlaTypesAll.tsv",header = T,fill = T)
dim(HLA.type)#469   8  
match(substr(names(re.fusions),1,12),HLA.type$patientBarcode)
table(substr(names(re.fusions),1,12)%in% HLA.type$patientBarcode)
table(substr(names(re.fusions),1,12)%in%  TCGA_SKCM_VCF$Patients)
HLA.type <- HLA.type[match(substr(names(re.fusions),1,12),
                           HLA.type$patientBarcode),]
HLA.type[1:5,]
names(re.fusions)[1:5]
#HLA.type <- HLA.type[HLA.type$patientBarcode %in% TCGA_SKCM_VCF$Patients,]
TCGA_SKCM_VCF[1:5,]
HLA.type[1:5,]
HLA.type <- HLA.type[,-2]
HLA.type[,-1]<-apply(HLA.type[,-1], 2, function(values){
  values <- as.character(values)
  values <- gsub("[*]","",values)
})

for (i in 1:(dim(HLA.type)[1]-1)){
  for (j in 1:(dim(HLA.type)[2]-1)){
    HLA.type[i,j] <- ifelse(HLA.type[i,j]==HLA.type[i,j+1],NA,HLA.type[i,j])
  }
}

HLA.type[1:5,]
match(names(re.fusions),paste(HLA.type$patientBarcode,"01A",sep = "-"))
HLA.type$patientBarcode <- names(re.fusions)#paste(HLA.type$patientBarcode,"01A",sep = "-")
HLA.type[1:5,1:5]
dim(HLA.type)#319

setwd("/public/home/lorihan/lrh/SKCM/")
output <- 1:nrow(HLA.type)
for(i in 1:nrow(HLA.type)){
  output[i] <- paste(HLA.type[i,-1],collapse = ",")
  output[i] <- gsub("NA,","",output[i])
}
system("mkdir /public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.fusion_results")
system("mkdir /public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.fusion_results/hla")

outPath <- "/public/home/lorihan/lrh/SKCM/Neoantigen/SKCM.fusion_results/hla" ##
out_fileName <- sapply(HLA.type[,1],function(x){
  paste(x, ".hla", sep='')}) ##
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ##
##
for(i in 1:nrow(HLA.type)){
  write.table(output[[i]], file=out_filePath[i], quote = F,row.names = F,col.names = F)
}
