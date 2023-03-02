#!/usr/local/bin/Rscript
myPaths <- .libPaths()  
new <- c('/public/home/lorihan/miniconda3/lib/R/library','/usr/local/lib64/R/library','/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1')#c('/public/home/lorihan/miniconda3/lib/R/library')#
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
library(stringr)
library(ABSOLUTE)
#setwd("/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE")
system('mkdir /public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE.new')
setwd("/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE.new")

# dat.SKCM <- read.delim('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE/focal_SKCM.input.seg.txt')
#dat.SKCM <- read.delim('/public/home/lorihan/lrh/SKCM/TCGA.SKCM.sampleMap%2FSNP6_nocnv_genomicSegment.gz')
dat.SKCM <- read.table('/public/home/lorihan/lrh/SKCM/ABSOLUTE.xena/TCGA.SKCM.hg38.SNP6_nocnv_genomicSegment',header = T)
str(dat.SKCM)
dat.SKCM$Sample<- as.factor(dat.SKCM$Sample)
#colnames(dat.SKCM) <- c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")
colnames(dat.SKCM) <- c("Sample","Chromosome","Start","End","Segment_Mean")

#dat.absolute <- list.files("/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE")
#dat.absolute <- dat.absolute[grep("solution.txt",dat.absolute)]
#table(substr(dat.absolute,1,28) %in% dat.SKCM$Sample)#516
#setdiff(dat.SKCM$Sample,substr(dat.absolute,1,28))#4
#diff <- setdiff(dat.SKCM$Sample,substr(dat.absolute,1,28))#269
#dat.diff <- dat.SKCM[dat.SKCM$Sample %in% diff,]
alpha <- NULL
for (i in levels(dat.SKCM$Sample)) {
  # ， 
  dat_i <- subset(dat.SKCM, Sample == i)
  sample_file <- paste(i, 'txt', sep = '.')
  write.table(dat_i, sample_file, sep = '\t', quote = FALSE, row.names = FALSE)
  # maf_file <- paste('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/MAF/',substr(i,1,16),'.maf',sep = '')
  maf_file <- paste('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/MAF/',substr(i,1,15),'*.maf',sep = '')
    RunAbsolute(
    seg.dat.fn = sample_file,    #"Chromosome、Start、End、Num_Probes、Segment_Mean"
    results.dir = i,    # ， 
    min.ploidy = 0.5,    # ploidy 
    max.ploidy = 10,     # ploidy  
    max.sigma.h = 0.2,     # 
    platform = 'Illumina_WES',    # 
    copy_num_type = 'total', # 
    sigma.p = 0,    # 
    primary.disease = 'SKCM',    #  
    sample.name = i,    # 
    max.as.seg.count = 30000,    # 
    max.non.clonal = 0,     # 
    maf.fn = maf_file,
    min.mut.af = 0.05,
    max.neg.genome = 0 )    # 
  load(paste(i, '/', i, '.ABSOLUTE.RData', sep = ''))
  alpha <- NULL
  alpha <- rbind(alpha, c({seg.dat$mode.res$mode.tab}[1,c("genome mass","alpha","tau")]))
  # ， 
  alpha <- data.frame(alpha, stringsAsFactors = FALSE)
  if (nrow(alpha)!=0) {
    names(alpha) <- c('Ploidy','tumorPurity','tumorPloidy')
    rownames(alpha) <- seg.dat$sample.name
    result=paste(i, 'solution.txt', sep = '_')
    write.table(alpha, result, row.names = TRUE, quote = FALSE, sep = '\t')
  }
}
