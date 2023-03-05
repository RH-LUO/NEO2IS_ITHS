# title: "SKCM PhenoGenotypes and modelling on neoantigens"
# author: "Ruihan Luo"
# date: "Aug 10th,2022"
.libPaths()  
myPaths <- .libPaths()  
new <- '/public/home/lorihan/miniconda3/lib/R/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
# rm(list = ls())
setwd("/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/")
load("TCGA-SKCM-summarize.RData")
# setwd("/public/home/lorihan/lrh/Neoantigen_TCGA/NSCLC")
options(scipen=1)
# load("TCGA-LUNG-summarize.RData")
ABSOLUTE.TCGA <- read.table("../TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                            sep="\t",fill = T,header = T)
#####---NEO2IS modelling---######
load('TCGA.LUNG.miceOutput.Rdata')
load('TCGA.SKCM.miceOutput.Rdata')
dat <- miceOutput
if (nrow(dat)>500) {
  dat$Age <- as.factor(ifelse(dat$Age>=65,"1","0"))
  dat$Sex <- as.factor(ifelse(dat$Sex=="Male","1","0"))
  dat$Stage <- factor(ifelse(dat$Stage=="I"|dat$Stage=="II","I/II","III/IV"))
} else{
  dat$Stage <- factor(ifelse(dat$Stage==1,"III/IV","0/I/II"))
  dat$Sex <- as.factor(dat$Sex)
}
table(dat$Stage)
dat$HLA_B2M.mut <- as.factor(dat$HLA_B2M.mut)
dat$T.stage <- as.factor(dat$T.stage)
dat$M.stage <- as.factor(dat$M.stage)
dat$N.stage <- as.factor(dat$N.stage)
shapiro.test(dat$CD8_Tcells)
table(is.na(dat$Purity))
dataCD <- dat
summary(dataCD$CD8.effector)
Count <- apply(dataCD,1,function(x){
  sum(as.numeric(x[c("No.Fusion","No.SNV","No.INDEL")]))
})
dataCD$Count.all <- Count
summary(dataCD$Count.all)#0-2154 (104,185,300)
var <- c("Score.Fusion.neo","Score.SNV.neo","Score.INDEL.neo")
library(factoextra)
library(ggplot2)
library(FactoMineR)
######
var <- c("No.Fusion", "No.INDEL", "No.SNV",
         "Score.Fusion","Score.SNV","Score.INDEL",
         "Score.Fusion.neo","Score.SNV.neo","Score.INDEL.neo")
dataCD[,var] <- apply(dataCD[,var], 2, function(x){
  log2(x+1)
})
###########################
set.seed(34245) 
dataCD$cytolytic <- apply(dataCD[,c("GZMA","PRF1")],1,function(x){
  mean(x)
})
library(tibble)
library(stringr)
if (nrow(dat) > 500) {
  load("/public/home/lorihan/lrh/NSCLC/ciber_perm1000_xena_ref_cluster.all.Rdata")
} else{
  load('/public/home/lorihan/lrh/SKCM/ciber_perm1000_GSE120575_cluster.all.new.Rdata')
}
dim(TME.results)
rownames(TME.results) <- gsub('[.]','-',rownames(TME.results))
rownames(TME.results) <- substr(rownames(TME.results),1,15)
table(rownames(TME.results) %in% dataCD$Samples)#T1007
dataCD <- dataCD[dataCD$Samples %in% rownames(TME.results),]
dim(dataCD)
if (nrow(TME.results)>500) {
  TME.results <- TME.results[rownames(TME.results) %in% Neoantigen.data_lung$Samples,]
} else{
  TME.results <- TME.results[rownames(TME.results) %in% Neoantigen.data_SKCM$Samples,]
}
library(dplyr)
library(tidyr)
colnames(TME.results)
dd1 <- TME.results %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:(ncol(TME.results)-2),
               names_to = "CellType",
               values_to = "Composition")
plot.info <- dd1[,c(5,1,6)]  
unique(plot.info$CellType)
plot.info$CellType <- str_split(plot.info$CellType,"- ",simplify = T)[,2]
plot.info$CellType <- gsub("CD8","CD8+",plot.info$CellType)
ylab <- expression('Proportion of exhausted CD8'^'+'*' T cell')
ylab <- expression('CD8'^'+'*' ')
#
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot","ggpubr","bslib","ggthemes")
lapply(pkgs, library, character.only = T)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",palette = pal,
  xlab = "",
  ylab = "CIBERSORTx-inferred cell fractions",
) + ggtitle("TCGA-LUNG cohort")+
  theme_base() +
  theme(legend.position='none',
        plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 14, angle = 40,
                                     hjust = 1,vjust = 1),
        axis.text.y   = element_text(color = 'black', size = 14, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 13.5, angle = 90),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x
  )
dat.Exh <- TME.results %>%
  as.data.frame() %>%
  rownames_to_column("Samples")
colnames(dat.Exh)
dat.Exh$Samples <- gsub("[.]","-",dat.Exh$Samples)
dat.Exh <- dat.Exh[dat.Exh$Samples %in% dataCD$Tumor_Sample_Barcode,]
dataCD$Exh_CD8 <- dat.Exh$`G11- Lymphocytes exhausted/cell-cycle`[match(dataCD$Samples,dat.Exh$Samples)]
descrCorr = dplyr::left_join(dataCD, dat.Exh[,1:(ncol(dat.Exh)-3)],by='Samples')

PCAdata <- descrCorr[,c('G11- Lymphocytes exhausted/cell-cycle',#'PDCD1',
                        'G9- Exhausted/HS CD8 T cells',
                        'G6- Exhausted CD8 T cells')]
rownames(PCAdata) <- dataCD$Samples
PCAdata[1:3,]
hla.pca <- PCA(PCAdata, graph = F)
get_eig(hla.pca)
ind <- get_pca_ind(hla.pca)
PC1 <- ind$coord[,1]
PC2 <- ind$coord[,2]
PC3 <- ind$coord[,3]
PC1[1:5]
corTest_PCA <- apply(PCAdata , 2, function(values){
  data.diff=cor.test(values,PC1,method = "spearman")
  cor = data.diff$estimate
  p.val = data.diff$p.value
  cor_DRR = t(c(cor,p.val))
})
rownames(corTest_PCA) = c("CC","P.Value")
corTest_PCA<-data.frame(t(corTest_PCA))
head(corTest_PCA)
dataCD$CD8Tex <-  (PC1[match(dataCD$Samples,names(PC1))])
######
LOH <- read.csv("/public/home/lorihan/lrh/Neoantigen_TCGA/TCGA_pan_HLA.LOH.csv",sep = "\t",header = T)
if (nrow(dataCD) >500) {
  Dat_LOH <- LOH[LOH$TCGA.Sample.ID %in% substr(Neoantigen.data_lung$Samples,1,12),]
} else{
  Dat_LOH <- LOH[LOH$TCGA.Sample.ID %in% substr(Neoantigen.data_SKCM$Samples,1,12),]
}
dataCD$LOH <- Dat_LOH$LOH.HLA[match(substr(dataCD$Samples,1,12),Dat_LOH$TCGA.Sample.ID)]
# table(dataCD$LOH)
table(dataCD$Age)
dataCD$LOH[is.na(dataCD$LOH)] <- "Unknown"
table(dataCD$LOH)
# -------
Dat_hla <- t(mut)[,c(grep("HLA-",rownames(mut)),grep("B2M",rownames(mut)))]
Dat_hla[1:5]
dim(Dat_hla)
if (length(unique(grepl("TFB2M",colnames(Dat_hla))))!=1) {
  Dat_hla <- Dat_hla[,-(grep("TFB2M",colnames(Dat_hla)))]
}
table(Dat_hla)
dataCD$B2M <- Dat_hla[match(dataCD$Samples,names(Dat_hla))]
dataCD$B2M <- ifelse(grepl("p.",dataCD$B2M),"Mut","Wt")
table(dataCD$B2M)
dat.new <-dataCD
dat.new[dat.new$B2M=="Mut",var] <- 0
# dataCD <- dat.new
dataCD.test <- dataCD
######
covariate <- c("CD80","TAM.effector","Str",
               'TAN',"CD4.effector",#'Th1'-->IL2-->CD8NaiveToCTL,'Th2'-->IL4-->BNaiveToBActivated
               'Tregs',
               "Purity",
               "Sex","Age",'HLA_B2M.mut'
)
var <- c("Score.Fusion.neo","Score.SNV.neo","Score.INDEL.neo"
)
covariates <- c(var,covariate)
covariates
length(covariates)#13
equation <- as.formula(paste0('CD8Tex~',paste0(covariates,collapse = "+")))
print(equation)
set.seed(321)
folds <- sample(nrow(dataCD),round(0.3*nrow(dataCD)))
fold_test<-dataCD[folds,]
fold_train<-dataCD[-folds,]
dim(fold_train)
dim(fold_test)
library(caret)
set.seed(123)
myControl <- trainControl(method = "cv",
                          number = 10,
                          verboseIter = T)
model_lin <- train(
  equation,
  fold_train,
  method = "lm",
  trControl = myControl)
model_svr <- train(
  equation,
  fold_train,
  method = "svmRadial",#"svmPoly"
  trControl = myControl)
model_gbm <- train(
  fold_train[,covariates], fold_train[,"CD8Tex"], 
  method = "gbm",#ntrees=ntrees,
  trControl = myControl)
summary(model_lin)
pred_lin <- predict(model_lin, fold_test[,covariates])
pred_svr <- predict(model_svr, fold_test[,covariates])
pred_gbm <- predict(model_gbm, fold_test[,covariates])

dim(fold_train)
dim(fold_test)
library(caret)
set.seed(123)
test <- train(
  equation.less,
  fold_train,
  method = "lm",
  trControl = myControl)
pred_test <- predict(test, fold_test[,covariates])
summary(test)
#####---feature selection using stepwise---#######
library(leaps)
library(MASS)
library(MASS, lib.loc = "/public/apps/R-4.1.0/library")
step.model.1 <- stepAIC(lm(equation, data=fold_train), direction = "both", 
                        trace = FALSE)
summary(step.model.1)
names(coef(step.model.1))
covariates <- gsub('Sex1','Sex',names(coef(step.model.1))[-1])
covariates <- gsub('HLA_B2M.mut1','HLA_B2M.mut',covariates)
equation <- as.formula(paste0('CD8Tex~',paste0(covariates,collapse = "+")))
print(equation)
########alternaltive models------
set.seed(123)
myControl <- trainControl(method = "cv",
                          number = 10,
                          verboseIter = T)
model_lin3 <- train(
  equation,
  fold_train,
  method = "lm",
  trControl = myControl)
summary(model_lin3)

model_svr3 <- train(
  equation,
  fold_train,
  method = "svmRadial",#"svmPoly"
  trControl = myControl)

model_gbm3 <- train(
  fold_train[,covariates], fold_train[,"CD8Tex"], 
  method = "gbm",#ntrees=ntrees,
  trControl = myControl)

summary(model_gbm3)
model_lin3$results$RMSE#
model_svr3$results$RMSE#
model_gbm3$results$RMSE#
pred_lin3 <- predict(model_lin3, fold_test[,covariates])
pred_svr3 <- predict(model_svr3, fold_test[,covariates])
pred_gbm3 <- predict(model_gbm3, fold_test[,covariates])

library(forecast)
#########---ME  RMSE   MAE   MPE MAPE---#########
cor.test(fold_test$CD8Tex,pred_lin3,method="spearman")#0.6556 
cor.test(fold_test$CD8Tex,pred_gbm3,method="spearman")#0.1505
cor.test(fold_test$CD8Tex,pred_svr3,method="spearman")#0.5861 
pdf("SKCM.fold_test_pred_lin3.pdf",width = 5,height = 5)
ggplot(data=fold_test, aes(pred_lin3,CD8Tex))+geom_point(colour = 'black', size = 1.5)+stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=as.data.frame(fold_test), method = "spearman",size=5,
           label.x = -4, label.y = 0.9)+theme(#panel.grid=element_blank(),
             axis.ticks=element_line(color='black'),
             axis.text=element_text(colour="black",size=13),
             # panel.background = element_rect(fill='transparent'),
             #panel.border=element_rect(fill='transparent'),
             panel.border=element_rect(fill='transparent',color='black'),
             axis.line=element_line(color='black',size=0.2),
             axis.title=element_text(color='black',size=15),#FC4E07","#00AFBB
             plot.title = element_text(size=14,hjust = 0.5))+
  #ggtitle( "Internal Validation (TCGA)" ) + 
  labs(x="Predicted value (linear model)",y="CD8-Tex")
dev.off()
pdf("SKCM.fold_test_pred_gbm3.pdf",width = 5,height = 5)
ggplot(data=fold_test, aes(pred_gbm3,CD8Tex))+geom_point(colour = 'black', size = 1.5)+stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=as.data.frame(fold_test), method = "spearman",size=5,
           label.x = -3.7, label.y = 0.98)+theme(#panel.grid=element_blank(),
             axis.ticks=element_line(color='black'),
             axis.text=element_text(colour="black",size=13),
             # panel.background = element_rect(fill='transparent'),
             #panel.border=element_rect(fill='transparent'),
             panel.border=element_rect(fill='transparent',color='black'),
             axis.line=element_line(color='black',size=0.2),
             axis.title=element_text(color='black',size=15),#FC4E07","#00AFBB
             plot.title = element_text(size=14,hjust = 0.5))+
  #ggtitle( "Internal Validation (TCGA)" ) + 
  labs(x="Predicted value (GBM model)",y="CD8-Tex")
dev.off()
pdf("SKCM.fold_test_pred_svr3.pdf",width = 5,height = 5)
ggplot(data=fold_test, aes(pred_svr3,CD8Tex))+geom_point(colour = 'black', size = 1.5)+stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=as.data.frame(fold_test), method = "spearman",size=5,
           label.x = -3.8, label.y = 0.98)+theme(#panel.grid=element_blank(),
             axis.ticks=element_line(color='black'),
             axis.text=element_text(colour="black",size=13),
             # axis.text.x=element_text(colour="black",size=15),
             # axis.text.y=element_text(colour="black",size=15),
             # panel.background = element_rect(fill='transparent'),
             #panel.border=element_rect(fill='transparent'),
             panel.border=element_rect(fill='transparent',color='black'),
             axis.line=element_line(color='black',size=0.2),
             axis.title=element_text(color='black',size=15),#FC4E07","#00AFBB
             plot.title = element_text(size=14,hjust = 0.5))+
  #ggtitle( "Internal Validation (TCGA)" ) +
  labs(x="Predicted value (SVR model)",y="CD8-Tex")
dev.off()

library(coin)
N <- nrow(dat.SKCM)
# Use custom color palette
# Add jitter points and change the shape by time
SKCMDf <- data.frame(id=factor(rep(dat.SKCM$Samples, times=3)),
                     Neo=c(dat.SKCM$No.SNV.mut, dat.SKCM$No.Fusion.mut, 
                           dat.SKCM$No.INDEL.mut),
                     gr=factor(rep(0:2, each=N), labels=c("SNV (/mut)","Fusion (/mut)", "Indel (/mut)")))
# Add p-values comparing groups
# Specify the comparisons you want
compar<-list(c("SNV (/mut)", "Indel (/mut)"),
             c("SNV (/mut)", "Fusion (/mut)"),c("Indel (/mut)", "Fusion (/mut)"))
compare_means(Neo~gr, data=SKCMDf,method="wilcox.test",p.adjust.method="BH")

# add boxplot with white fill color
ylab <- "Neoantigen load"
p <- ggboxplot(SKCMDf, x = "gr", y = "Neo",
               color = "gr", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "gr")+xlab("")+ylab(ylab)+ ggtitle("TCGA SKCM cohort") 

pdf("TCGA_SKCM_SNV_INDEL_Fusion.diff_boxplot.pdf",width=5.5,height = 5)
p + #stat_compare_means(comparisons = compar,label.x = 1.5, label.y = c(10,6,8),
  #                        method = "wilcox.test") +
  theme(legend.position='none')+
  theme(
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    axis.title.y  = element_text(color = 'black', size = 14, angle = 90),
    legend.title  = element_text(color = 'black', size  = 16),
    legend.text   = element_text(color = 'black', size   = 16),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
  )
dev.off()
model <- lm(equation, data=fold_train)
cor<-summary(model_lin3)
coef <- round(coef(cor)[,1],4)
tvalue<-round(cor$coefficients[,3],4)
Pvalue<-round(cor$coefficients[,4],4)
SE <- round(cor$coefficients[,2],4)
OR<-round(exp(coef(model)),4)
OR<-round(exp(coef),4)
rows<-rownames(cor$coefficients)
CIL<-paste0(round(exp(confint(model))[,1],2))
CIU<-paste0(round(exp(confint(model))[,2],2))
CI<-paste0(CIL,'-',CIU)                          
logit_result<-data.frame('characteristics'=rows,
                         'Coef' = coef,#coef(cor),#
                         'SE' = SE,
                         't value' = tvalue,
                         'Odds Ratio'=OR,
                         'CI95'=CI,
                         'P Value'=Pvalue)
logit_result.lin3 <- logit_result
write.table(logit_result.lin3,"logit_result_all.txt",quote = F,row.names = F)

dataCD.new <- Neoantigen.data_SKCM
var <- c("Score.Fusion","Score.SNV","Score.INDEL")
var <- c("Score.Fusion.neo","Score.SNV.neo","Score.INDEL.neo",
         "No.Fusion","No.SNV","No.INDEL")
dataCD.new[,var] <- apply(dataCD.new[,var], 2, function(x){
  log2(x+1)
})
logit_result.3 <- logit_result.lin3[logit_result.lin3$characteristics %in% var,]
logit_result.3 <- logit_result.3[logit_result.3$P.Value<0.05,]
sigExpressionProfile <- dataCD.new[,logit_result.3$characteristics]
coef <- logit_result.3$Coef
Neoantigen.model <- data.matrix(sigExpressionProfile) %*% coef
names(Neoantigen.model) <- dataCD.new$Samples
Neoantigen.data_SKCM$Neoantigen.model <- as.numeric(Neoantigen.model)#data.matrix(sigExpressionProfile) %*% coef
# Neoantigen.model is NEO2IS 
table(is.na(Neoantigen.data_SKCM$Neoantigen.model))
summary(Neoantigen.model)
Neoantigen.data_SKCM$B2M <- dataCD$B2M[match(Neoantigen.data_SKCM$Samples,
                                             dataCD$Samples)]
Neoantigen.data_SKCM$CD8Tex <- dataCD$CD8Tex[match(Neoantigen.data_SKCM$Samples,
                                                   dataCD$Samples)]


#########---NEOITH score calculation---##########
options("scipen"=100, "digits"=4)
library(stringr)
library(ABSOLUTE)
setwd("/public/home/lorihan/lrh/SKCM/ABSOLUTE")
dat.skcm <- read.delim('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/ABSOLUTE/focal_SKCM.input.seg.txt')
colnames(dat.skcm) <- c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")
setwd("/public/home/lorihan/lrh/SKCM")
dat.SKCM <- read.delim('/public/home/lorihan/lrh/SKCM/TCGA.SKCM.sampleMap%2FSNP6_nocnv_genomicSegment.gz')
table(unique(dat.SKCM$sample) %in% substr(unique(dat.cn$Sample),1,15))
dat.skcm.19 <- dat.SKCM[,c(2:4)]
system("mkdir /public/home/lorihan/lrh/SKCM/ABSOLUTE.xena")
setwd("/public/home/lorihan/lrh/SKCM/ABSOLUTE.xena")
write.table(dat.skcm.19,"skcm_hg19.seg.bed",col.names = F,row.names = F,quote = F)
system("/public/home/lorihan/bin/liftOver skcm_hg19.seg.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz skcm_hg38.bed skcm_hg38_unMapped.bed")

dat.skcm.38 <- read.delim("skcm_hg38.bed",header = F)
hg38_unmapped <- read.table("skcm_hg38_unMapped.bed")
unique(hg38_unmapped$V2)
unique(hg38_unmapped$V3)
hg38_unmapped.data <- hg38_unmapped$V3
hg38_unmapped.data <- hg38_unmapped$V2

hg38_unmapped.data <- paste(hg38_unmapped$V1,hg38_unmapped$V3,sep = ":")
tmp.unmapped <- lapply(hg38_unmapped.data,function(x){  
  tmp.unmapped <- grep(x,paste(dat.skcm.19$chr,dat.skcm.19$end,sep = ":"))
})
count <- lapply(tmp.unmapped,function(x){  
  count <- length(x)
})
dup <- which(count==2)
tmp.unmapped[dup]
dim(dat.skcm.19)
dup.19 <- dat.skcm.19[unlist(tmp.unmapped[dup]),]
rownames(dup.19)[dup.19$start %in% hg38_unmapped$V2]
table(unlist(count))
hg38_unmapped.data[3]
tmp.unmapped[[3]]
dat.skcm.19[119,]
hg19_unmapped <- unique(unlist(tmp.unmapped))
hg19_unmapped <- hg19_unmapped[!hg19_unmapped %in% rownames(dup.19)[!dup.19$start %in% hg38_unmapped$V2]]
dat.skcm.19 <- dat.skcm.19[-hg19_unmapped,]
write.table(dat.skcm.19,"skcm_hg19.seg.new.bed",col.names = F,row.names = F,quote = F)
system("/public/home/lorihan/bin/liftOver skcm_hg19.seg.new.bed  /public/home/lorihan/lrh/LUAD/fusion/hg19ToHg38.over.chain.gz skcm_hg38_new.bed skcm_hg38_unMapped_new.bed")

dat.skcm.hg38 <- read.delim("skcm_hg38_new.bed",header = F)
head(dat.SKCM)
head(dat.skcm.19)
dat.SKCM <- dat.SKCM[-hg19_unmapped,]
dat.SKCM$start <- dat.skcm.hg38$V2
dat.SKCM$end <- dat.skcm.hg38$V3
head(dat.skcm.19)
head(dat.SKCM)
colnames(dat.SKCM) <- c("Sample","Chromosome","Start","End","Segment_Mean")
write.table(dat.SKCM,"TCGA.SKCM.hg38.SNP6_nocnv_genomicSegment",col.names = T,row.names = F,quote = F)
# -------
setwd("/public/home/lorihan/lrh/SKCM/ABSOLUTE.xena")
dat.cn <- read.table('TCGA.SKCM.hg38.SNP6_nocnv_genomicSegment',header = T)
table(substr(unique(dat.skcm$Sample),1,15) %in% unique(dat.cn$Sample))
dat.cn$Sample <- as.factor(dat.cn$Sample)
for (i in levels(dat.cn$Sample)) {
  #按样本，拆分成独立文件
  dat_i <- subset(dat.cn, Sample == i)
  sample_file <- paste(i, 'seg', sep = '.')
  write.table(dat_i, sample_file, sep = '\t', quote = FALSE, row.names = FALSE)
}
setdiff(Neoantigen_SKCM$Barcode,substr(levels(dat.cn$Sample),1,15))#104-->2
########
Tumor_Sample_Barcode <- unique(maf$Tumor_Sample_Barcode)
SKCM_samples <- unique(dat.cn$Sample)#unique(c(dat.lusc$Sample,dat.luad$Sample))
Samples <- substr(SKCM_samples,1,15)
setdiff(Samples, Tumor_Sample_Barcode)#"TCGA-ER-A19T-06" "TCGA-ER-A2NF-06"
#-->5
setdiff(Tumor_Sample_Barcode,Samples)#104
#-->2
Samples <- intersect(Samples,Tumor_Sample_Barcode)#365
setdiff(Neoantigen.data_SKCM$Samples,substr(SKCM_samples,1,15))#2
# "TCGA-XV-AB01-06" "TCGA-WE-A8ZR-06"
SKCM_samples <- SKCM_samples[substr(SKCM_samples,1,15) %in% Samples]
# write.table(SKCM_samples,col.names = F, row.names = F,quote = F,
#             "/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/TCGA_samples_maf_cn.txt")
# write.table(SKCM_samples,col.names = F, row.names = F,quote = F,
#             "/public/home/lorihan/lrh/SKCM/TCGA_samples_maf_cn.txt")
#######
path="/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/CCF.new"
fileNames <- dir(path) 
fileNames  <- fileNames[grep(".tsv",fileNames)]
filePath.abridged <- sapply(fileNames,function(x){
  paste(x,sep = "/")
}) ##csv
matched <- gsub(".tsv","",filePath.abridged)
table(substr(matched,1,15) %in% unique(Neoantigen.data_SKCM$Sample))
matched[!substr(matched,1,15) %in% unique(Neoantigen.data_SKCM$Sample)]
filePath  <- sapply(filePath.abridged, function(x){
  paste(path,x,sep='/')}) ##output path
data_maf.SKCM <- lapply(filePath, function(x){
  if (file.info(x)$size!=0) {
    read.csv(x,header = T,sep = "\t",fill = T)}})  
data_maf.SKCM <- do.call('rbind',data_maf.SKCM)
Tumor_Sample_Barcode <-  str_split(data_maf.SKCM$mutation_id,"[:]",simplify = T)[,2]
Tumor_Sample_Barcode <- str_split(rownames(data_maf.SKCM),".tsv",simplify = T)[,1]
Tumor_Sample_Barcode <- substr(Tumor_Sample_Barcode,1,15)
length(unique(Tumor_Sample_Barcode))#365-->467
#------
Neoantigen_SKCM$Barcode <- substr(Neoantigen_SKCM$Barcode,1,15)
setdiff(unique(Tumor_Sample_Barcode), Neoantigen_SKCM$Barcode)#"TCGA-EB-A5VV-06"
setdiff(Neoantigen_SKCM$Barcode,unique(Tumor_Sample_Barcode))#104-->2
Symbol <-  str_split(data_maf.SKCM$mutation_id,"[:]",simplify = T)[,1]
#vaf <- data_maf.SKCM
data_maf.SKCM$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
data_maf.SKCM$Score <- data_maf.SKCM$var_counts/(data_maf.SKCM$ref_counts+data_maf.SKCM$var_counts)
table(is.na(data_maf.SKCM$all_effects))
table(data_maf.SKCM$all_effects %in% maf$all_effects)
table( unique(maf$Tumor_Sample_Barcode) %in% data_maf.SKCM$Tumor_Sample_Barcode)
#####
vaf <- maf[,c("Tumor_Sample_Barcode","Score","Hugo_Symbol","Chromosome",
              "Start_Position","End_Position","HGVSp_Short","all_effects")]
vaf$Tumor_Sample_Barcode <- substr(vaf$Tumor_Sample_Barcode,1,15)
table( vaf$all_effects %in% data_maf.SKCM$all_effects)
table(unique(vaf$Tumor_Sample_Barcode) %in% data_maf.SKCM$Tumor_Sample_Barcode)
data_maf.SKCM$Mutation <- paste(data_maf.SKCM$Tumor_Sample_Barcode,
                                data_maf.SKCM$all_effects, sep = ":")
setdiff(unique(vaf$Tumor_Sample_Barcode),
        unique(data_maf.SKCM$Tumor_Sample_Barcode))###104-->2
Mutation <- paste(vaf$Tumor_Sample_Barcode,
                  vaf$all_effects, sep = ":")
table( Mutation %in% data_maf.SKCM$Mutation)
table(data_maf.SKCM$Mutation %in% paste(vaf$Tumor_Sample_Barcode,vaf$all_effects, sep = ":"))
table(is.na(match(Mutation,data_maf.SKCM$Mutation)))
#F 536480 T 392067
#F 535109 391084
vaf$CN <- data_maf.SKCM$major_cn[match(Mutation,data_maf.SKCM$Mutation)]
table(is.na(vaf$CN))
length(unique(data_maf.SKCM$Tumor_Sample_Barcode))#365-->467
#########---first filling---#########
vaf$Mutation <- paste(vaf$Tumor_Sample_Barcode,vaf$Hugo_Symbol,
                      str_split(vaf$HGVSp_Short,"[*]",simplify = T)[,1]
                      , sep = " ")
Mutation <- paste(substr(Neoantigen_SKCM$Barcode, 1,15),
                  str_split(Neoantigen_SKCM$Mutation,"[*]",simplify = T)[,1], sep = " ")
table(is.na(match(Mutation,vaf$Mutation)))
#F 247078  T 15563
Neoantigen_SKCM$VAF <- vaf$Score[match(Mutation,vaf$Mutation)]
Neoantigen_SKCM$CN <- vaf$CN[match(Mutation,vaf$Mutation)]
table(is.na(Neoantigen_SKCM$VAF))#1650-->1407           -->15563
table(is.na(Neoantigen_SKCM$CN))#47912 -->8416 -->108301-->116781
table(is.na(Neoantigen_SKCM$VAF[-grep("p.",Neoantigen_SKCM$Mutation)]))
#T 3500
Neoantigen_SKCM$VAF[-grep("p.",Mutation)] <- 1
Neoantigen_SKCM$CN[-grep("p.",Neoantigen_SKCM$Mutation)] <- 2
table(is.na(Neoantigen_SKCM$VAF))#12063
table(is.na(Neoantigen_SKCM$CN))#113281 
#---second filling----
cn.vaf <- vaf#maf#vaf#data_maf.SKCM
cn.missed <- Neoantigen_SKCM[is.na(Neoantigen_SKCM$VAF),]#12063    11
#cn.missed <- cn.missed[grep("p.",cn.missed$Mutation),]
gene <- str_split(cn.missed$Mutation," ",simplify = T)[,1]
table(unique(gene) %in% vaf$Hugo_Symbol)
###round 1-----
load("/public/home/lorihan/lrh/NSCLC/gene_info.RData")
cn.missed.1 <- cn.missed[!gene %in% cn.vaf$Hugo_Symbol,]
cn.vaf$Mutation <- paste(cn.vaf$Tumor_Sample_Barcode,#cn.vaf$Hugo_Symbol
                         cn.vaf$HGVSp_Short
                         ,sep = " ")
Mutation <- paste(cn.missed.1$Barcode,#cn.missed$Mutation
                  str_split(cn.missed.1$Mutation," ",simplify = T)[,2]
                  ,sep = " ")
match <- match(Mutation,cn.vaf$Mutation)
table(is.na(match))#F9076(5576)    T89
cn.missed.1$type <- paste(cn.vaf$Hugo_Symbol[match],
                          cn.vaf$HGVSp_Short[match],sep = " ")
cn.missed.1$VAF <- cn.vaf$Score[match]
cn.missed.1$CN <- cn.vaf$CN[match]
cn.missed.2 <- cn.missed.1[!duplicated(cn.missed.1$Mutation),]
variant.logi <- 1:nrow(cn.missed.2)
variant.logi <- vector("list",nrow(cn.missed.2))
variant.logi <- apply(cn.missed.2, 1, function(x){
  gene1 <- str_split(x["Mutation"]," ",simplify = T)[,1]
  gene2 <- str_split(x["type"]," ",simplify = T)[,1]
  logi <- grep(gene2,gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
  return(identical(logi, integer(0)))
})
table(variant.logi)#F2110   T289
cn.missed.2$logi <-  !variant.logi
table(cn.missed.2$logi)
cn.missed.1$logi <- cn.missed.2$logi[match(cn.missed.1$Mutation,cn.missed.2$Mutation)]
cn.missed.1$VAF[cn.missed.1$logi=="FALSE"] <- NA
cn.missed.1$CN[cn.missed.1$logi=="FALSE"] <- NA
#########---round 2---########
cn.missed.3 <- cn.missed[gene %in% vaf$Hugo_Symbol,]
cn.missed.2 <- cn.missed.3[!duplicated(cn.missed.3$Mutation),]
variant.missed <- apply(cn.missed.2, 1, function(x){
  sample <- x["Barcode"]
  variant <- str_split(x["Mutation"]," ",simplify = T)[,1]
  aa <- str_split(x["Mutation"]," ",simplify = T)[,2]
  if (grepl("_",aa)) {
    aa <- str_split(aa,"_",simplify = T)[,1]
    bb <- substring(aa,4,nchar(aa))  
  }else{
    bb <- substring(aa,4,nchar(aa)-1)  
  }
  gene.row <- (variant==cn.vaf$Hugo_Symbol
               & sample==cn.vaf$Tumor_Sample_Barcode)
  variant.row <- which(gene.row)[grepl(bb,as.vector(toupper(
    cn.vaf[gene.row,]$all_effects)))|grepl(bb,as.vector(toupper(
      cn.vaf[gene.row,]$HGVSp_Short)))]
  if(length(variant.row)>1){
    bb <- substring(aa,4,nchar(aa))
    variant.row <- which(gene.row)[grepl(bb,as.vector(toupper(
      cn.vaf[gene.row,]$all_effects)))|grepl(bb,as.vector(toupper(
        cn.vaf[gene.row,]$HGVSp_Short)))]
  }
  if(length(variant.row)>1){
    gene1 <- variant
    variant.logi <- apply(cn.vaf[variant.row,], 1, function(z){
      gene2 <- z["Hugo_Symbol"]
      logi <- grep(gene2,gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
      return(identical(logi, integer(0)))
    })
    variant.row <- variant.row[!variant.logi]
  }
  if(length(variant.row)>1){
    bb <- substring(aa,4,nchar(aa))
    variant.row <- variant.row[which.max(cn.vaf$Score[variant.row])]
  }
  return(variant.row)
})
variant.row <- unlist(variant.missed)
cn.missed.2$VAF <- cn.vaf$Score[variant.row[
  match(rownames(cn.missed.2),names(variant.row))]]
cn.missed.2$CN <- cn.vaf$CN[variant.row[
  match(rownames(cn.missed.2),names(variant.row))]]
cn.missed.3$VAF <- cn.missed.2$VAF[match(cn.missed.3$Mutation,cn.missed.2$Mutation)]
cn.missed.3$CN <- cn.missed.2$CN[match(cn.missed.3$Mutation,cn.missed.2$Mutation)]

cn.missed.2 <- rbind(cn.missed.1[,-ncol(cn.missed.1)],cn.missed.3)
table(is.na(cn.missed.1$VAF))#F4979  686
table(is.na(cn.missed.3$VAF))#4815 1583
# setwd('/public/home/lorihan/lrh/Neoantigen_TCGA/')
# save(cn.missed.1,cn.missed.2,cn.missed.3,file="SKCM.cn.missed.Rdata")
load('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM.cn.missed.Rdata')
table(is.na(cn.missed.2$VAF))#F9794 2269
dim(cn.missed)
dim(cn.missed.1)#
dim(cn.missed.3)
missed <- is.na(Neoantigen_SKCM$VAF)
Neoantigen_SKCM$CN[missed][!gene %in% cn.vaf$Hugo_Symbol] <- cn.missed.1$CN
Neoantigen_SKCM$VAF[missed][!gene %in% cn.vaf$Hugo_Symbol] <- cn.missed.1$VAF
Neoantigen_SKCM$CN[missed][gene %in% cn.vaf$Hugo_Symbol] <- cn.missed.3$CN
Neoantigen_SKCM$VAF[missed][gene %in% cn.vaf$Hugo_Symbol] <- cn.missed.3$VAF
table(is.na(Neoantigen_SKCM$VAF))# T2269
Neoantigen_SKCM.new <- Neoantigen_SKCM
#########---Third filling---#########
cn.vaf <- vaf
cn.missed <- Neoantigen_SKCM[is.na(Neoantigen_SKCM$VAF),]
cn.missed.2 <- cn.missed[!duplicated(cn.missed$Mutation),]
dim(cn.missed.2)#944
variant.missed <- apply(cn.missed.2, 1, function(x){
  sample <- x["Barcode"]
  variant <- str_split(x["Mutation"]," ",simplify = T)[,1]
  aa <- str_split(x["Mutation"]," ",simplify = T)[,2]
  if (grepl("_",aa)) {
    aa <- str_split(aa,"_",simplify = T)[,1]
    bb <- substring(aa,4,nchar(aa))  
  }else{
    bb <- substring(aa,4,nchar(aa)-1)  
  }
  gene.row <- (#variant==cn.vaf$Hugo_Symbol&
    sample==cn.vaf$Tumor_Sample_Barcode)
  variant.row <- which(gene.row)[grepl(bb,as.vector(toupper(
    cn.vaf[gene.row,]$all_effects)))|grepl(bb,as.vector(toupper(
      cn.vaf[gene.row,]$HGVSp_Short)))]
  if(length(variant.row)>1){
    gene1 <- variant
    variant.logi <- apply(cn.vaf[variant.row,], 1, function(z){
      gene2 <- z["Hugo_Symbol"]
      logi <- grep(gene2,gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
      return(identical(logi, integer(0)))
    })
    variant.row <- variant.row[!variant.logi]
  }
  if(length(variant.row)>1){
    variant.row <- variant.row[which.max(cn.vaf$Score[variant.row])]
  }
  return(variant.row)
})
variant.row <- unlist(variant.missed)
match <- match(rownames(cn.missed.2),names(variant.row))
table(is.na(match))#F371   573
cn.missed.2$VAF <- cn.vaf$Score[variant.row[match]]
cn.missed.2$CN <- cn.vaf$CN[variant.row[match]]
cn.missed.2$type <- paste(cn.vaf$Hugo_Symbol[variant.row[match]],
                          cn.vaf$HGVSp_Short[variant.row[match]],sep = " ")
cn.missed.2$type[is.na(match)] <- NA
cn.missed$VAF <- cn.missed.2$VAF[match(cn.missed$Mutation,cn.missed.2$Mutation)]
cn.missed$CN <- cn.missed.2$CN[match(cn.missed$Mutation,cn.missed.2$Mutation)]
# ----
load('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM.vaf.missed.Rdata')
Neoantigen_SKCM <- Neoantigen_SKCM.new
table(is.na(cn.missed$VAF))
missed <- is.na(Neoantigen_SKCM$VAF)
table(missed)
Neoantigen_SKCM$CN[missed] <- cn.missed$CN
Neoantigen_SKCM$VAF[missed]<- cn.missed$VAF
table(is.na(Neoantigen_SKCM$VAF))#1387
# save(missed,cn.missed.2,cn.missed,Neoantigen_SKCM,file="SKCM.vaf.missed.Rdata")
#######
variant.logi <- 1:nrow(cn.missed.2)
variant.logi <- vector("list",nrow(cn.missed.2))
variant.logi <- apply(cn.missed.2, 1, function(x){
  gene1 <- str_split(x["Mutation"]," ",simplify = T)[,1]
  gene2 <- str_split(x["type"]," ",simplify = T)[,1]
  logi <- grep(gene2,gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
  return(identical(logi, integer(0)))
})
table(variant.logi)#ALL FALSE
cn.missed.2$logi <-  !variant.logi
table(cn.missed.2$logi)#944
cn.missed$logi <- cn.missed.2$logi[match(cn.missed$Mutation,cn.missed.2$Mutation)]
table(cn.missed$logi)#2269:all TRUE
table(is.na(Neoantigen_SKCM$VAF))#1387
####final determine----
cn.missed <- Neoantigen_SKCM[is.na(Neoantigen_SKCM$VAF),]
cn.missed.2 <- cn.missed[!duplicated(cn.missed$Mutation),]
dim(cn.missed.2)#573
variant.missed <- apply(cn.missed.2, 1, function(x){
  sample <- as.character(x["Barcode"])
  variant <- str_split(x["Mutation"]," ",simplify = T)[,1]
  aa <- str_split(x["Mutation"]," ",simplify = T)[,2]
  if (grepl("_",aa)) {
    aa <- str_split(aa,"_",simplify = T)[,1]
    bb <- substring(aa,4,nchar(aa))  
  }else{
    bb <- substring(aa,4,nchar(aa))  
  }
  gene.row <- (#variant==cn.vaf$Hugo_Symbol&
    sample==cn.vaf$Tumor_Sample_Barcode)
  variant.row <- which(gene.row)[grepl(bb,as.vector(toupper(
    cn.vaf[gene.row,]$all_effects)))|grepl(bb,as.vector(toupper(
      cn.vaf[gene.row,]$HGVSp_Short)))]
  if(length(variant.row)>1){
    gene1 <- variant
    variant.logi <- apply(cn.vaf[variant.row,], 1, function(z){
      gene2 <- as.character(z["all_effects"])
      gene2 <- unlist(str_split(gene2,',;'))
      gene2 <- gene2[gene2!=""]
      gene2 <- unique(str_split(gene2,',',simplify = T)[,1])
      gene2 <- gene2[gene2!=""]
      logi <- rep("FALSE",length(gene2))
      for (i in 1:length(gene2)) {
        gene.logi <- grep(gene2[i],gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
        logi[i] <- identical(gene.logi, integer(0))
      }
      logi <- unique(logi)
      if (length(logi)==1) {
        return(logi=="FALSE")
      } else{
        return(length(logi)==2)}
    })
    variant.row <- variant.row[variant.logi]
  }
  return(variant.row)
})
#save(cn.missed,cn.missed.2,variant.missed,file="SKCM.variant.missed.Rdata")
load('/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM.variant.missed.Rdata')
variant.row <- unlist(variant.missed)
match <- match(rownames(cn.missed.2),names(variant.row))
table(is.na(match))#F2   T9 -->324   249
cn.missed.2$VAF <- cn.vaf$Score[variant.row[match]]
cn.missed.2$CN <- cn.vaf$CN[variant.row[match]]
table(is.na(cn.missed.2$VAF))
cn.missed.2$type <- paste(cn.vaf$Hugo_Symbol[variant.row[match]],
                          cn.vaf$HGVSp_Short[variant.row[match]],sep = " ")
cn.missed.2$type[is.na(match)] <- NA
cn.missed$VAF <- cn.missed.2$VAF[match(cn.missed$Mutation,cn.missed.2$Mutation)]
cn.missed$CN <- cn.missed.2$CN[match(cn.missed$Mutation,cn.missed.2$Mutation)]
table(is.na(cn.missed$VAF))
missed <- is.na(Neoantigen_SKCM$VAF)
table(missed)#1387 
Neoantigen_SKCM$CN[missed] <- cn.missed$CN
Neoantigen_SKCM$VAF[missed]<- cn.missed$VAF
table(is.na(Neoantigen_SKCM$VAF))#605
cn.missed.3 <- cn.missed
missed.3 <- missed

cn.vaf <- vaf
cn.missed <- Neoantigen_SKCM[is.na(Neoantigen_SKCM$VAF),]
cn.vaf$Mutation <- paste(cn.vaf$Tumor_Sample_Barcode,#cn.vaf$Hugo_Symbol,cn.vaf$HGVSp_Short
                         str_split(cn.vaf$HGVSp_Short,"[*]",simplify = T)[,1]
                         ,sep = " ")
Mutation <- paste(cn.missed$Barcode,
                  #str_split(cn.missed$Mutation,"[ ]",simplify = T)[,2]
                  #cn.missed$Mutation
                  str_split(str_split(cn.missed$Mutation,"[ ]",simplify = T)[,2],"[*]",simplify = T)[,1]
                  ,sep = " ")
table(is.na(match(Mutation,cn.vaf$Mutation)))#F586 T19
match <- match(Mutation,cn.vaf$Mutation)
cn.missed$type <- paste(cn.vaf$Hugo_Symbol[match],
                        cn.vaf$HGVSp_Short[match],sep = " ")
cn.missed$CN <- cn.vaf$CN[match(Mutation,cn.vaf$Mutation)]
cn.missed$VAF <- cn.vaf$Score[match(Mutation,cn.vaf$Mutation)]
table(is.na(cn.missed$VAF))#45-->19
table(is.na(cn.missed$CN))#709-->299
missed  <- is.na(Neoantigen_SKCM$VAF)
Neoantigen_SKCM[missed,]$CN <- cn.missed$CN#[match]
Neoantigen_SKCM[missed,]$VAF <- cn.missed$VAF#[match]

variant.logi <- 1:nrow(cn.missed.2)
variant.logi <- vector("list",nrow(cn.missed.2))
variant.logi <- apply(cn.missed.2, 1, function(x){
  gene1 <- str_split(x["Mutation"]," ",simplify = T)[,1]
  gene2 <- str_split(x["type"]," ",simplify = T)[,1]
  logi <- grep(gene2,gene_info$geneAlias[grep(gene1,gene_info$geneAlias)])
  return(identical(logi, integer(0)))
})
table(variant.logi)#F8   565-->F249   324
cn.missed.2$logi <-  !variant.logi
table(cn.missed.2$logi)#F565     8-->F324   249
cn.missed.3$logi <- cn.missed.2$logi[match(cn.missed.3$Mutation,cn.missed.2$Mutation)]
table(cn.missed.3$logi)#T36-->605
cn.missed.3$VAF[cn.missed.3$logi=="FALSE"] <- NA
cn.missed.3$CN[cn.missed.3$logi=="FALSE"] <- NA

# Neoantigen_SKCM.old <- Neoantigen_SKCM
Neoantigen_SKCM <- Neoantigen_SKCM.old
table(Neoantigen_SKCM$Neoantigen %in% Neoantigen_SKCM.p5$Neoantigen)
logi.p5 <- paste(Neoantigen_SKCM$Barcode,Neoantigen_SKCM$Neoantigen,sep = ":") %in% paste(substr(Neoantigen_SKCM.p5$Barcode,1,15),Neoantigen_SKCM.p5$Neoantigen,sep = ":")
table(logi.p5)#131899 130742
Neoantigen_SKCM <- Neoantigen_SKCM[logi.p5,]
table(is.na(Neoantigen_SKCM$VAF))#10
table(is.na(Neoantigen_SKCM$CN))#52244
####
table(is.na(match(substr(Neoantigen_SKCM$Barcode, 1,15),
                  Neoantigen.data_SKCM$Samples)))
table(is.na(match(substr(Neoantigen_SKCM$Barcode, 1,15),
                  miceOutput_skcm$Samples)))

Neoantigen_SKCM$purity <- miceOutput_skcm$Purity[
  match(substr(Neoantigen_SKCM$Barcode,1,15),miceOutput_skcm$Samples)]
Neoantigen_SKCM$ploidy <- miceOutput_skcm$Ploidy[
  match(substr(Neoantigen_SKCM$Barcode,1,15),miceOutput_skcm$Samples)]
Purity.na.samples <- Neoantigen_SKCM$Barcode[is.na(Neoantigen_SKCM$purity)]
table(is.na(match(substr(Purity.na.samples, 1,15), miceOutput$Samples)))
unique(Purity.na.samples)#1:TCGA-XV-AB01-06#TCGA-EE-A3AE-06-->11
ABSOLUTE.SKCM$purity[match(unique(Purity.na.samples),ABSOLUTE.SKCM$array)]
Neoantigen_SKCM$purity[is.na(Neoantigen_SKCM$purity)] <- tumorPurity$tumorPurity[match(Purity.na.samples,substr(rownames(tumorPurity),1,15))]
Neoantigen_SKCM$ploidy[is.na(Neoantigen_SKCM$ploidy)] <- tumorPurity$tumorPloidy[match(Purity.na.samples,substr(rownames(tumorPurity),1,15))]
table(is.na(Neoantigen_SKCM$purity))#0-->10
table(is.na(Neoantigen_SKCM$CN))#F156206 106435
unique(Neoantigen_SKCM$Barcode[is.na(Neoantigen_SKCM$CN)])

Neoantigen_SKCM.all <- Neoantigen_SKCM[grep(' p.',Neoantigen_SKCM$Mutation),]
Neoantigen_SKCM.uniq <- Neoantigen_SKCM.all[!duplicated(Neoantigen_SKCM.all$Mutation),]
summary(Neoantigen_SKCM.uniq$VAF)
summary(Neoantigen_SKCM$CN)
Neoantigen_SKCM$CN[is.na(Neoantigen_SKCM$CN)]=2
table(is.na(Neoantigen_SKCM$CN))
table(is.na(Neoantigen_SKCM$purity))
table(is.na(Neoantigen_SKCM$VAF))
Neoantigen_SKCM$CCF <- Neoantigen_SKCM$VAF*Neoantigen_SKCM$CN/Neoantigen_SKCM$purity
summary(Neoantigen_SKCM$CCF)#NA: 10
VAF <- Neoantigen_SKCM$VAF
CN <- Neoantigen_SKCM$CN
purity <- Neoantigen_SKCM$purity
m=VAF/purity*(purity*CN+2*(1-purity))#Multiplicity of a mutation: the number of DNA copies bearing a mutation m
summary(m)
summary(round(m))
hist(m, col = rgb(1,0,0,0.2),freq = F)
lines(density(m), col = "red")
hist(round(m), col = rgb(1,0,0,0.2),freq = F)
lines(density(m), col = "red")
m <- ifelse(round(m)>0,round(m),1)
summary(m)
CCF=VAF/(m*purity)*(purity*CN+2*(1-purity))
summary(CCF)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.048   0.596   0.860   0.822   1.055   1.500      10
num1 <- which(substr(Neoantigen_SKCM$Barcode,14,15)=="01")
num2 <- which(substr(Neoantigen_SKCM$Barcode,14,15)=="06")
site <- rep("Metastatic",nrow(Neoantigen_SKCM))
site[num1] <- "Primary"
site[num2] <- "Metastatic"
summary(Neoantigen_SKCM$CCF)
table(is.na(Neoantigen_SKCM$CCF))
Neoantigen_SKCM$CCF <- CCF
table(is.na(Neoantigen_SKCM$CCF))#8414-->0-->20
#-->45-->10
chimer_neofus$CCF <- ifelse(chimer_neofus$genetype=="   ",0.5,1)
Neoantigen_SKCM$CCF[-grep(".p",Neoantigen_SKCM$Mutation)]<- chimer_neofus$CCF[match(Neoantigen_SKCM$Mutation[-grep(".p",Neoantigen_SKCM$Mutation)], chimer_neofus$Fusion_pair)]
summary(Neoantigen_SKCM$CCF)
table(Neoantigen_SKCM$Mutation %in% chimer_neofus$Fusion_pair)
table(Neoantigen_SKCM$Mutation[-grep(".p",Neoantigen_SKCM$Mutation)] 
      %in% chimer_neofus$Fusion_pair)

Neoantigen_SKCM$CCF <- ifelse(Neoantigen_SKCM$CCF>=1,1,Neoantigen_SKCM$CCF)
summary(Neoantigen_SKCM$CCF)
Neoantigen_SKCM$clone <- ifelse(Neoantigen_SKCM$CCF>=0.84,"Clone","Subclone")
##----
Neoantigen_SKCM.ITH <- Neoantigen_SKCM#[grep(".p",Neoantigen_SKCM$Mutation),]
summary(Neoantigen_SKCM.ITH$CCF)
table(is.na(Neoantigen_SKCM.ITH$purity))#2127
table(is.na(Neoantigen_SKCM.ITH$CCF))#2018
Neoantigen_SKCM.ITH <- Neoantigen_SKCM.ITH[!is.na(Neoantigen_SKCM.ITH$CCF),]
hist(Neoantigen_SKCM.ITH$CCF, col = rgb(1,0,0,0.2),freq = F)
lines(density(Neoantigen_SKCM.ITH$CCF), col = "red")
table(Neoantigen_SKCM.ITH$clone)
Neoantigen_SKCM.ITH$Neoantigen.model <- Neoantigen.data_SKCM$Neoantigen.model[match(Neoantigen_SKCM.ITH$Barcode,rownames(Neoantigen.data_SKCM))
]
library(ggplot2)
library(scales)
library(reshape2)
ITH <- dcast(Neoantigen_SKCM.ITH,Barcode~clone, length)
rownames(ITH) <- ITH$Barcode
dim(ITH)#468   3
ITH$Samples <- substr(ITH$Barcode,1,15)
setdiff(Neoantigen.data_SKCM[!is.na(Neoantigen.data_SKCM$Purity),]$Samples,ITH$Samples)
#"TCGA-EB-A5VV-06"
#ITH <- ITH[ITH$NA==0,]
table(is.na(Neoantigen.data_SKCM$Purity))
table(is.na(match(Neoantigen.data_SKCM$Samples, substr(ITH$Barcode,1,15))))
Neoantigen.data_SKCM$Clone <- ITH$Clone[match(Neoantigen.data_SKCM$Samples,
                                              substr(ITH$Barcode,1,15))]
Neoantigen.data_SKCM$Subclone <- ITH$Subclone[match(Neoantigen.data_SKCM$Samples,
                                                    substr(ITH$Barcode,1,15))]
########---merge phenotypes & NEO2IS & NEOITH score---#######
colnames(clin.all)
table(Neoantigen.data_SKCM$Samples %in% clin.all$Tumor_Sample_Barcode)
dataAB.SKCM<-clin.all[match(Neoantigen.data_SKCM$Samples,clin.all$Tumor_Sample_Barcode),-(17:18)]
table(dataAB.SKCM$Stage)
table(dataAB.SKCM$Outcome)
table(dataAB.SKCM$Age)
dim(dataAB.SKCM)##469   17
table(is.na(dataAB.SKCM$Age))#28
num1 <- which(substr(dataAB.SKCM$Tumor_Sample_Barcode,14,15)=="01")
num1
num2 <- which(substr(dataAB.SKCM$Tumor_Sample_Barcode,14,15)=="06")
num2
dataAB.SKCM$site <- rep("Metastatic",nrow(dataAB.SKCM))
dataAB.SKCM$site[num1] <- "Primary"
dataAB.SKCM$site[num2] <- "Metastatic"
table(dataAB.SKCM$site)

dataAB.SKCM<-cbind(dataAB.SKCM,Neoantigen.data_SKCM)
table(is.na(dataAB.SKCM$Purity))#102
table(dataAB.SKCM$Samples %in% data_maf.SKCM$Tumor_Sample_Barcode)
table(is.na(dataAB.SKCM$Clone))
dataAB.SKCM$ITH.score <- dataAB.SKCM$Subclone/(dataAB.SKCM$Clone+dataAB.SKCM$Subclone)
summary(dataAB.SKCM$ITH.score)#NA:1
dataAB.SKCM$ITH.score[is.na(dataAB.SKCM$ITH.score)] <- 1
summary(dataAB.SKCM$ITH.score)#0
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.339   0.437   0.485   0.611   1.000 (CCF>=0.84)

#########---merge TCR diversity---#########
ABSOLUTE.TCGA <- read.table("../TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                            sep="\t",fill = T,header = T)
ABSOLUTE.SKCM <- ABSOLUTE.TCGA[ABSOLUTE.TCGA$array %in% substr(dataAB$Samples,1,15),]
TCGA.info <- read.table('/public/home/lorihan/lrh/Neoantigen_TCGA/mitcr_sampleStatistics_20160714.tsv',sep = "\t",header = T)
library(openxlsx)
TCGA.ITH.all <- read.xlsx('/public/home/lorihan/lrh/Neoantigen_TCGA/TCGA-ITH-info.xlsx')#,1)
table(is.na(match(substr(dataAB.SKCM$Samples,1,12),
                  TCGA.ITH.all$TCGA.Participant.Barcode)))
TCGA.ITH.SKCM <- TCGA.ITH.all[ TCGA.ITH.all$TCGA.Participant.Barcode %in% substr(dataAB.SKCM$Samples,1,12),]
colnames(TCGA.ITH.SKCM)[1] <- colnames(dataAB.SKCM)[1]
TCGA.ITH.SKCM$Tumor_Sample_Barcode <- dataAB.SKCM$Samples[match(TCGA.ITH.SKCM$Tumor_Sample_Barcode,substr(dataAB.SKCM$Samples,1,12))]
dataAB.SKCM <- dplyr::left_join(dataAB.SKCM,TCGA.ITH.SKCM[,c(1,26:28)],by=colnames(dataAB.SKCM)[1])
table(is.na(match(dataAB.SKCM$Samples,
                  substr(TCGA.info$AliquotBarcode,1,15))))
dataAB.SKCM$Shanno <- TCGA.info$shannon[match(dataAB.SKCM$Samples,
                                              substr(TCGA.info$AliquotBarcode,1,15))]
dataAB.SKCM$Purity<- ABSOLUTE.SKCM$purity[match(dataAB.SKCM$Samples,ABSOLUTE.SKCM$array)]
dataAB.SKCM$Ploidy<- ABSOLUTE.SKCM$ploidy[match(dataAB.SKCM$Samples,ABSOLUTE.SKCM$array)]
Purity.na.samples <- dataAB.SKCM$Samples[is.na(dataAB.SKCM$Purity)]
dataAB.SKCM$Purity[is.na(dataAB.SKCM$Purity)] <- tumorPurity$tumorPurity[match(Purity.na.samples,substr(rownames(tumorPurity),1,15))]
dataAB.SKCM$Ploidy[is.na(dataAB.SKCM$Ploidy)] <- tumorPurity$tumorPloidy[match(Purity.na.samples,substr(rownames(tumorPurity),1,15))]

#########---cor scRNA-deconvolution data---#########
load("/public/home/lorihan/lrh/SKCM/ciber_perm1000_skcm.gsea.CD8T.ctrl.Rdata")
colnames(TME.results)
rownames(TME.results) <- substr(rownames(TME.results),1,15)
#colnames(TME.results) <- gsub("CD8Tex","CD8_Tex",colnames(TME.results))
TME.results <- TME.results[rownames(TME.results) %in% Neoantigen.data_SKCM$Samples,]
library(dplyr)
library(tidyr)
colnames(TME.results)
CD8__T_cell <- apply(TME.results[,c(1,3)],1,function(x){
  sum(x)
})
TME.results[,1] <- CD8__T_cell
dd1 <- TME.results[,-3] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:(ncol(TME.results)-1-2),
               names_to = "CellType",
               values_to = "Composition")
colnames(dd1)
plot.info <- dd1[,c(5,1,6)]  

plot.info$CellType <- gsub(" Exhausted CD8"," exhausted CD8+",plot.info$CellType)
plot.info$CellType <- gsub("cells","cell",plot.info$CellType)
plot.info$CellType <- gsub("Terminal","Terminally",plot.info$CellType)
library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)]
vc_cols =RColorBrewer::brewer.pal(n = 12, name = 'Paired')#mycol
pal <- vc_cols[c(1,3,5,6,8,9)]
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot","ggpubr","bslib","ggthemes")
lapply(pkgs, library, character.only = T)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",palette = pal,
  xlab = "",
  ylab = ""# "CIBERSORTx-inferred cell fractions",
) + #ggtitle("CD8+ T cell subsets (TCGA-SKCM)")+
  theme_base() +
  theme(legend.position='none',
        plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
        plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
        plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
        axis.text.x   = element_text(color = 'black', size = 14, angle = 60,
                                     hjust = 1,vjust = 1),
        axis.text.y   = element_text(color = 'black', size = 14, angle = 0),
        axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
        axis.title.y  = element_text(color = 'black', size = 13.5, angle = 90),
        axis.line.y = element_line(color = 'black', linetype = 'solid'), # y
        axis.line.x = element_line (color = 'black',linetype = 'solid'), # x
  )
ggsave("SKCM468.CD8Tex.pdf",width = 7.5,height = 6)

TME.results.SKCM <- TME.results
library(tidyr)
library(tibble)
dat.CD8Exh <- TME.results %>%
  as.data.frame() %>%
  rownames_to_column("Samples")
colnames(dat.CD8Exh)
dat.CD8Exh$Samples <- gsub("[.]","-",dat.CD8Exh$Samples)
dat.CD8Exh <- dat.CD8Exh[dat.CD8Exh$Samples %in% dataCD.SKCM$Tumor_Sample_Barcode,]
colnames(dat.CD8Exh) <- gsub(' ','_',colnames(dat.CD8Exh))
colnames(dat.CD8Exh) <- gsub('[+]','_',colnames(dat.CD8Exh))
colnames(dat.CD8Exh) <- gsub("[)]","",gsub("[(]","",colnames(dat.CD8Exh)))
table(dataCD.SKCM$Samples %in% dat.CD8Exh$Samples)
descrCorr = dplyr::left_join(dataCD.SKCM, dat.CD8Exh[,1:(ncol(dat.CD8Exh)-3)],by='Samples')
descrCorr.var <- descrCorr[,c(colnames(dat.CD8Exh)[2:(ncol(dat.CD8Exh)-3)],'ITH.score','Neoantigen.model',
                              'CD8.effector',"CD8Tex",'Immunoediting.score','TCR.Shannon','TCR.Richness','TCR.Evenness')]
library(Hmisc)
res2<-rcorr(as.matrix(descrCorr.var))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
cor.res <- flattenCorrMatrix(res2$r,res2$P)

colnames(descrCorr)
cor.test(descrCorr$TCR.Shannon,descrCorr$Naive_CD8__T_cell)
cor.test(descrCorr$TCR.Shannon,descrCorr$Effector_CD8__memory_T_Tem_cell)
cor.test(descrCorr$ITH.score,descrCorr$Effector_CD8__memory_T_Tem_cell)
cor.test(descrCorr$ITH.score,descrCorr$Naive_CD8__T_cell)
descrCorr.SKCM <- descrCorr
descrCorr.var$site <- dataAB.SKCM$site[match(descrCorr$Samples,dataAB.SKCM$Samples)]
ggscatter(descrCorr.SKCM, x = "Terminal_Exhausted_CD8_T_cells", y = "TCR.Shannon",
          color = "#E5C494", size = 1, #color = "site", palette = "nature",#
          add = "reg.line", add.params = list(color = "#B2DF8A", fill = "lightgray"), # 
          conf.int = TRUE, # 
          cor.coef = TRUE, # 
          cor.coeff.args = list(method = "pearson", label.x = 0.45,label.y=5.6, label.sep =", ") #"\n")#
)+ ylab("TCR diversity") + xlab("Terminally exhausted CD8+ T cells")#+  stat_cor(aes(color = site), label.x = 3)
ggsave(file="/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/SKCM_TCR_terminalCD8Tex.pdf",width = 4,height = 4)
ggscatter(descrCorr.SKCM, x = "Effector_CD8__memory_T_Tem_cell", y = "TCR.Shannon",
          color = "#A6CEE3", size = 1, #color = "site", palette = "nature",
          add = "reg.line", add.params = list(color = "salmon", fill = "lightgray"), # 
          conf.int = TRUE, # 
          cor.coef = TRUE, # 
          cor.coeff.args = list(method = "pearson", label.x = 0.6,label.y=2.2, label.sep =", ") #"\n")#
)+ ylab("TCR diversity") + xlab("Effector/memory CD8+ T cell")#+  stat_cor(aes(color = site), label.x = 3)
ggsave(file="/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/SKCM_TCR_effectorCD8.pdf",width = 4,height = 4)
