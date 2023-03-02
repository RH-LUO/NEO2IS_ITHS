# title: "SKCM PhenoGenotypes"
# author: "Ruihan Luo"
# date: "Jul 5th,2022"
setwd("/public/home/lorihan/lrh/SKCM")
.libPaths()  
myPaths <- .libPaths()  
new <- '/public/home/lorihan/miniconda3/lib/R/library'
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
.libPaths()
###########
TCGA_SKCM_VCF_Info <- read.table("/public/home/lorihan/lrh/SKCM/gdc_sample_sheet.2022-03-23.tsv",sep = "\t",fill = T,header = T)
table(TCGA_SKCM_VCF_Info$Data.Type)
TCGA_SKCM_VCF_Info.raw <- TCGA_SKCM_VCF_Info[TCGA_SKCM_VCF_Info$Data.Type=="Raw Simple Somatic Mutation",]
unique(TCGA_SKCM_VCF_Info.raw$Sample.ID)#470 caseID, 472 sampleID
###########
getIDtable <- function(metaFile){
  suppressPackageStartupMessages(library("jsonlite"))
  meta = jsonlite::fromJSON(metaFile)
  IDtable = meta[, c("submitter_id","file_size", "file_name", "file_id","associated_entities")]
  TCGA_barcode = sapply(IDtable$associated_entities, function(x) paste(x[[1]],collapse = ","))
  IDtable[, 5] = TCGA_barcode
  colnames(IDtable)[5] = "barcode"
  return(IDtable)
}
TCGA_SKCM_VCF_meta <- getIDtable("/public/home/lorihan/lrh/SKCM/metadata.cart.2022-03-23.json")
unique(TCGA_SKCM_VCF_meta$barcode)#472
TCGA_SKCM_VCF <- TCGA_SKCM_VCF_meta#TCGA_SKCM_VCF_files
#TCGA_SKCM_VCF$Barcode <- TCGA_SKCM_VCF_meta$submitter_id[match(TCGA_SKCM_VCF$filename,TCGA_SKCM_VCF_meta$file_name)]
TCGA_SKCM_VCF$Barcode <- str_split(TCGA_SKCM_VCF$submitter_id,"[_]",simplify = T)[,1]
TCGA_SKCM_VCF$Matched <- str_split(TCGA_SKCM_VCF$submitter_id,"[_]",simplify = T)[,2]
TCGA_SKCM_VCF$Sample <- TCGA_SKCM_VCF_Info$Sample.ID[match(TCGA_SKCM_VCF$file_name,TCGA_SKCM_VCF_Info$File.Name)]
TCGA_SKCM_VCF$Patients <- substr(TCGA_SKCM_VCF$Sample,1,12)
TCGA_SKCM_VCF$Workflow <-  str_split(TCGA_SKCM_VCF$submitter_id,"[_]",simplify = T)[,3]
unique(TCGA_SKCM_VCF$Barcode)#472
table(substr(TCGA_SKCM_VCF$Barcode,14,16))
# 01A  06A  06B 
#832 2928   16 
table(TCGA_SKCM_VCF$Workflow)
#table(TCGA_SKCM_VCF$file_name %in% TCGA_SKCM$filename)#618 mutect2
length(TCGA_SKCM_VCF$Barcode[grep("-01A",TCGA_SKCM_VCF$Barcode)])
length(TCGA_SKCM_VCF$Barcode[grep("-06A",TCGA_SKCM_VCF$Barcode)])
unique(TCGA_SKCM_VCF$barcode)#472
#unique(str_split(TCGA_SKCM_VCF$barcode,"[,]",simplify = T)[,1])#575
#unique(str_split(TCGA_SKCM_VCF$barcode,"[,]",simplify = T)[,2])#572

#filter for RawVcfs
TCGA_SKCM_VCF<- TCGA_SKCM_VCF[TCGA_SKCM_VCF$file_name %in% TCGA_SKCM_VCF_Info.raw$File.Name,]
unique(TCGA_SKCM_VCF$Barcode)#472
unique(TCGA_SKCM_VCF$Sample)#472
#filter for Workflow
TCGA_SKCM_VCF<- TCGA_SKCM_VCF[TCGA_SKCM_VCF$Workflow=="mutect",]#113
# TCGA_SKCM_VCF<- TCGA_SKCM_VCF.all[TCGA_SKCM_VCF.all$Workflow=="varscan",]#107
# TCGA_SKCM_VCF<- TCGA_SKCM_VCF.all[TCGA_SKCM_VCF.all$Workflow=="somaticsniper",]#116
# TCGA_SKCM_VCF<- TCGA_SKCM_VCF.all[TCGA_SKCM_VCF.all$Workflow=="muse",]#116

####
length(unique(substr(TCGA_SKCM_VCF$Barcode,1,16)))#472
length(unique(TCGA_SKCM_VCF$Sample))#472
length(grep("-06B",TCGA_SKCM_VCF$Barcode))
unique(TCGA_SKCM_VCF$Sample[-grep("-06B",TCGA_SKCM_VCF$Barcode)])
#470
unique(TCGA_SKCM_VCF$Patients[-grep("-06B",TCGA_SKCM_VCF$Barcode)])#468
unique(TCGA_SKCM_VCF$Patients[-grep("-06B",TCGA_SKCM_VCF$Sample)])#468
table(substr(TCGA_SKCM_VCF$Barcode,14,16))
#01A 06A 06B 
#104 366   2 
grep("-06B",TCGA_SKCM_VCF$Sample)
TCGA_SKCM_VCF$Sample <- substr(TCGA_SKCM_VCF$Barcode,1,16)
length(unique(TCGA_SKCM_VCF$Sample))#472
length(unique(TCGA_SKCM_VCF$Patients))#470
TCGA_SKCM_VCF[duplicated(TCGA_SKCM_VCF$Patients),]
#TCGA-ER-A19T-01/6A TCGA-ER-A2NF-01/6A
TCGA_SKCM_VCF <- TCGA_SKCM_VCF[!duplicated(TCGA_SKCM_VCF$Patients),]
dim(TCGA_SKCM_VCF)
####
#filter for patients
table(substr(TCGA_SKCM_VCF$Barcode,14,16))
#01A 06A 06B 
#104 366-2   2 
dim(TCGA_SKCM_VCF)#[1] 472-2  10
unique(TCGA_SKCM_VCF$Patients)#470
length(grep("-01A-",TCGA_SKCM_VCF$Barcode))#104
# tmp=by(TCGA_SKCM_VCF,TCGA_SKCM_VCF$Patients,function(x) rownames(x)[grepl("01A",x$Barcode)])
# length(unique(TCGA_SKCM_VCF$Sample))#535
# probes = as.character(tmp)
# #TCGA_SKCM_VCF=TCGA_SKCM_VCF[rownames(TCGA_SKCM_VCF) %in% probes,]
unique(TCGA_SKCM_VCF$Barcode)#472-2
unique(TCGA_SKCM_VCF$Patients)#470
######
SKCM.clin <- read.csv("/public/home/lorihan/lrh/SKCM/TCGA-SKCM.GDC_phenotype.tsv.gz",sep = '\t',header = TRUE)
SKCM.survival<-read.csv('survival%2FSKCM_survival.txt',sep = '\t',header = T)
SKCM.pheno<-read.table("TCGA.SKCM.sampleMap%2FSKCM_clinicalMatrix",sep = '\t',header = TRUE)
table(unique(substr(TCGA_SKCM_VCF$Sample,1,15)) %in% SKCM.pheno$sampleID)
table(unique(substr(TCGA_SKCM_VCF$Sample,1,15)) %in% SKCM.survival$sample)
table(unique(TCGA_SKCM_VCF$Sample) %in% SKCM.clin$submitter_id.samples)
#SKCM.pheno[,c("days_to_last_followup","RFS.time","OS.time")]

colnames(SKCM.pheno)
clinical_var<-c('sampleID','CDE_ID_3226963','CIMP','MSI_updated_Oct62011','X_PANCAN_DNAMethyl_PANCAN','braf_gene_analysis_result','kras_mutation_found',
                'gender','days_to_birth','additional_pharmaceutical_therapy',
                'additional_radiation_therapy','history_of_colon_polyps',
                'hypermutation','days_to_initial_pathologic_diagnosis','days_to_last_followup',
                'lymphatic_invasion','venous_invasion','pathologic_M','pathologic_N','pathologic_T','pathologic_stage',
                'BRAF','CTNNB1','Canonical_mut_in_KRAS_EGFR_ALK','EGFR','ERBB2','ERBB4','KRAS','MET','PIK3CA','STK11',
                'X_primary_disease','tobacco_smoking_history_indicator','days_to_death','number_pack_years_smoked','years_smoked.exposures',
                'OS_STATUS','DFS_STATUS','lymphatic_metastasis','distant_metastasis','RFS.year','OS.year')

clinical_var
paste(clinical_var,collapse="','")
colnames(SKCM.pheno)
length(clinical_var)
table(colnames(SKCM.pheno) %in% clinical_var)
clinical_var[!(clinical_var %in% colnames(SKCM.pheno))]

clinical<-SKCM.pheno[,colnames(SKCM.pheno) %in% clinical_var]
table(SKCM.survival$sample %in% clinical$sampleID)
table(SKCM.survival$OS)
head(SKCM.survival)
survival<-SKCM.survival[,c('sample','OS','PFI','OS.time','PFI.time')]
colnames(survival)<-c('sampleID','OS_STATUS','DFS_STATUS','OS.time','DFS.time')
table(is.na(match(clinical$sampleID,SKCM.survival$sample)))
table(unique(substr(TCGA_SKCM_VCF$Sample,1,15)) %in% clinical$sampleID)

clinical<-merge.data.frame(clinical,survival,by = "sampleID",all.x = TRUE)
table(is.na(clinical$OS_STATUS))
table(is.na(clinical$pathologic_N))
str(clinical)
table(clinical$pathologic_stage)
###clear up clinical info###
stage<-as.character(clinical$pathologic_stage)
table(stage)
table(is.na(stage))
names(stage)<-clinical$sampleID
stage[grep("Stage IV",stage)]="III/IV" 
stage[grep("Stage III",stage)]="III/IV" 
stage[grep("Stage II",stage)]="I/II" 
stage[grep("Stage I",stage)]="I/II" 
stage[grep("I/II NOS",stage)]="I/II" 
stage[grep("Stage 0",stage)]="0" 
#stage[grep("[Discrepancy]",stage)]="Unknown" 
stage[stage==""]="Unknown" 
levels(factor(stage))
table(stage)

tumor_size<-as.character(clinical$pathologic_T)
table(tumor_size)
table(is.na(tumor_size))
names(tumor_size)<-clinical$sampleID
tumor_size[grep("T4",tumor_size)]="T3/T4" 
tumor_size[grep("T3",tumor_size)]="T3/T4" 
tumor_size[grep("T2",tumor_size)]="T1/T2" 
tumor_size[grep("T1",tumor_size)]="T1/T2" 
#tumor_size[grep("Tis",tumor_size)]="T1/T2" 
tumor_size[grep("TX",tumor_size)]<-""
tumor_size[tumor_size==""]="Unknown" 
levels(factor(tumor_size))
table(tumor_size)

lymphatic_metastasis<-as.character(clinical$pathologic_N)
table(lymphatic_metastasis)
table(is.na(lymphatic_metastasis))
names(lymphatic_metastasis)<-clinical$sampleID
lymphatic_metastasis[grep("N1",lymphatic_metastasis)]="Present" 
lymphatic_metastasis[grep("N2",lymphatic_metastasis)]="Present"
lymphatic_metastasis[grep("N3",lymphatic_metastasis)]="Present" 

lymphatic_metastasis[grep("N0",lymphatic_metastasis)]="Absent" 
lymphatic_metastasis[grep("NX",lymphatic_metastasis)]="Unknown" 
lymphatic_metastasis[lymphatic_metastasis==""]="Unknown" 
levels(factor(lymphatic_metastasis))
table(lymphatic_metastasis)
##
Relapse<-as.character(clinical$DFS_STATUS)
names(Relapse)<-clinical$sampleID
table(Relapse)
table(is.na(Relapse))
Relapse[grep("0",Relapse)]="Relapse-free" 
Relapse[grep("1",Relapse)]="Relapse" 
Relapse[is.na(Relapse)]="Unknown" 
table(Relapse)
Recurrence<-as.character(SKCM.pheno$new_tumor_event_after_initial_treatment)
names(Recurrence)<-SKCM.pheno$sampleID
table(Recurrence)
table(is.na(Recurrence))
Recurrence[grep("NO",Recurrence)]="Recurrence-free" 
Recurrence[grep("YES",Recurrence)]="Recurrence" 
Recurrence[Recurrence==""]="Unknown" 
Recurrence[Recurrence=="[Discrepancy]"]="Unknown" 
table(Recurrence)

##
distant_metastasis<-as.character(clinical$pathologic_M)
table(distant_metastasis)
table(is.na(distant_metastasis))
names(distant_metastasis)<-clinical$sampleID
distant_metastasis[grep("M1",distant_metastasis)]="Present" 
distant_metastasis[grep("M0",distant_metastasis)]="Absent" 
distant_metastasis[grep("MX",distant_metastasis)]="Unknown" 
distant_metastasis[distant_metastasis==""]="Unknown" 
levels(factor(distant_metastasis))
table(distant_metastasis)
##
Age<- as.numeric(-clinical$days_to_birth/365)
table(is.na(Age))
table(Age=="")
Age<-ifelse(!is.na(Age),ifelse(Age>=65,">=65","<65"),"Unknown")
table(Age)

Sex <- as.character(clinical$gender)
table(is.na(Sex))
table(Sex=="")
table(Sex)
Sex<-ifelse(Sex!="",
            ifelse(Sex=="FEMALE","Female","Male"),"Unknown")
table(Sex)
table(is.na(clinical$OS_STATUS))
table(clinical$OS_STATUS)
table(is.na(clinical$DFS_STATUS))
table(clinical$DFS_STATUS)
##
# response<-as.character(SKCM.pheno$primary_therapy_outcome_success)
# table(response)
# table(is.na(response))
# grep("Discrepancy",response)
# response[grep("Discrepancy",response)]<-""
# table(response)
# levels(factor(response))
####
colnames(clinical)
paste(colnames(clinical),collapse = "','")
Dat<-clinical[,c('sampleID','X_primary_disease','days_to_birth','gender','pathologic_M','pathologic_N','pathologic_T',
                 'pathologic_stage','OS_STATUS','DFS_STATUS','OS.time','DFS.time')]
colnames(Dat)<-c('sampleID','Histology','Age','Sex','M.stage','N.stage','T.stage',
                 'Stage','OS_STATUS','DFS_STATUS','OS.month','DFS.month')
table(Dat$Histology)
table(is.na(Dat$Histology))
Dat$Age<-Age;Dat$Sex<-Sex
Dat$M.stage<-distant_metastasis;Dat$N.stage<-lymphatic_metastasis;
Dat$T.stage<-tumor_size;Dat$Stage<-stage; 
table(Dat$Pack_years_smoked)
table(is.na(Dat$Pack_years_smoked))
#Dat$Response<-response
colnames(Dat)
Dat$OS.month<-Dat$OS.month/30; Dat$DFS.month<-Dat$DFS.month/30
table(Dat$Response)
table(Dat$DFS_STATUS)
Dat <- within(Dat, {
#  Response<-factor(Response ,labels = c("Unknown","Complete Remission","Partial Remission","Progressive Disease","Stable Disease"))
  OS_STATUS <- factor(OS_STATUS , labels = c("Alive","Dead"))
  DFS_STATUS <- factor(DFS_STATUS , labels = c("Relapse-free","Relapse"))
})
str(Dat)#481  12
table(unique(substr(TCGA_SKCM_VCF$Sample,1,15)) %in% Dat$sampleID)

######
load("GDC_TCGA_SKCM_clinical_df.Rdata")
clin_SKCM<-clinical.patient.SKCM
unique(clin_SKCM$bcr_patient_barcode)
table(unique(substr(TCGA_SKCM_VCF$Sample,1,12)) %in% clin_SKCM$bcr_patient_barcode)
table(clinical.drug.SKCM$drug_name)
str(clinical.drug.SKCM)
class(clinical.drug.SKCM)
clinical.drug.SKCM[,-1]<-lapply(clinical.drug.SKCM[,-1], as.character)
clinical.drug<- clinical.drug.SKCM
str(clinical.drug)
table(clinical.drug$therapy_types)
clinical.drug<-clinical.drug[match(unique(clinical.drug$bcr_drug_barcode),clinical.drug$bcr_drug_barcode),]
#clinical.drug <- clinical.drug[!duplicated(clinical.drug$bcr_patient_barcode),]
table(clinical.drug$therapy_types)
drug.table <- table(clinical.drug$drug_name)
drug_table <- table(clinical.drug.SKCM$drug_name)

######
unique(clinical.drug.SKCM$bcr_drug_barcode)#265
patients1<-unique(substr(Dat$sampleID,1,12))#471
patients2<-unique(substr(clin_SKCM$bcr_patient_barcode,1,12))#470

patient.drug<-unique(clinical.drug$bcr_patient_barcode)#146
patient.diff<-setdiff(patients1,patient.drug)#325
clin.fake <- rbind(clinical.drug.SKCM,clinical.drug)
clin.fake <- head(clin.fake,length(patient.diff))
clin.fake$bcr_patient_barcode<-patient.diff
clin.fake[,-1]<-"Unknown"

clin.fake<-rbind(clinical.drug,clin.fake)
match(patients1,patients2)
table(patients1 %in% patients2)
table(patients1 %in% patient.drug)##TRUE:146
Dat$Tumor_Sample_Barcode<-substr(Dat$sampleID,1,12)
table(is.na(match(Dat$Tumor_Sample_Barcode,patient.drug)))##FALSE:153

clin_meta<-clinical.drug[clinical.drug$bcr_patient_barcode %in% Dat$Tumor_Sample_Barcode,]
clin_meta1<-Dat[match(clin_meta$bcr_patient_barcode,Dat$Tumor_Sample_Barcode),]
clin_meta2<-cbind(clin_meta,clin_meta1)
clin_meta2$sampleID[!clin_meta2$sampleID  %in% unique(substr(TCGA_SKCM_VCF$Sample,1,15))]

clin_meta2$Tumor_Sample_Barcode<-substr(clin_meta2$sampleID,1,15)
remove(clin_meta1)
clin_meta<-clin_meta2
length(unique(clin_meta$bcr_patient_barcode))#323 --> 146
table(clin_meta$Stage)
#0    I/II  III/IV Unknown 
#4     130     110      21 
#str(clin_meta)
table(clin_meta$measure_of_response)
clin_meta <- within(clin_meta, {
  measure_of_response<-factor(measure_of_response,labels = c("Unknown","Progressive Disease","Complete Response","Partial Response","Stable Disease"))
})
table(clin_meta$measure_of_response)
table(clin_meta$Response)
unique(clin_meta$Tumor_Sample_Barcode)##146
unique(clin_meta$bcr_patient_barcode)##146
######
na.omit(match(Dat$Tumor_Sample_Barcode,clin.fake$bcr_patient_barcode))
unique(Dat$Tumor_Sample_Barcode)
##471
clin_res<-clin.fake[clin.fake$bcr_patient_barcode %in% Dat$Tumor_Sample_Barcode,]
clin_res1<-Dat[match(clin_res$bcr_patient_barcode,Dat$Tumor_Sample_Barcode),]
clin_res2<-cbind(clin_res,clin_res1)
table(unique(substr(TCGA_SKCM_VCF$Sample,1,15)) %in% clin_res1$sampleID)
table(unique(substr(TCGA_SKCM_VCF$Sample,1,12)) %in% clin_res$bcr_patient_barcode)

clin_res2$Tumor_Sample_Barcode<-substr(clin_res1$sampleID,1,15)
remove(clin_res1)
clin_res<-clin_res2
remove(clin_res2)
table(clin_res$Stage)
#0    I/II  III/IV Unknown 
#10     287     245      48 
length(unique(clin_res$bcr_patient_barcode))#471
table(Dat$M.stage,Dat$Stage)
length(unique(clin_res$bcr_patient_barcode))#471
table(clin_res$Stage)
# 0    I/II  III/IV Unknown 
#10     287     245      48  
str(clin_res)
table(clin_res$measure_of_response)
clin_res$measure_of_response<-gsub("Unknown","",clin_res$measure_of_response)
table(clin_res$measure_of_response)
clin_res <- within(clin_res, {
   measure_of_response<-factor(measure_of_response,labels = c("Unknown","Progressive Disease","Complete Response","Partial Response","Stable Disease"))
 })
 table(clin_res$measure_of_response)
 table(clin_res$Response)
unique(clin_res$Tumor_Sample_Barcode)##471
unique(clin_res$bcr_patient_barcode)##471
#write.table(clin_res,"TCGA_LUNG.clin.txt",sep="\t",quote=FALSE,row.names = F)
#write.table(clin_meta,"TCGA_LUNG.clin_drug.meta.txt",sep="\t",quote=FALSE,row.names = F)

##
colnames(Dat)
clinAB<-clin_meta
colnames(clinAB)
clinAB<-clinAB[,c("Tumor_Sample_Barcode","DFS_STATUS","DFS.month","OS_STATUS","OS.month",
                  "Histology","Sex","Age","T.stage","N.stage",
                  "M.stage","Stage","days_to_drug_therapy_start","days_to_drug_therapy_end",
                  "bcr_patient_barcode","therapy_types","drug_name","measure_of_response")]
table(clinAB$therapy_types)
#clinAB<-clinAB[clinAB$therapy_types!="Ancillary",]
#clinAB<-clinAB[clinAB$therapy_types!="Vaccine",]
#clinAB<-clinAB[-grep("Other,",clinAB$therapy_types),]
table(clinAB$therapy_types)
table(clinAB$drug_name)
drug<-clinAB$drug_name
drug<-tolower(drug)
table(drug)
library(stringr)
drug<-sapply(str_split(drug,"[, ]"),'[',1)
drug<-sapply(str_split(drug,"[/]"),'[',1)
table(drug)
drug[drug=="dacabarzine"]="dacarbazine"
drug[drug==""]="unknown"
# drug[drug=="erlotonib"]="erlotinib"
# drug[drug=="vinorelbin"]="vinorelbine"
# drug[drug=="carboplatinum"]="carboplatin"
# drug[drug=="paraplatin"]="carboplatin"
# drug[drug=="premetrexed"]="pemetrexed"
# drug[drug=="cisplatinum"]="cisplatin"
# drug[drug=="ironotecan"]="irinotecan"
# drug[drug=="cpt-11"]="irinotecan"
# drug[drug=="docetoxel"]="docetaxel"
# drug[drug=="doxetaxol"]="docetaxel"
drug[drug=="actinomycin"]="actinomycin-d"
drug[drug=="recmage-"]="mage-3"
drug[drug=="mage"]="mage-3"

drug[drug=="temozolomide"]="temodal"
drug[drug=="temodar"]="temodal"
drug[drug=="vincristine"]="vinblastine"
drug[drug=="yervoy"]="ipilimumab"

table(drug)
clinAB$drug_name<-drug
table(clinAB$drug_name)
table(clinAB$therapy_types)

therapy_types<-clinAB$therapy_types
table(therapy_types)
grep("nib",clinAB$drug_name)
grep("mab",clinAB$drug_name)
drug[grep("mab",clinAB$drug_name)]
drug[grep("nib",clinAB$drug_name)]

#therapy_types[grep("alimta",clinAB$drug_name)]="Chemotherapy"
therapy_types[grep("nib",clinAB$drug_name)]="Targeted Molecular therapy"
therapy_types[grep("mab",clinAB$drug_name)]="Immunotherapy"
therapy_types[grep("avastin",clinAB$drug_name)]
therapy_types[grep("tarceva",clinAB$drug_name)]
# therapy_types[grep("iressa",clinAB$drug_name)]
# therapy_types[grep("ras",clinAB$drug_name)]

therapy_types[grep("avastin",clinAB$drug_name)]="Targeted Molecular therapy"
#therapy_types[grep("tarceva",clinAB$drug_name)]="Targeted Molecular therapy"
#therapy_types[grep("iressa",clinAB$drug_name)]="Targeted Molecular therapy"
#therapy_types[grep("ras",clinAB$drug_name)]="Targeted Molecular therapy"

drug[therapy_types=="Targeted Molecular therapy"]
drug[grep("Targeted Molecular therapy",therapy_types)]
table(therapy_types)

clinAB$drug_name[clinAB$therapy_types==""]

therapy_types[clinAB$drug_name=="alimta"]="Chemotherapy"
therapy_types[clinAB$drug_name=="carboplatin"]="Chemotherapy"
table(therapy_types)
clinAB[therapy_types=="",]
clinAB$therapy_types<-therapy_types
table(clinAB$therapy_types)

clinCD<-clinAB[clinAB$therapy_types=="Immunotherapy",]
length(unique(clinCD$Tumor_Sample_Barcode))#76
clin.immuno <- clinAB[clinAB$drug_name %in% c("ipilimumab","pembrolizumab","nivolumab"),]
length(unique(clin.immuno$Tumor_Sample_Barcode))#25
clinEF<-clinAB[clinAB$Tumor_Sample_Barcode %in% clinCD$Tumor_Sample_Barcode,]
#write.csv(clinE,file="TCGA.SKCM.targeted.drug.csv",row.names = F)
table(clinCD$drug_name)
table(clinEF$drug_name)
#####
clin.target.drug <- read.csv("/public/home/lorihan/lrh/NSCLC/NSCLC_cnv/TCGA-LUNG.targeted.drug.csv")
colnames(clin.target.drug)

######
clin.chemo.drug<-clinAB
unique(clin.chemo.drug$Tumor_Sample_Barcode)#16
table(clin.chemo.drug$therapy_types)
clin.chemo.drug[clin.chemo.drug$therapy_types=="",]
clin.chemo.drug<-clin.chemo.drug[clin.chemo.drug$therapy_types!="",]
table(clin.chemo.drug$Stage)
table(clin.chemo.drug$M.stage)
str(clin.chemo.drug)##[ 145  18]
clin.chemo.drug[,-(1:5)]<-lapply(clin.chemo.drug[,-(1:5)], as.character)
table(clin.chemo.drug$measure_of_response)
table(clin.chemo.drug$Response)
#####
colnames(clin.chemo.drug)
table(is.na(clin.chemo.drug$days_to_drug_therapy_start))
table(is.na(clin.chemo.drug$days_to_drug_therapy_start)& is.na(clin.chemo.drug$days_to_drug_therapy_end))
na.time <- is.na(clin.chemo.drug$days_to_drug_therapy_start)&is.na(clin.chemo.drug$days_to_drug_therapy_end)#& clin.chemo.drug$drug_name=="unknown"
table(na.time)
clin.chemo.drug[na.time,]
na.sample <- clin.chemo.drug$Tumor_Sample_Barcode[na.time]
clinAB[clinAB$Tumor_Sample_Barcode %in% na.sample,c("bcr_patient_barcode","therapy_types","drug_name","days_to_drug_therapy_end")]
na.patient <- clin.chemo.drug$bcr_patient_barcode[na.time]
clinical.drug[clinical.drug$bcr_patient_barcode %in% na.patient,]
clin.chemo.drug <- clin.chemo.drug[!na.time,]

clin.chemo.drug[is.na(clin.chemo.drug$days_to_drug_therapy_start),]$days_to_drug_therapy_start="0"
##一线治疗：
temp = by(clin.chemo.drug,clin.chemo.drug$bcr_patient_barcode,function(x) rownames(x)[which.min(x$days_to_drug_therapy_start)])
skip = as.character(temp)
clin.chemo.drug=clin.chemo.drug[rownames(clin.chemo.drug) %in% skip,]
na.omit(match(skip,rownames(clin.chemo.drug)))
clin.chemo.drug=clin.chemo.drug[na.omit(match(skip,rownames(clin.chemo.drug))),]
table(clin.chemo.drug$drug_name)
table(clin.chemo.drug$therapy_types)
unique(na.sample)
unique(clin.chemo.drug$Tumor_Sample_Barcode)#143
table(clin.chemo.drug$Tumor_Sample_Barcode %in% na.sample)
#####
table(clin.chemo.drug$Response)
clin.chemo.drug$Response <- clin.chemo.drug$measure_of_response
clin.chemo.drug$Response<-gsub("Remission","Response",clin.chemo.drug$Response)
table(clin.chemo.drug$measure_of_response)
clin.chemo.drug[clin.chemo.drug$Response=="Unknown",]$Response=
  clin.chemo.drug[clin.chemo.drug$Response=="Unknown",]$measure_of_response
clin.chemo.drug[clin.chemo.drug$measure_of_response=="Unknown",]$measure_of_response=
  clin.chemo.drug[clin.chemo.drug$measure_of_response=="Unknown",]$Response

table(clin.chemo.drug$measure_of_response)
table(clin.chemo.drug$Response)
clin.chemo.drug[clin.chemo.drug$measure_of_response=="Unknown",]$DFS_STATUS
Response<-clin.chemo.drug$measure_of_response
table(Response)
clin.chemo.survival <- SKCM.survival[match(clin.chemo.drug$Tumor_Sample_Barcode,SKCM.survival$sample),c("sample","PFI.time","OS.time")]
chemo.logi.1<- Response=="Unknown" & clin.chemo.drug$days_to_drug_therapy_start < clin.chemo.survival$PFI.time
chemo.logi.2<- Response=="Unknown" & clin.chemo.drug$days_to_drug_therapy_start >= clin.chemo.survival$PFI.time 
table(chemo.logi.2)
table(chemo.logi.1)
clin.chemo.drug[chemo.logi.2,]
Response[chemo.logi.1] <- ifelse(clin.chemo.drug$DFS_STATUS[chemo.logi.1]=="Relapse-free","Stable Disease","Progressive Disease")
table(Response)
Response[chemo.logi.2] <- ifelse(clin.chemo.drug$OS_STATUS[chemo.logi.2]=="Alive","Stable Disease","Progressive Disease")
table(Response)
clin.chemo.drug$measure_of_response <- Response
str(clin.chemo.drug)
#####
table(is.na(clin.chemo.survival$OS.time))
table(is.na(clin.chemo.survival$PFI.time))
# clin.chemo.survival$PFI.time[is.na(clin.chemo.survival$PFI.time)]=0
# clin.chemo.survival$OS.time[is.na(clin.chemo.survival$OS.time)]=0
clin.chemo.drug[,c("days_to_drug_therapy_start","days_to_drug_therapy_end")] <- 
  lapply(clin.chemo.drug[,c("days_to_drug_therapy_start","days_to_drug_therapy_end")], as.numeric)

TTP <- clin.chemo.drug$days_to_drug_therapy_end-clin.chemo.drug$days_to_drug_therapy_start
table(is.na(TTP))#27
table(TTP==0)
table(clin.chemo.drug$DFS_STATUS)
table(clin.chemo.drug$Response)

TTP[clin.chemo.drug$measure_of_response =="Progressive Disease"]
TTP.logi.1 <-  clin.chemo.drug$measure_of_response !="Progressive Disease"  & clin.chemo.drug$days_to_drug_therapy_start >= clin.chemo.survival$PFI.time
table(TTP.logi.1)#T: 17
clin.chemo.drug[TTP.logi.1,]

TTP.logi.2 <- clin.chemo.drug$measure_of_response!="Progressive Disease"  & clin.chemo.drug$days_to_drug_therapy_start > clin.chemo.survival$OS.time
clin.chemo.survival[TTP.logi.2,]
clin.chemo.drug[TTP.logi.2,]
table(TTP.logi.2)#T: 2

TTP.logi.3 <-  clin.chemo.drug$days_to_drug_therapy_start < clin.chemo.survival$PFI.time #clin.chemo.drug$measure_of_response!="Progressive Disease" & clin.chemo.drug$DFS_STATUS =="Relapse" 
table(TTP.logi.3)#T: 267-->256 
clin.chemo.drug[TTP.logi.3,]
clin.chemo.survival[TTP.logi.3,]

TTP.logi.4 <- clin.chemo.drug$measure_of_response!="Progressive Disease" &  clin.chemo.drug$DFS_STATUS =="Relapse" & clin.chemo.drug$days_to_drug_therapy_start < clin.chemo.survival$PFI.time 
table(TTP.logi.4)#T:60-->59
clin.chemo.drug[TTP.logi.4,]
clin.chemo.drug$measure_of_response[TTP.logi.4] = "Progressive Disease"
table(clin.chemo.drug$measure_of_response)
TTP.logi.5 <- clin.chemo.drug$measure_of_response=="Progressive Disease" & clin.chemo.drug$days_to_drug_therapy_start < clin.chemo.survival$PFI.time & clin.chemo.drug$DFS_STATUS !="Relapse" 
table(TTP.logi.5)#T: 1
clin.chemo.drug[TTP.logi.5,]
clin.chemo.drug$measure_of_response[TTP.logi.5] = "Stable Disease"
TTP[TTP.logi.5]

####
TTP.logi.6 <- clin.chemo.drug$measure_of_response=="Progressive Disease" & clin.chemo.drug$days_to_drug_therapy_start >= clin.chemo.survival$PFI.time & is.na(clin.chemo.drug$days_to_drug_therapy_end)#&  clin.chemo.drug$OS_STATUS =="Dead" 
clin.chemo.drug[TTP.logi.6,]
TTP[TTP.logi.6]

TTP[TTP.logi.1] = (clin.chemo.survival$OS.time-clin.chemo.drug$days_to_drug_therapy_start)[TTP.logi.1]

TTP[TTP.logi.2] = (clin.chemo.drug$days_to_drug_therapy_end-clin.chemo.drug$days_to_drug_therapy_start)[TTP.logi.2]

TTP[TTP.logi.3] = (clin.chemo.survival$PFI.time-clin.chemo.drug$days_to_drug_therapy_start)[TTP.logi.3]
TTP

TTP[TTP.logi.6] = (clin.chemo.survival$OS.time-clin.chemo.drug$days_to_drug_therapy_start)[TTP.logi.6]

clin.chemo.drug[TTP.logi.6,]
(clin.chemo.survival$PFI.time-clin.chemo.drug$days_to_drug_therapy_start)[TTP.logi.6]

table(is.na(TTP))#2
table(TTP==0)
clin.chemo.drug[TTP==0,]
clin.chemo.drug$TTP <- TTP
table(clin.chemo.drug$measure_of_response)

table(clin.chemo.drug$Tumor_Sample_Barcode%in% clin.target.drug$Tumor_Sample_Barcode)
table(clin.chemo.drug[clin.chemo.drug$Tumor_Sample_Barcode%in% clin.target.drug$Tumor_Sample_Barcode,]$Stage)
table(clin.target.drug$Stage)

######
colnames(clin.chemo.drug)
paste(colnames(clin.chemo.drug),collapse = "','")
clin.drug <- clin.chemo.drug[,c('Tumor_Sample_Barcode','DFS_STATUS','DFS.month','OS_STATUS','OS.month','Histology','Sex','Age','T.stage','N.stage','M.stage','Stage',
                                'bcr_patient_barcode','therapy_types','drug_name','measure_of_response', 'days_to_drug_therapy_start','days_to_drug_therapy_end','TTP')]
colnames(clin.drug)<-colnames(clin.target.drug)[1:20][-9]
head(clin.target.drug)
table(clin.drug$Outcome)
clin.drug$Outcome[clin.drug$Outcome=="Progressive Disease"]="PD"
clin.drug$Outcome[clin.drug$Outcome=="Complete Response"]="CR"
clin.drug$Outcome[clin.drug$Outcome=="Stable Disease"]="SD"
clin.drug$Outcome[clin.drug$Outcome=="Partial Response"]="PR"

#####
dim(clin.drug)#143
table(clin.drug$Outcome)
table(clin.drug$Therapy_types)
setdiff(unique(clin_meta$Tumor_Sample_Barcode),clin.drug$Tumor_Sample_Barcode)#T:3 
#####
setdiff(clin.drug$Tumor_Sample_Barcode,na.sample)
intersect(clin.drug$Tumor_Sample_Barcode,na.sample)

na.sample <- unique(na.sample[ !( na.sample %in%clin.drug$Tumor_Sample_Barcode)])
unique(clin_meta$bcr_patient_barcode)
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% clin.drug$Tumor_Sample_Barcode)

clin.undrug <- clin_res[!(clin_res$bcr_patient_barcode %in% clin_meta$bcr_patient_barcode), ]
dim(clin.undrug)#325
clin.undrug <- clin_res[!(clin_res$bcr_patient_barcode %in% clin.drug$bcr_patient_barcode), ]
dim(clin.undrug)#331
unique(clin.undrug$bcr_patient_barcode)#325
clin.undrug <- clin.undrug[match(unique(clin.undrug$bcr_patient_barcode),clin.undrug$bcr_patient_barcode),]
intersect(clin.undrug$Tumor_Sample_Barcode,na.sample)
intersect(clin.undrug$bcr_patient_barcode,clin.drug$bcr_patient_barcode)
intersect(clin.undrug$bcr_patient_barcode,unique(clin_meta$bcr_patient_barcode))

clin.undrug <- clin.undrug[,c("Tumor_Sample_Barcode","DFS_STATUS","DFS.month","OS_STATUS","OS.month",
                              "Histology","Sex","Age","T.stage","N.stage",
                              "M.stage","Stage","days_to_drug_therapy_start","days_to_drug_therapy_end",
                              "bcr_patient_barcode","therapy_types","drug_name","measure_of_response")]
table(clin.undrug$days_to_drug_therapy_start)
clin.undrug$therapy_types <- as.character(clin.undrug$therapy_types)
clin.undrug$therapy_types[clin.undrug$days_to_drug_therapy_start!="Unknown"] <- "Unspecific"
clin.undrug$therapy_types[is.na(clin.undrug$days_to_drug_therapy_start)] <- "Unspecific"

clin.undrug$DFS_STATUS<-as.character(clin.undrug$DFS_STATUS)
clin.undrug$OS_STATUS<-as.character(clin.undrug$OS_STATUS)

clin.undrug[clin.undrug$measure_of_response=="Unknown",]$DFS_STATUS
Response <- as.character(clin.undrug$measure_of_response)
clin.undrug.survival <- SKCM.survival[match(clin.undrug$Tumor_Sample_Barcode,SKCM.survival$sample),c("sample","PFI.time","OS.time")]
table(Response)
Response[clin.undrug$measure_of_response=="Unknown"] <- ifelse(!is.na(clin.undrug$DFS_STATUS[clin.undrug$measure_of_response=="Unknown"]),
                                                    ifelse(clin.undrug$DFS_STATUS[clin.undrug$measure_of_response=="Unknown"]=="Relapse-free",
                                                           "Stable Disease","Progressive Disease"),"Unknown")
table(Response)
res.logi.1 <- Response!="Unknown" & Response!="Progressive Disease" & clin.undrug$DFS_STATUS =="Relapse"
table(res.logi.1)
Response[res.logi.1]="Progressive Disease"
clin.undrug[res.logi.1,]

res.logi.2 <- Response!="Unknown" & Response=="Progressive Disease" & clin.undrug$DFS_STATUS !="Relapse"
table(res.logi.2)
Response[res.logi.2]="Stable Disease"
clin.undrug[res.logi.2,]

table(Response)
clin.undrug$Response <- Response
str(clin.undrug)
clin.undrug$TTP<-clin.undrug.survival$PFI.time
clin.undrug[,c("DFS.month","TTP")]

colnames(clin.undrug)
colnames(clin.chemo.drug)
colnames(clin.drug)
clin.undrug <- clin.undrug[,c('Tumor_Sample_Barcode','DFS_STATUS','DFS.month','OS_STATUS','OS.month','Histology','Sex','Age','T.stage','N.stage','M.stage','Stage',
                              'bcr_patient_barcode','therapy_types','drug_name','Response', 'days_to_drug_therapy_start','days_to_drug_therapy_end','TTP')]
colnames(clin.undrug)<-colnames(clin.target.drug)[1:20][-9]
head(clin.target.drug)
table(clin.undrug$Outcome)
clin.undrug$Outcome[clin.undrug$Outcome=="Progressive Disease"]="PD"
clin.undrug$Outcome[clin.undrug$Outcome=="Complete Remission"]="CR"
clin.undrug$Outcome[clin.undrug$Outcome=="Stable Disease"]="SD"
clin.undrug$Outcome[clin.undrug$Outcome=="Partial Remission"]="PR"

######
clin.all<-rbind(clin.undrug,clin.drug)
table(clin.all$Therapy_types)
table(clin.all$Stage)
clin.all[clin.all$Stage=="",c("T.stage","N.stage","M.stage","Stage")]

######
stage.x<-SKCM.pheno[SKCM.pheno$sampleID %in% clin.all$Tumor_Sample_Barcode,c("sampleID","pathologic_T","pathologic_N","pathologic_M","pathologic_stage")]
stage.y<-SKCM.pheno[SKCM.pheno$sampleID %in% clin.all$Tumor_Sample_Barcode,]
stage<-as.character(stage.x$pathologic_stage)
names(stage)=stage.x$sampleID
table(stage)
stage.x[stage=="Stage II",]
stage.x[grep("Stage II",stage),]
table(stage)
stage.x[stage=="[Discrepancy]",]
stage.x[stage=="",]
which(stage=="")
intersect <- intersect(which(stage==""), grep("T4",stage.x$pathologic_T))
stage.y[intersect,]
stage[intersect] <- ifelse(stage.x[intersect,]$pathologic_M !="M1",
                           ifelse(stage.x[intersect,]$pathologic_N == "N0",
                                  "Stage IIB","Stage III"), "Stage IV")
intersect <- intersect(which(stage==""),grep("T3",stage.x$pathologic_T))
stage.y[intersect,]
#stage[intersect]<-ifelse(stage.x[intersect,]$pathologic_N == "N0","IIB","III")

stage.x[stage=="",]
stage[stage==""]<-ifelse(stage.x[stage=="",]$pathologic_N != "N1",
                  ifelse(stage.x[stage=="",]$pathologic_T == "T2b","Stage IIA","Stage I"),"Stage IIB")

table(stage)
stage[stage==""]="Unknown"
stage.x[stage=="Unknown",]
# stage[intersect(grep("Unknown",stage), grep("T2b",stage.x$pathologic_T))]<-ifelse(
#   stage.x[intersect(grep("Unknown",stage), grep("T2b",stage.x$pathologic_T)),]$pathologic_N == "N0","IIB","IIA")
stage.0 <- clin.all$Stage
table(stage.0)
table(stage)
stage[grep("Stage 0",stage)]="0" 
stage[grep("I/II NOS",stage)]="I/II" 
table(stage)
stage[stage=="Stage IV"]="III/IV" 
stage[stage=="Stage III"]="III/IV" 
stage[stage=="Stage II"]="I/II" 
stage[stage=="Stage I"]="I/II" 
stage[grep("Stage III",stage)]="III/IV"
stage[grep("Stage II",stage)]="I/II"
stage[grep("Stage I",stage)]="I/II"
table(stage)
clin.all$Stage<-stage[match(clin.all$Tumor_Sample_Barcode,names(stage))]
table(clin.all$Tumor_Sample_Barcode %in% SKCM.pheno[SKCM.pheno$history_of_neoadjuvant_treatment=="No",]$sampleID)
table(SKCM.pheno[SKCM.pheno$sampleID %in% clin.all$Tumor_Sample_Barcode,]$history_of_neoadjuvant_treatment)

clin.all[,c("DFS.month","OS.month","TTP")] <- lapply(clin.all[,c("DFS.month","OS.month","TTP")],as.numeric)
clin.all[,-match(c("DFS.month","OS.month","TTP"),colnames(clin.all))] <- lapply(clin.all[,-match(c("DFS.month","OS.month","TTP"),colnames(clin.all))],as.character)
table(clin.all$Stage)
table(clin.all$Sex)
#write.table(clin.all,"TCGA_LUNG.clin_meta.txt",sep="\t",row.names = F,quote = F)
###
#####CLINICAL ALL INFO#####
library(data.table)
id_mapped <- read.table("gencode.v22.annotation.gene.probeMap",header = T)
count<-fread("TCGA.SKCM.sampleMap%2FHiSeqV2.gz",sep = '\t',header = TRUE)
#count<-fread("TCGA-SKCM.htseq_counts.tsv.gz",sep = '\t',header = TRUE)
count=as.data.frame(count)
count[1:4,1:4] 
rownames(count)=count[,1]
count=count[,-1]
#genes=rownames(count)
count[1:4,1:4]
count=2^count-1
count[1:4,1:4] 
exprSet<-count
exprSet[1:4,1:4] 
##nonstaged
tmp=apply(exprSet,1,function(x){
  sum(x==0) < 10
}) #CYP2D6 was filtered because of this step earlier
table(tmp)
save=exprSet[tmp,]
save[1:5,1:4]
save = ceiling(save)
dim(save)
dim(exprSet)
#########
immune_checkpoint<-c("CD274","PDCD1LG2","PDCD1","CTLA4","HAVCR2","LAG3","VTCN1","TNFRSF8")
##HAVCR2 for TIM3
Tcell_receptor<-c("CD27", "GRAP2", "LCK", "PTPRCAP", "CCL5", "IL2RB", "C10orf54",
                  "IKZF3", "CD3G", "CD74", "CD3D", "CD8A", "CD4", "TIGIT")
TME<-c("IDO1", "PTGS2", "IL1B", "IL18", "IL6", "IL12A", "TNF", "NT5E") 
##:NT5E for CD73
match(immune_checkpoint,rownames(exprSet))
match(Tcell_receptor,rownames(exprSet))

T_effector_INFγ_gene<-c("GBP1", "IFI16", "IFI30", "IFNG", "IRF1", "STAT1", "TAP1", "TAP2",
                        "FAS", "PSMB9", "IL15RA", "GZMA", "GZMB", "EOMES", "CXCL10", "CXCL9",
                        "CXCL11", "TBX21", "PRF1")
antigen_presentation<-c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-DPA1","HLA-DPB1","HLA-DRA","B2M")

immune_gene<-c(immune_checkpoint,Tcell_receptor,TME,T_effector_INFγ_gene,antigen_presentation)
table(immune_gene %in% rownames(exprSet))
table(immune_gene %in% rownames(save))

immune_gene_set <- read.csv("/public/home/lorihan/lrh/Neoantigen_NSCLC/immume_gene_set.csv",header = T)
table(immune_gene_set$Metagene %in% rownames(exprSet))
table(immune_gene_set$Metagene %in% rownames(save))

filtered.immune <- immune_gene[!(immune_gene %in% rownames(save))]
#exprSet=exprSet[tmp,]
exprSet[1:4,1:4]
exprSet[filtered.immune,1:10]
##########
exprSet[1:5,1:4]
exprSet <- ceiling(exprSet)
#exprSet <- rbind(save,exprSet[filtered,])
remove(save);remove(count)
expr<-as.matrix(exprSet)
library(stringr)
colnames(expr)<-str_replace_all(colnames(expr),"[.]","-")
dim(expr)#20530   474
expr[1:5,1:5]    
class(expr)
num1<-which(as.numeric(substr(colnames(expr),14,15))>=10)
num2<-which(as.numeric(substr(colnames(expr),14,15))<10)
expr1<-as.matrix(expr[,num1])
colnames(expr1) <- colnames(expr)[num1]#"TCGA-GN-A4U8-11"
## expr1<-expr[,as.numeric(substr(colnames(expr),14,15))>10]
expr2<-expr[,num2]
length(unique(SKCM.pheno[,1]))##481
remove(expr)
#####
## ExprSet
exprSet<-cbind(expr1,expr2)
sample<-substr(colnames(exprSet),1,12)
######Normalization
library("edgeR")
require("limma")
num2<-which(as.numeric(substr(colnames(exprSet),14,15))<10)
y <- DGEList(counts=expr2)

keep <- rowSums(cpm(y)>1) >= 2#
y$samples$lib.size <- colSums(y$counts)

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
#y$samples
logcpm <- cpm(y, prior.count=2, log=TRUE)
logcpm[1:5,1:5]
dim(logcpm)
logcpm[1:5,1:5]
options(bitmapType='cairo')
#load("/public/home/lorihan/lrh/SKCM/SKCM_202203_result.Rdata")
exprSet[1:5,1:5];dim(exprSet)
logcpm[1:5,1:5];dim(logcpm)
#load("/public/home/lorihan/lrh/NSCLC/ensemble2name.Rdata")
#remove(expr1);remove(expr2)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
library(stringr)
gene_length <- read.csv("/public/home/lorihan/lrh/Neoantigen_NSCLC/gene_length_2.csv")
colnames(gene_length) <- c("gene_id","gene_length")
eff_length <- gene_length
length(unique(str_split(eff_length$gene_id,'[.]',simplify = T)[,1]))
eff_length <- eff_length[na.omit(match(id_mapped$id,eff_length$gene_id)),]
id_symbol <- id_mapped[match(eff_length$gene_id,id_mapped$id),]
eff_length[1:5,]
id_symbol[1:5,1:2]
eff_length$gene_symbol <- id_symbol$gene
table(rownames(exprSet) %in% eff_length$gene_symbol)#13101  
tpms <- exprSet[na.omit(match(eff_length$gene_symbol,rownames(exprSet))),]
eff_length <- eff_length[match(rownames(tpms),eff_length$gene_symbol),]
####
tpms <- apply(tpms,2,countToTpm, effLen=eff_length$gene_length)
tpms[1:3,1:5]
colSums(tpms)
table(colSums(tpms)=="1e+06")
#colSums(tpms)[colSums(tpms)!=1e+06]
#which(colSums(tpms)=="1e+06")
shapiro.test(tpms["CD8A",])
shapiro.test(logcpm["IFNG",])
expr_tpms <- log2(tpms+1)
library(limma) 
logtpm=normalizeBetweenArrays(expr_tpms)#expr_tpms#
# boxplot(expr_tpms,outline=FALSE, notch=T,col=group_neo, las=2)
# boxplot(logtpm,outline=FALSE, notch=T,col=group_neo, las=2)
ex <- logtpm
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
print("log2 transform unfinished")}else{print("log2 transform not needed")}
#boxplot(logtpm,col = "blue",xaxt = "n",outline = F)
shapiro.test(logtpm["CD8A",])
shapiro.test(logtpm["IFNG",])
# boxplot(logtpm[rownames(logtpm) %in% T_effector_INFγ_gene,1:100],col = "blue",xaxt = "n",outline = F)
# boxplot(logcpm[rownames(logcpm) %in% T_effector_INFγ_gene,1:100],col = "blue",xaxt = "n",outline = F)
boxplot(logtpm["IFNG",])
summary(logtpm["CD274",])
summary(logcpm["CD274",])
summary(exprSet["CD274",])
dim(exprSet)#20530   474
######
table(is.na(match(clin.all$Tumor_Sample_Barcode,SKCM.clin$submitter_id.samples)))#0
#Age <- CLIN$age_at_initial_pathologic_diagnosis[match(Neoantigen.data_lung$Samples,CLIN$sampleID)]
Age <- (-clin_SKCM$age_at_initial_pathologic_diagnosis[match(
  substr(clin.all$Tumor_Sample_Barcode,1,12),clin_SKCM$bcr_patient_barcode)])/365
table(is.na(Age))
Age <- (-SKCM.clin$days_to_birth[match(clin.all$Tumor_Sample_Barcode,
              substr(SKCM.clin$submitter_id.samples,1,15))])/365
table(is.na(Age))
table(Age[!is.na(Age)]>=65)
table(clin.all$Age)
colnames(clin.all)
##########
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% clin.all$Tumor_Sample_Barcode)
unmatch <- TCGA_SKCM_VCF[!substr(TCGA_SKCM_VCF$Sample,1,15) %in%
                           clin.all$Tumor_Sample_Barcode,]
clin.all[! clin.all$Tumor_Sample_Barcode %in% substr(TCGA_SKCM_VCF$Sample,1,15),]
clin.all$Tumor_Sample_Barcode[clin.all$bcr_patient_barcode %in% unmatch$Patients]
paste(clin.all$bcr_patient_barcode[clin.all$bcr_patient_barcode %in% unmatch$Patients],"06",sep = "-")
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% clin.all$Tumor_Sample_Barcode)
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% Dat$sampleID)
table(substr(TCGA_SKCM_VCF$Sample,1,15) %in% SKCM.pheno$sampleID)
TCGA_SKCM_VCF[!substr(TCGA_SKCM_VCF$Sample,1,15) %in% 
                clin.all$Tumor_Sample_Barcode,]

table(TCGA_SKCM_VCF$Sample %in% SKCM.clin$submitter_id.samples)
table(TCGA_SKCM_VCF$Patients %in% clin.all$bcr_patient_barcode)
TCGA_SKCM_VCF$Patients[!TCGA_SKCM_VCF$Patients %in% clin.all$bcr_patient_barcode]
table(TCGA_SKCM_VCF$Patients %in% SKCM.pheno$bcr_patient_barcode)
#
#SKCM_VCF <- read.table("/public/home/lorihan/lrh/SKCM/TCGA_SKCM_VCF_MuTect2.txt",header = T)
#match(SKCM_VCF$file_name,TCGA_SKCM_VCF$file_name)
HLA.type <- read.table("/public/home/lorihan/lrh/SKCM/TCGA-SKCM-hlaTypesAll.tsv",header = T,fill = T)
dim(HLA.type)#469   8  
dim(TCGA_SKCM_VCF)
match(TCGA_SKCM_VCF$Patients,HLA.type$patientBarcode)
table(TCGA_SKCM_VCF$Patients %in% HLA.type$patientBarcode)
TCGA_SKCM_VCF$Patients[!TCGA_SKCM_VCF$Patients %in% HLA.type$patientBarcode]
##"TCGA-GN-A269"
TCGA_SKCM_VCF <- TCGA_SKCM_VCF[na.omit(match(HLA.type$patientBarcode,
                                             TCGA_SKCM_VCF$Patients)),]
HLA.type <- HLA.type[HLA.type$patientBarcode %in% TCGA_SKCM_VCF$Patients,]
TCGA_SKCM_VCF[1:5,]
HLA.type[1:5,]
table(HLA.type$patientBarcode %in% TCGA_SKCM_VCF$Patients)
table(HLA.type$patientBarcode %in% clin.all$bcr_patient_barcode)
unique(HLA.type$patientBarcode[!HLA.type$patientBarcode %in% TCGA_SKCM_VCF$Patients])
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

match(TCGA_SKCM_VCF$Sample,paste(HLA.type$patientBarcode,"01A",sep = "-"))
HLA.type$patientBarcode <- TCGA_SKCM_VCF$Sample #paste(HLA.type$patientBarcode,"01A",sep = "-")
HLA.type[1:5,1:5]
dim(HLA.type)#469
length(unique(TCGA_SKCM_VCF$Patients))#469
table(TCGA_SKCM_VCF$Sample %in% HLA.type$patientBarcode)#469
TCGA_SKCM_VCF$Sample[is.na(match(TCGA_SKCM_VCF$Sample,HLA.type$patientBarcode))]
write.table(HLA.type,sep = "\t",row.names = F,col.names = F,
          quote = F,"/public/home/lorihan/lrh/SKCM/HLA_I.txt")

setwd("/public/home/lorihan/lrh/SKCM/")
filename <- list.files("/public/home/lorihan/lrh/SKCM/skcm_mutect2/")
table(TCGA_SKCM_VCF$file_id %in% filename)
TCGA_SKCM_VCF <- TCGA_SKCM_VCF[TCGA_SKCM_VCF$file_id %in% filename,]
table(TCGA_SKCM_VCF$Patients %in% clin.immuno$bcr_patient_barcode)
#load("Mutect2_SKCM.mutation.Rdata")
#unique(maf$Tumor_Sample_Barcode)##467
########
write.table(TCGA_SKCM_VCF[,c("Sample","Barcode","file_id","file_name","Matched","file_size")],"/public/home/lorihan/lrh/SKCM/TCGA_SKCM_VCF_MuTect2.txt",quote = FALSE, row.names = FALSE,sep = "\t")
save(TCGA_SKCM_VCF,clin.all,exprSet,tpms, HLA.type,logtpm,logcpm,clin.immuno,
     file="/public/home/lorihan/lrh/SKCM/SKCM_202204_result.Rdata")