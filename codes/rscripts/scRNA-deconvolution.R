# title: "Sade-Feldman scRNA-seq dataset"
# author: "Ruihan Luo"
# date: "Jul 19th,2022"
setwd("/public/home/lorihan/lrh/scRNA-seq")
library(data.table)
scRNA_skcm <- fread("/public/home/lorihan/lrh/scRNA-seq/SKCM/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
                    stringsAsFactor=F, header=T,sep = "\t",fill = T)
scRNA_skcm <- as.data.frame(scRNA_skcm)
scRNA_skcm[1:5,1:5]
dim(scRNA_skcm)
length(scRNA_skcm[,1])#23686
hg=scRNA_skcm[,1]#rownames(scRNA_skcm)#
scRNA_skcm=scRNA_skcm[,-ncol(scRNA_skcm)]
scRNA_skcm=scRNA_skcm[,-1]
treatment <- scRNA_skcm[1,]
unique(as.character(treatment[1,]))
dim(scRNA_skcm)
scRNA_skcm=scRNA_skcm[-1,]
rownames(scRNA_skcm)=hg[-1]
scRNA_skcm[1:5,1:5]
hg[grepl('^MT-',hg)]
meta=as.data.frame(colnames(scRNA_skcm))
colnames(meta)=c('title')
head(meta)
annot <- read.csv("SKCM/GSE120575_patient_ID.txt",header = T,
                  sep = "\t",fill = T, skip = 1,strip.white = TRUE, blank.lines.skip = TRUE,
                  comment.char = "")
annot <- annot[,1:7]
annot[1:5,1:5]
meta = dplyr::left_join(meta, annot,by='title')
rownames(meta)=colnames(scRNA_skcm)
table(is.na(match(annot$title,colnames(scRNA_skcm))))
markergene <- read.csv('SKCM/markergene.csv', header = T,
                       sep = ",",fill = T, skip = 0,strip.white = TRUE, blank.lines.skip = TRUE,
                       comment.char = "")
ref.split<-NULL
ref.split.new<-NULL
for(i in seq(5,55,5)){
  ref.split <- markergene[,(i-4):i]
  colnames(ref.split) <- colnames(markergene)[1:5]
  ref.split.new <<- rbind(ref.split.new,ref.split)
}
str(ref.split.new)
ref.split.new <- na.omit(ref.split.new)
unique(ref.split.new$CellType)
celltypes <- unique(ref.split.new$CellType)
annot_cluster <- read.csv("SKCM/GSE120575_clusterdata.txt",header = T,
                          sep = "\t",fill = T,strip.white = TRUE, blank.lines.skip = TRUE,
                          comment.char = "")
colnames(annot_cluster)[1] <- 'title'
names(celltypes) <- 1:11
annot_cluster$title[is.na(match(colnames(scRNA_skcm),annot_cluster$title))] <- 
  colnames(scRNA_skcm)[is.na(match(annot_cluster$title,colnames(scRNA_skcm)))]
annot_cluster$cluster <- celltypes[match(annot_cluster$Cluster.number,
                                         names(celltypes))]

CD8sub <- read.csv("SKCM/CD8fine_cluster.csv",header = T)
match(CD8sub$Cell.Name,annot_cluster$title)
table(CD8sub$Cell.Name %in% annot_cluster$title)
# colnames(CD8sub) <- c("title","CD8cluster")
annot_cluster$CD8cluster <- CD8sub$Cluster[match(annot_cluster$title,CD8sub$Cell.Name)]
meta = dplyr::left_join(meta, annot_cluster,by='title')
rownames(meta)=colnames(scRNA_skcm)
cd_genes <- read.table("SKCM/cd-genes.txt",sep = "\t",header=T)
skcm.seu <- scRNA_skcm[rownames(scRNA_skcm) %in% cd_genes$Gene.Name,]
table(is.na(match(annot_cluster$title,colnames(skcm.seu))))
table(is.na(match(rownames(meta),colnames(skcm.seu))))

dim(skcm.seu)
gkeep <- rowSums(skcm.seu > 4.5) >= 10;
table(gkeep)#T 16179; F 4466
gkeep.2 <- rowSums(skcm.seu > 12) >= 1;
table(gkeep.2)#T 20338 F 307
table(gkeep & gkeep.2)

skcm.seu <- skcm.seu[gkeep,]
sceset <- skcm.seu#scRNA_skcm[rownames(scRNA_skcm) %in% cd_genes$Gene.Name,]
sceset=apply(sceset,2,as.numeric)
#skcm.seu <- as.matrix(skcm.seu)
rownames(sceset) <- rownames(skcm.seu)
skcm.seu <- sceset
library(Seurat)
pbmc  <- CreateSeuratObject(counts = skcm.seu,#assays(sceset_mon_norm)$logcounts,#
                            meta.data = meta,
                            min.cells = 10, min.features = 200,project = 'Smart-seq2')

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")  #
table(pbmc@meta.data$percent.mt<5)
sum(pbmc[["RNA"]]@data[,1])
TPM <- expm1(pbmc[["RNA"]]@data)
head(colSums(TPM))
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 4000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes,verbose = FALSE)#check.for.norm=F
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),verbose = FALSE) 
DimPlot(object = pbmc, reduction = "pca")
ElbowPlot(pbmc)
library(harmony)
seu.harmony <- RunHarmony(pbmc, group.by.vars = "orig.ident")
dim.use <- 1:10
seu.harmony <- FindNeighbors(seu.harmony, dims = dim.use, reduction = "harmony")
seu.harmony <- FindClusters(seu.harmony, resolution = 0.4, reduction = "harmony")
seu.harmony <- FindClusters(seu.harmony, resolution = 0.6, reduction = "harmony")
seu.harmony <- RunTSNE(seu.harmony, dims = dim.use, do.fast = TRUE, reduction = "harmony")

# UMAP
Idents(seu.harmony)="RNA_snn_res.0.4"
table(Idents(seu.harmony))
harmony.markers <- FindAllMarkers(seu.harmony, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                                  min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                                  assay = "RNA")
saveRDS(seu.harmony,file = "SKCM/GSE120575_seu.harmony.rds")
save(harmony.markers,file="SKCM/GSE120575_harm_res.0.4_cluster_logfc0.Rdata")
DimPlot(seu.harmony, reduction = "tsne", group.by = "cluster", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
# ctrl <- pbmc
# save(ctrl,file='SKCM/tSNE_after_GSE120575.pbmc.Rdata')
save(scRNA_skcm,pbmc,seu.harmony,file = 'SKCM/tSNE_after_GSE120575.seurat.Rdata')
# load('SKCM/tSNE_after_GSE120575.seurat.Rdata')
# ------
Idents(pbmc) 
logFCfilter=0.25
adjPvalFilter=0.05
Idents(pbmc)="RNA_snn_res.0.5"
table(Idents(pbmc))
pbmc.markers <- FindAllMarkers(pbmc, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                               min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                               assay = "RNA")
save(pbmc.markers,file="SKCM/GSE120575_pbmc_res.0.5_cluster_logfc0.Rdata")
# ----
ctrl <- readRDS('SKCM/GSE120575_seu.harmony.rds')
table(ctrl$cluster)
table(ctrl$ref_gsea)
ctrl$cluster <- str_split(ctrl$cluster,"- ",simplify = T)[,2]
table(Idents(ctrl))
Idents(ctrl) <- "Cluster.number"
pdf("/public/home/lorihan/lrh/scRNA-seq/All.refclusters_GS120575_SKCM.pdf",height = 6,width = 9.5)
DimPlot(ctrl,reduction = "tsne",label=T, group.by = 'cluster')+ #NoAxes() + 
  ggtitle("GSE120575 (Melanoma)")
dev.off()
pdf("/public/home/lorihan/lrh/scRNA-seq/All.refclusters_nolabel_SKCM.pdf",height = 6,width = 8)
DimPlot(ctrl,reduction = "tsne",label=F, group.by = 'cluster')+ #NoAxes() + 
  ggtitle("GSE120575 (Melanoma)")
dev.off()

cowplot::plot_grid(ncol = 2, DimPlot(ctrl, label = T, group.by = "Cluster.number") + 
                     NoAxes(), DimPlot(ctrl, label = T, group.by = "RNA_snn_res.0.4") + NoAxes())

###########---CD8+ subclusters---#########
load('SKCM/tSNE_after_GSE120575.seurat.Rdata')
ctrl <- seu.harmony[,!is.na(seu.harmony$CD8cluster)]
ctrl <- FindVariableFeatures(ctrl, 
                             selection.method = "vst", nfeatures = 2000)  
# load("Lung/GSE179994_tSNE_seurat.Rdata") #from scRNA-lung.R
# UMAP_CD8_annot <- read.table('/public/home/lorihan/lrh/scRNA-seq/Lung/UMAP_CD8_3clusters.txt',header = T)
# colnames(UMAP_CD8_annot)[3] <- 'CD8_cluster' 
# CD8_annot <- annot[annot$celltype %in% "CD8",]
# CD8_annot <- cbind(CD8_annot,UMAP_CD8_annot)
# ctrl = pbmc[,!is.na(pbmc$CD8cluster)]
# seu.harmony <- readRDS('Lung/GSE179994_seu.harmony.rds')
# ctrl <- subset(seu.harmony, subset = celltype=="CD8")
# ctrl=CreateSeuratObject(counts = ctrl@assays$RNA@counts,
#                         meta.data = ctrl@meta.data) 

# ctrl <- NormalizeData(ctrl, normalization.method =  "LogNormalize",
#                       scale.factor = 1e4)
# GetAssay(ctrl,assay = "RNA")
# ctrl <- FindVariableFeatures(ctrl, 
#                              selection.method = "vst", nfeatures = 3000)  

ctrl <- ScaleData(ctrl) 
ctrl <- RunPCA(object = ctrl, pc.genes = VariableFeatures(ctrl)) 
# DimHeatmap(ctrl, dims = 1:12, cells = 100, balanced = TRUE)
# ElbowPlot(ctrl) 
ctrl <- FindNeighbors(ctrl, dims = 1:10)
ctrl <- FindClusters(ctrl, resolution = 0.4)
table(ctrl@meta.data$RNA_snn_res.0.4)  
set.seed(123)
ctrl <- RunTSNE(object = ctrl, dims = 1:10, do.fast = TRUE)
DimPlot(ctrl,reduction = "tsne",label=T)
ctrl <- RunUMAP(object = ctrl, dims = 1:10, do.fast = TRUE)
DimPlot(ctrl,reduction = "umap",label=T) 
head(ctrl@meta.data)
DimPlot(ctrl,reduction = "umap",label=T,group.by = 'sample')
ctrl_CD8T <- ctrl

.libPaths()  
myPaths <- .libPaths()  
new <- c('/public/apps/R-4.1.0/library','/public/home/lorihan/R/x86_64-pc-linux-gnu-library/4.1','/usr/local/lib64/R/library',
         '/public/home/lorihan/miniconda3/lib/R/library','/public/home/lorihan/R/x86_64-conda-linux-gnu-library/4.1')#
myPaths <- c(new,myPaths) 
.libPaths(myPaths) 
library(harmony)
library(Seurat)
ctrl <- ctrl_CD8T
seuratObj <- RunHarmony(ctrl, "sample")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:10, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=T )  
ctrl=seuratObj
remove(seuratObj)
ctrl <- FindNeighbors(ctrl, reduction = "harmony",
                      dims = 1:10)
ctrl <- FindClusters(ctrl, resolution = 0.1)
table(ctrl@meta.data$RNA_snn_res.0.1)
ctrl <- FindClusters(ctrl, resolution = 0.15)
table(ctrl@meta.data$RNA_snn_res.0.15)#seurat_clusters)
DimPlot(ctrl,reduction = "umap",label=T)
ctrl <- RunTSNE(object = ctrl, dims = 1:10, do.fast = TRUE)
DimPlot(ctrl,reduction = "tsne",label=T)
# ggsave(filename = 'harmony_umap_recluster_by_0.1.pdf') 
DimPlot(ctrl,reduction = "umap",label=T)#, group.by = 'orig.ident') 
# ggsave(filename = 'harmony_umap_ctrl_recluster_by_orig.ident.pdf') 
DimPlot(ctrl,reduction = "tsne",label=T, group.by = 'CD8cluster') 
DimPlot(ctrl_CD8T,reduction = "tsne",label=T, group.by = 'CD8cluster') 
ctrl.harmony <- ctrl
saveRDS(ctrl.harmony,file = "SKCM/CD8Tcluster.rds") 
save(ctrl_CD8T,ctrl.harmony,file = 'SKCM/CD8Tcluster_GSE120575.Rdata')
# save(ctrl_CD8T,ctrl.harmony,file = 'Lung/CD8Tcluster_newGSE179994.Rdata')
dim(ctrl.harmony)
meta <- ctrl.harmony@meta.data
meta <- dplyr::left_join(meta,CD8_annot,by='cellid')
rownames(meta)=colnames(ctrl.harmony)
ctrl.harmony@meta.data <-meta
markergene <- read.csv('SKCM/CD8markergene.csv', header = T,
                       sep = ",",fill = T, skip = 0,strip.white = TRUE, blank.lines.skip = TRUE,
                       comment.char = "")
ref.split<-NULL
ref.split.new<-NULL
for(i in seq(5,ncol(markergene),5)){
  ref.split <- markergene[,(i-4):i]
  colnames(ref.split) <- colnames(markergene)[1:5]
  ref.split.new <<- rbind(ref.split.new,ref.split)
}
str(ref.split.new)
ref.split.new <- na.omit(ref.split.new)
unique(ref.split.new$CellType)
ref.split.new[match(terminal_CD8,ref.split.new$GeneName),c(1,5)]
ref.split.new[match(progenitor_CD8,ref.split.new$GeneName),c(1,5)]
ref.split.CD8 <- ref.split.new
save(ref.split.CD8,ref.split.new,file="SKCM/markergene.Rdata")

########---scPred will train a classifier based on all principal components---########
#First, getFeatureSpace will create a scPred object stored in the @misc slot where it extracts the PCs that best separates the different celltypes. Then trainModel will do the actual training for each celltype.
# reference <- getFeatureSpace(reference, "cell_type")
# # ●  Extracting feature space for each cell type...
# ## DONE!
# reference <- trainModel(reference)
# ##  Training models for each cell type...
# ## maximum number of iterations reached 0.000116588 -0.0001156614DONE!
# #We can then print how well the training worked for the different celltypes by printing the number of PCs used for each, the ROC value and Sensitivity/Specificity.
# get_scpred(reference)
# ctrl <- scPredict(ctrl, reference)
# DimPlot(ctrl, group.by = "scpred_prediction", label = T, repel = T) + NoAxes()
# #Now plot how many cells of each celltypes can be found in each cluster.
# ggplot(ctrl@meta.data, aes(x = RNA_snn_res.0.5, fill = scpred_prediction)) + geom_bar() + 
#   theme_classic()
# crossTab(ctrl, "predicted.id", "scpred_prediction")

#--------method 2-------
# load('Lung/CD8Tcluster_GSE179994.Rdata')
load("SKCM/CD8Tcluster_GSE120575.Rdata")
ctrl <- ctrl.harmony
table(Idents(ctrl))
ctrl <- SetIdent(ctrl, value = "RNA_snn_res.0.4")
ctrl <- SetIdent(ctrl, value = "RNA_snn_res.0.1")
table(Idents(ctrl))
DimPlot(ctrl,reduction = "umap",label=T, group.by = 'cluster') 

DimPlot(ctrl,reduction = "tsne",label=T, group.by = 'CD8cluster') 
cowplot::plot_grid(ncol = 2, DimPlot(ctrl, reduction = "tsne",label = T, group.by = "CD8cluster") + 
                     NoAxes(), DimPlot(ctrl,reduction = "tsne", label = T, group.by = "RNA_snn_res.0.1") + NoAxes())

DGE_table <- FindAllMarkers(ctrl, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                            assay = "RNA")
save(DGE_table,file="SKCM/ctrl_logfc0.Rdata")
# save(DGE_table,file="Lung/ctrl_0.15_CD8Tall_logfc0.Rdata")
logFCfilter <- 0.25
adjPvalFilter <- 0.05
sig.markers=DGE_table[(abs(as.numeric(as.vector(DGE_table$avg_log2FC)))>
                         logFCfilter & as.numeric(as.vector(
                           DGE_table$p_val_adj))<adjPvalFilter),]
table(abs(as.numeric(as.vector(DGE_table$avg_log2FC)))>
        logFCfilter)
#table(DGE_table$p_val_adj<0.05)
sig.markers=DGE_table[abs(as.numeric(as.vector(DGE_table$avg_log2FC)))>
                        logFCfilter,]

load("SKCM/GSE120575_harm_res.0.4_cluster_logfc0.Rdata")
# load("Lung/GSE179994_harm_res.0.6_cluster_logfc0.Rdata")
DGE_table <- harmony.markers
DGE_table <- pbmc.markers
# split into a list
DGE_table <- sig.markers
load("Lung/ctrl_0.15_CD8Tall_logfc0.Rdata")
load("SKCM/ctrl_logfc0.Rdata")
DGE_list <- split(DGE_table, DGE_table$cluster)
#DGE_list <- split(sig.markers, sig.markers$cluster)
unlist(lapply(DGE_list, nrow))
##    0    1    2    3    4    5    6    7    8    9   10 
## 3153 2483 3394 2837 2573 3956 2150 3753 2465 2142 3342
# 0    1    2    3    4    5    6    7    8 
# 2311 3380 3177 1734 2848 3543 3194 1822 2858
# Compute differential gene expression in reference dataset (that has cell annotation)
# reference <- SetIdent(reference, value = "cell_type")
# 
# reference_markers <- FindAllMarkers(reference, min.pct = 0.1, min.diff.pct = 0.2, 
#                                     only.pos = T, max.cells.per.ident = 20, return.thresh = 1)
# 
# # Identify the top cell marker genes in reference dataset select top 50 with highest foldchange among top 100 signifcant genes.
# reference_markers <- reference_markers[order(reference_markers$avg_log2FC, decreasing = T), ]
# top50_cell_selection <- reference_markers %>% group_by(cluster) %>% top_n(-100, p_val) %>% 
#   top_n(50, avg_log2FC)
# 
# # Transform the markers into a list
# ref_list = split(top50_cell_selection$gene, top50_cell_selection$cluster)

########---GSEA on DEGs---##########
#####refering from SKCM marker gene----
load('SKCM/markergene.Rdata')
# ref_list = split(ref.split.new$GeneName, ref.split.new$CellType)
ref_list = split(ref.split.CD8$GeneName, ref.split.CD8$CellType)
unlist(lapply(ref_list, length))
suppressPackageStartupMessages(library(fgsea))
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea(pathways = ref_list, stats = gene_rank, nperm = 10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)
res_gsea <- res
res <- res_gsea
# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
  x[x$pval < 0.1, ]
})
res <- lapply(res, function(x) {
  x[x$size > 2, ]
})
res <- lapply(res, function(x) {
  x[order(x$NES, decreasing = T), ]
})
res

new.cluster.ids <- unlist(lapply(res, function(x) {
  as.data.frame(x)[1, 1]
}))
ctrl$ref_gsea <- new.cluster.ids[as.character(ctrl@active.ident)]
cowplot::plot_grid(ncol = 2, DimPlot(ctrl, label = T, group.by = "cluster") + 
                     NoAxes(), DimPlot(ctrl, label = T, group.by = "ref_gsea") + NoAxes())

ctrl.harmony$ref_gsea_CD8 <- new.cluster.ids[as.character(ctrl.harmony@active.ident)]
seu.harmony$ref_gsea <- new.cluster.ids[as.character(seu.harmony@active.ident)]
pdf("Cross_cluster_skcm_gsea_GSE179994_harm_res.0.6.pdf",height = 8,width = 18)
cowplot::plot_grid(ncol = 2, DimPlot(seu.harmony, label = T, group.by = "cluster") + 
                     NoAxes(), DimPlot(seu.harmony, label = T, group.by = "ref_gsea") + NoAxes())
dev.off()
# saveRDS(seu.harmony,file = "Lung/GSE179994_seu.harmony.rds")

cowplot::plot_grid(ncol = 2, DimPlot(ctrl, label = T, group.by = "cluster") + 
                     NoAxes(), DimPlot(ctrl, label = T, group.by = "predicted.id") + NoAxes())

cowplot::plot_grid(ncol = 3, DimPlot(ctrl, label = T, group.by = "ref_gsea") + NoAxes() + 
                     ggtitle("GSEA"), DimPlot(ctrl, label = T, group.by = "predicted.id") + NoAxes() + 
                     ggtitle("LabelTransfer"), DimPlot(ctrl, label = T, group.by = "scpred_prediction") + 
                     NoAxes() + ggtitle("scPred"))
########---download genesets from CellMarker---#######
# Download gene marker list
if (!dir.exists("reference_data/CellMarker_list/")) {
  dir.create("reference_data/CellMarker_list")
  download.file(url = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt", 
                destfile = "./reference_data/CellMarker_list/Human_cell_markers.txt")
  # download.file(url = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt", 
  #               destfile = "./data/CellMarker_list/Mouse_cell_markers.txt")
}
# Load the human marker table
SingleCellMarker <- read.delim("reference_data/CellMarker_list/Single_cell_markers.txt")
SingleCellMarker <- SingleCellMarker[SingleCellMarker$speciesType == "Human", ]
#SingleCellMarker <- SingleCellMarker[SingleCellMarker$cancerType == "Normal", ]
CellMarker <- read.delim("reference_data/CellMarker_list/all_cell_markers.txt")
CellMarker <- CellMarker[CellMarker$speciesType == "Human", ]
# CellMarker <- CellMarker[CellMarker$cancerType == "Normal", ]
# CellMarker <- CellMarker[CellMarker$tissueType == "Lung", ]

# Filter by tissue (to reduce computational time and have tissue-specific
# classification) sort(unique(markers$tissueType))
# grep('blood',unique(markers$tissueType),value = T) markers <- markers [
# markers$tissueType %in% c('Blood','Venous blood', 'Serum','Plasma',
# 'Spleen','Bone marrow','Lymph node'), ]
CD8Tex <- CellMarker[1:2,]
CD8Tex$cellName <- c("Progenitor Exhausted CD8 T cells","Terminally Exhausted CD8 T cells")
CD8Tex$tissueType <- c("skin","skin")
CD8Tex$cancerType <- c("Melanoma","Melanoma")
CD8Tex$cellType <- c("Cancer cell","Cancer cell")
CD8Tex$PMID <- c('30388456','30388456')

CD8Tex$geneSymbol[1] <- paste(progenitor_CD8,', ',collapse  = "")
CD8Tex$geneSymbol[2] <- paste(terminal_CD8, ', ',collapse  = "")
CD8Tex <- rbind(CD8Tex,CellMarker)
CD8Tex[grepl("CD8",CD8Tex$cellName),]
# CD8Tex <- CD8Tex[grepl("CD8",CD8Tex$cellName)|grepl("CD8A",CD8Tex$geneSymbol),]
# CD8Tex <- CD8Tex[grepl("CD8",CD8Tex$cellName),]
CD8Tex <- CD8Tex[grepl("T ",CD8Tex$cellName),]
# remove strange characters etc.
celltype_list <- lapply(unique(CD8Tex$cellName), function(x) {
  x <- paste(CD8Tex$geneSymbol[CD8Tex$cellName == x], sep = ",")
  x <- gsub("[[]|[]]| |-", ",", x)
  x <- unlist(strsplit(x, split = ","))
  x <- unique(x[!x %in% c("", "NA", "family")])
  x <- casefold(x, upper = T)
})
names(celltype_list) <- unique(CD8Tex$cellName)
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )

# celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
# celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 5]

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank, nperm = 10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)
res_gsea <- res
res <- res_gsea
# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
  x[x$pval < 0.05, ]
})
res <- lapply(res, function(x) {
  x[x$size > 5, ]
})
res <- lapply(res, function(x) {
  x[order(x$NES, decreasing = T), ]
})

# show top 3 for each cluster.
lapply(res, head, 3)
#可视化基因集富集分析注释到的细胞类型。

new.cluster.ids <- unlist(lapply(res, function(x) {
  as.data.frame(x)[1, 1]
}))
ctrl$cellmarker_gsea <- new.cluster.ids[as.character(ctrl@active.ident)]
table(ctrl$CD8cluster[is.na(ctrl$cellmarker_gsea)])
table(ctrl$RNA_snn_res.0.1[is.na(ctrl$cellmarker_gsea)])
table(ctrl$CD8cluster[ctrl$RNA_snn_res.0.1 %in% c(7,8)])
cowplot::plot_grid(ncol = 2, DimPlot(ctrl, reduction = "tsne",label = T, group.by = "RNA_snn_res.0.1") + 
                     NoAxes(), DimPlot(ctrl,reduction = "tsne", label = T, group.by = "cellmarker_gsea") + NoAxes())
ctrl$cellmarker_gsea[ctrl$RNA_snn_res.0.1 %in% c(5)] <- 'Exhausted CD8+ T cell'
ctrl$cellmarker_gsea[ctrl$RNA_snn_res.0.1 %in% c(7,8)] <- 'Effector CD8+ memory T (Tem) cell'
ctrl$cellmarker_gsea[ctrl$RNA_snn_res.0.1 %in% c(3)] <- 'Terminal Exhausted CD8 T cells'

cowplot::plot_grid(ncol = 2, DimPlot(ctrl, reduction = "tsne",label = T, group.by = "RNA_snn_res.0.1") + 
                     NoAxes(), DimPlot(ctrl,reduction = "tsne", label = T, group.by = "CD8cluster") + NoAxes())

#ctrl$cellmarker_gsea[is.na(ctrl$cellmarker_gsea)] <- 'Terminal Exhausted CD8 T cells'
save(DGE_table,celltype_list, ctrl,file = '/public/home/lorihan/lrh/scRNA-seq/SKCM/ctrl_CD8T_newcelltype.Rdata') 

load("SKCM/CD8Tcluster_GSE120575.Rdata")
load("/public/home/lorihan/lrh/scRNA-seq/SKCM/ctrl_CD8T_newcelltype.Rdata")
table(Idents(ctrl.harmony))
table(ctrl.harmony$cellmarker_gsea)
table(ctrl$CD8cluster)
table(Idents(ctrl))
table(ctrl$cellmarker_gsea)
cowplot::plot_grid(ncol = 2, DimPlot(ctrl, reduction = "tsne",label = T, group.by = "CD8cluster") + 
                     NoAxes(), DimPlot(ctrl, reduction = "tsne",label = T, group.by = "cellmarker_gsea") + NoAxes())
ctrl$cellmarker_gsea[ctrl$cellmarker_gsea=="Naive CD8+ T cell"] <- 'Effector CD8+ memory T (Tem) cell'
ctrl$cellmarker_gsea <- gsub(" Exhausted CD8"," exhausted CD8+",ctrl$cellmarker_gsea)
ctrl$cellmarker_gsea <- gsub("cells","cell",ctrl$cellmarker_gsea)
pdf("/public/home/lorihan/lrh/scRNA-seq/CD8.GSEA_GS120575_SKCM.old.pdf",height = 5,width = 7.7)
DimPlot(ctrl,reduction = "tsne",label=T, group.by = 'cellmarker_gsea')+ NoLegend()+#NoAxes() + 
  ggtitle("CD8+ T cell subclustering (Melanoma)")
dev.off()
pdf("/public/home/lorihan/lrh/scRNA-seq/CD8.GSEA_GSE120575_nolabel_SKCM.pdf",height = 5,width = 7.7)
DimPlot(ctrl,reduction = "tsne",label=F, group.by = 'cellmarker_gsea')+ #NoAxes() + 
  ggtitle("CD8+ T cell subclustering (Melanoma)")
dev.off()
library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)]
pal <- vc_cols[c(1,3,5,6,8,9)]
plot6 <- DimPlot(ctrl,reduction = "tsne",label=T, pt.size = 1,cols = pal,
                 group.by = 'cellmarker_gsea')+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "CD8+ T cell subclustering (Melanoma)") + 
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot6
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  

plot7 <- DimPlot(ctrl,reduction = "tsne", label = F, pt.size = 0.7,#label.size = 5,
                 group.by = 'cellmarker_gsea',cols = c(brewer.pal(7,"Set2")))+
  NoLegend()+labs(x = "tSNE 1", y = "tSNE 2",title =  "CD8+ T cell subclustering (Melanoma)") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
pdf("/public/home/lorihan/lrh/scRNA-seq/CD8.GSEA_GS120575_SKCM.pdf",height = 5,width = 6)
plot7
dev.off()

plot_grid(plot6,plot7)
#########---Prepare for Cibersort---#########
table(Idents(seu.harmony))
table(seu.harmony$cluster)
table(grepl("Exhausted",seu.harmony$cluster))
table(seu.harmony$cluster.new)
group.by <- "cluster" 
group.by <- "cluster.new" 
mat <- as.data.frame(t(as.matrix(GetAssayData(seu.harmony, assay = "RNA", slot = "data")))) 
mat <- aggregate(mat, by=list(seu.harmony@meta.data[[group.by]]), FUN="mean") 
group.by <- "CD8cluster" 
group.by <- "cellmarker_gsea" 
mat <- as.data.frame(t(as.matrix(GetAssayData(ctrl.harmony, assay = "RNA", slot = "data")))) 
mat <- as.data.frame(t(as.matrix(GetAssayData(ctrl, assay = "RNA", slot = "data")))) 
#mat <- expm1(x = mat) #using log1p， expm1 v.s. log1p
mat <- aggregate(mat, by=list(ctrl.harmony@meta.data[[group.by]]), FUN="mean") 
# -----
rownames(mat) <- mat$Group.1 
mat <- t(mat[,-1]) 
head(mat)
#----CD8T----
X.sub <- mat[rownames(mat) %in%  VariableFeatures(object = ctrl),]
dim(X.sub)
write.table(X.sub,"/public/home/lorihan/lrh/scRNA-seq/SKCM/ctrl.CD8.gsea.txt",sep = "\t",col.names = T,row.names = T)     #sig.txt
write.table(mat,"/public/home/lorihan/lrh/scRNA-seq/SKCM/sig.harmony.cluster.all.txt",sep = "\t",col.names = T,row.names = T) #sig.txt

###########---bulk RNAseq SKCM---#########
dim(logtpm)#18767   474
exp2 = as.data.frame(logtpm[match(unique(rownames(logtpm)),
                                  rownames(logtpm)),])
exp2 = rownames_to_column(exp2)
write.table(exp2,file = "/public/home/lorihan/lrh/SKCM/logtpm_skcm474.txt",row.names = F,quote = F,sep = "\t")

#!/usr/local/bin/Rscript
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library('e1071')            
source("/public/home/lorihan/R/Rscripts_test/CIBERSORTx.R")  #activate function

setwd("/public/home/lorihan/lrh/SKCM")
if(T){
  TME.results = CIBERSORT("/public/home/lorihan/lrh/scRNA-seq/SKCM/ctrl.CD8.gsea.txt",
                          "/public/home/lorihan/lrh/SKCM/logtpm_skcm474.txt" ,
                          perm = 1000,
                          QN = FALSE)
  save(TME.results,file = "/public/home/lorihan/lrh/SKCM/ciber_perm1000_skcm.gsea.CD8T.ctrl.Rdata")
}

if(T){
  TME.results = CIBERSORT("/public/home/lorihan/lrh/scRNA-seq/SKCM/sig.harmony.cluster.all.txt",
                          "/public/home/lorihan/lrh/SKCM/logtpm_skcm474.txt" ,
                          perm = 1000,
                          QN = FALSE)
  save(TME.results,file = "/public/home/lorihan/lrh/SKCM/ciber_perm1000_GSE120575_cluster.all.new.Rdata")
}
