library(Seurat)
library(cowplot)
setwd('/lustre/fyh/lab_other_work/singlecell/output/workspace')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
library(dplyr)

BATCH_deci1_all=rep('deci1',ncol(DATA)) 
BATCH_deci1_all[c((ncol(deci1)+1):(ncol(deci1)+ncol(deci2)))]='deci2'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)))]='deci3'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)))]='deci4'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)))]='deci5'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)))]='deci6'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)))]='deci7'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)))]='deci8'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)))]='deci9'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)))]='pla1'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)))]='pla2'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)))]='pla3'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)))]='pla4'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+ncol(pla5)))]='pla5'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+ncol(pla5)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+ncol(pla5)+ncol(pla6)))]='pla6'
BATCH_deci1_all[c((ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+ncol(pla5)+ncol(pla6)+1):(ncol(deci1)+ncol(deci2)+ncol(deci3)+ncol(deci4)+ncol(deci5)+ncol(deci6)+ncol(deci7)+ncol(deci8)+ncol(deci9)+ncol(pla1)+ncol(pla2)+ncol(pla3)+ncol(pla4)+ncol(pla5)+ncol(pla6)+ncol(pla7)))]='pla7'
deci1_1.data<- deci1.data
deci2_1.data<- deci2.data
deci3_1.data<- deci3.data
deci4_1.data<- deci4.data
deci5_1.data<- deci5.data
deci6_1.data<- deci6.data
deci7_1.data<- deci7.data
deci8_1.data<- deci8.data
deci9_1.data<- deci9.data
pla1_1.data<- pla1.data
pla2_1.data<- pla2.data
pla3_1.data<- pla3.data
pla4_1.data<- pla4.data
pla5_1.data<- pla5.data
pla6_1.data<- pla6.data
pla7_1.data<- pla7.data

colnames(deci1_1.data)=paste0('deci1_',colnames(deci1_1.data))
colnames(deci2_1.data)=paste0('deci2_',colnames(deci2_1.data))
colnames(deci3_1.data)=paste0('deci3_',colnames(deci3_1.data))
colnames(deci4_1.data)=paste0('deci4_',colnames(deci4_1.data))
colnames(deci5_1.data)=paste0('deci5_',colnames(deci5_1.data))
colnames(deci6_1.data)=paste0('deci6_',colnames(deci6_1.data))
colnames(deci7_1.data)=paste0('deci7_',colnames(deci7_1.data))
colnames(deci8_1.data)=paste0('deci8_',colnames(deci8_1.data))
colnames(deci9_1.data)=paste0('deci9_',colnames(deci9_1.data))
colnames(pla1_1.data)=paste0('pla1_',colnames(pla1_1.data))
colnames(pla2_1.data)=paste0('pla2_',colnames(pla2_1.data))
colnames(pla3_1.data)=paste0('pla3_',colnames(pla3_1.data))
colnames(pla4_1.data)=paste0('pla4_',colnames(pla4_1.data))
colnames(pla5_1.data)=paste0('pla5_',colnames(pla5_1.data))
colnames(pla6_1.data)=paste0('pla6_',colnames(pla6_1.data))
colnames(pla7_1.data)=paste0('pla7_',colnames(pla7_1.data))
> pdf('2.pdf')
> DimPlot(pbmc_batch, reduction.use='umap', group.by='orig.ident', pt.size=0.1)
> dev.off()
> pdf('deci1.pdf',width=16,height=16)
> DimPlot(pbmc_batch, reduction = "umap", split.by = "orig.ident",pt.size=2, label.size=3,do.label=T,ncol = 4)
> dev.off()

pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
pdf('3.pdf')
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)  
dev.off() 

VEC=pbmc@reductions$umap@cell.embeddings
N=20
set.seed(123)
K=kmeans(VEC,centers=N)
CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)
pdf('4_.pdf')
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
dev.off() 

pdf('5.pdf')
ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
dev.off() 
used.cells <- CellSelector(plot = ppp)

pbmc@meta.data$celltype=rep(NA,length(pbmc@meta.data$batch))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
pdf('6.pdf')
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
dev.off() 
all.markers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
all.markers <- FindAllMarkers(object = pbmc_batch, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top2 <- all.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)%>%print(n = 50)

HB.genes <- HB.genes[!is.na(HB.genes)]
pbmc_batch1[["percent.HB"]]<-PercentageFeatureSet(pbmc_batch1,features=HB.genes)

# write  nGene_nUMI_mito_HB
head(pbmc_batch1@meta.data)[,c(2,3,4,5)]
VlnPlot(pbmc_batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 4)#+scale_color_npg() 
pdf('7_Feature.pdf',width=16,height=16)
VlnPlot(pbmc_batch1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), ncol = 2)
dev.off()

library(ggsci)
pdf('8_Feature.pdf',width=16,height=8)
plot1 <- FeatureScatter(pbmc_batch1, feature1 = "nCount_RNA", feature2 = "percent.mt")+scale_color_npg()
plot2 <- FeatureScatter(pbmc_batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+scale_color_npg()
plot3 <- FeatureScatter(pbmc_batch1, feature1 = "nCount_RNA", feature2 = "percent.HB")+scale_color_npg()
CombinePlots(plots = list(plot1, plot2,plot3),legend="none")
dev.off()


pdf('9_dimheatmap.pdf', width=8,height=16)
DimHeatmap(pbmc_batch1, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf('10_.pdf',width=16,height=8)
pbmc_batch2 <- JackStraw(pbmc_batch1, num.replicate = 100)
pbmc_batch3 <- ScoreJackStraw(pbmc_batch2, dims = 1:20)
plot1<-JackStrawPlot(pbmc_batch3, dims = 1:15)
plot2<-ElbowPlot(pbmc_batch3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

all.markers %>% group_by(orig.ident) %>% top_n(n = 2, wt = avg_logFC)

pbmc_batch1 <- subset(pbmc_batch1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc_batch1 <- NormalizeData(pbmc_batch1, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_batch1 <- FindVariableFeatures(pbmc_batch1, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc_batch1), 10)
# plot variable features with and without labels
pdf('12_highly_variable_genes.pdf', width=16,height=8)
plot1 <- VariableFeaturePlot(pbmc_batch1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()


pdf('13_feature_plot.pdf', width=16,height=8)
FeaturePlot(object = pbmc_batch1, features = c("IGHG3" "IGHG4" "TFPI2" "IGLC3" "IGHG1" "PSG1"))
dev.off()

pdf('13_featureplot_top2.pdf', width=16,height=16)
FeaturePlot(object = pbmc_batch1, features = c("CCL5", "CD2", "IGHG4", "IGHG3", "CP", "IFI44L", "MMRN1","MYH11", "GNLY", "NR4A2","SPARCL1", "VCAN", "HBB", "RNA28S5", "RNASE1", "C1QC", "HBG1","HBG2", "IGHG4", "IGHG3","IGFBP3", "IL1RL1","HLA-G", "PRG2", "SPP1","PEG10", "IGF2", "TNFSF10", "PEG10", "MT-ATP6"))
dev.off()



pbmc_batch1@meta.data$celltype=rep(NA,length(pbmc_batch1@meta.data$batch))
pbmc_batch1@meta.data$celltype[which(pbmc_batch1@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('14_DoHeatmap.pdf', width=16,height=8)
DoHeatmap(object = pbmc_batch1, features = top10$gene) + NoLegend()
dev.off()

subpbmc_deci1<-subset(x = pbmc_batch1,idents="deci1")
subpbmc_deci2<-subset(x = pbmc_batch1,idents="deci2")
subpbmc_deci3<-subset(x = pbmc_batch1,idents="deci3")
subpbmc_deci4<-subset(x = pbmc_batch1,idents="deci4")
subpbmc_deci5<-subset(x = pbmc_batch1,idents="deci5")
subpbmc_deci6<-subset(x = pbmc_batch1,idents="deci6")
subpbmc_deci7<-subset(x = pbmc_batch1,idents="deci7")
subpbmc_deci8<-subset(x = pbmc_batch1,idents="deci8")
subpbmc_deci9<-subset(x = pbmc_batch1,idents="deci8")
subpbmc_pla1<-subset(x = pbmc_batch1,idents="pla1")
subpbmc_pla2<-subset(x = pbmc_batch1,idents="pla2")
subpbmc_pla3<-subset(x = pbmc_batch1,idents="pla3")
subpbmc_pla4<-subset(x = pbmc_batch1,idents="pla4")
subpbmc_pla5<-subset(x = pbmc_batch1,idents="pla5")
subpbmc_pla6<-subset(x = pbmc_batch1,idents="pla6")
subpbmc_pla7<-subset(x = pbmc_batch1,idents="pla7")

DefaultAssay(subpbmc_deci1) <- "integrated"

# Run the standard workflow for visualization and clustering
subpbmc_deci1 <- ScaleData(subpbmc_deci1, verbose = FALSE)
subpbmc_deci1 <- RunPCA(subpbmc_deci1, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
subpbmc_deci1 <- RunUMAP(subpbmc_deci1, reduction = "pca", dims = 1:20)
subpbmc_deci1 <- FindNeighbors(subpbmc_deci1, reduction = "pca", dims = 1:20)
subpbmc_deci1 <- FindClusters(subpbmc_deci1, resolution = 0.5)
# Visualization
pdf('17_DimplotDECI1.pdf', width=16,height=8)
p1 <- DimPlot(subpbmc_deci1, reduction = "umap")
p2 <- DimPlot(subpbmc_deci1, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

deci1.markers <- FindMarkers(object = subpbmc_deci1, min.pct = 0.25)
head(x = deci1.markers, n = 5)

VEC_subpbmc_deci1=subpbmc_deci1@reductions$umap@cell.embeddings
N=20
set.seed(123)
K=kmeans(VEC,centers=N)
CLUST=K$cluster
#subpbmc_deci1@meta.data$clust=as.character(CLUST)
pdf('4_deci1.pdf')
DimPlot(subpbmc_deci1, reduction.use='umap', pt.size=0.5,label=TRUE)
dev.off() 

VEC=pbmc_batch1@reductions$umap@cell.embeddings
N=20
set.seed(123)
K=kmeans(VEC,centers=N)
CLUST=K$cluster
pbmc_batch1@meta.data$clust=as.character(CLUST)
pdf('4_splitby.pdf')
DimPlot(pbmc_batch1, reduction.use='umap', group.by='clust', split.by='clust', pt.size=0.5,label=TRUE)
dev.off() 


library(Seurat)
library(cowplot)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('/lustre/fyh/lab_other_work/singlecell/output/workspace')
library(reticulate)
#use_python("/lustre/tianlab/tools/anaconda3/bin/python3.6m")
.set_python("/lustre/tianlab/tools/anaconda3/bin/python3.6m")
py_config()

library(reticulate)
source("/lustre/fyh/lab_other_work/singlecell/output/workspace/zf_hfz_test_worksp/BEER/BEER.R")


# 2019.8.5

pbmc_enhence <- mybeer$seurat
PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
pbmc_enhence=BEER.combat(pbmc_enhence) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=20, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('NB20_enhence.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()

pbmc_enhence <- mybeer$seurat
PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=20, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('NB20_enhence_nocombat.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()



pbmc_enhence <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc_enhence) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=3, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
 
saveRDS(pbmc_enhence, file='seurat.enh.RDS')





PCUSE=mybeer$select
#PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))




PCUSE=mybeer$select
#PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
pdf('PLOT.pdf',width=7,height=7)
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
dev.off()






pbmc_enhence <- mybeer$seurat
#PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
PCUSE=mybeer$select   
pbmc_enhence=BEER.combat(pbmc_enhence) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=5, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('test_enhence.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()


pbmc_enhence <- mybeer$seurat
#PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
PCUSE=mybeer$select   
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=3, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('test_enhence_nocombat.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()


pbmc_enhence <- mybeer$seurat
#PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
PCUSE=mybeer$select   
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=3, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('featureplot_nocombat.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()



pbmc_enhence <- mybeer$seurat
#PCUSE=c(1:ncol(pbmc_enhence@reductions$pca@cell.embeddings))
PCUSE=mybeer$select   
#pbmc_enhence=BEER.combat(pbmc_enhence) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc_enhence, PCUSE, NB=3, NT=10)
pbmc_enhence@reductions$umap@cell.embeddings=umap
pdf('NB3_enhence.pdf')
DimPlot(pbmc_enhence, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()

# pla-CTBs
pdf('marker_pla_CTBs.pdf')
FeaturePlot(pbmc, features=c('CDH1','EGFR','GATA3','RRM2','CDK1','PAGE4','SMAGP','EFEMP1','ERVFRD-1','INSL4','TBX3'))
dev.off()
# pla-STB
pdf('marker_pla_STB.pdf')
FeaturePlot(pbmc, features=c('CGB','HLA-G', 'CSH1', 'CSH2','HSD3B1', 'HSD17B4', 'HSD11B2', 'HSD17B12', 'HSD17B1'))
dev.off()

# pla-EVTs
pdf('marker_pla_EVTS.pdf', width=8,height=16)
FeaturePlot(pbmc, features=c('HLA-G','MMP2','HLA-G','MMP2','ASCL2','BZW2','C12orf75','EFNA1','FSTL3','KRT19','KRT8','MFAP5','PGF','QSOX1','SLC16A3'))
dev.off()
# pla-MSCs
pdf('marker_pla_MSCs.pdf')
FeaturePlot(pbmc, features=c('5B5','CD105','CD73','CD90','DLK1','COL1A1','COL1A2','HGF','SERPINH1','WNT2','MEST','MEG3','PLAC9','COL6A1'))
dev.off()
# pla-mac
pdf('marker_pla_mac.pdf')
FeaturePlot(pbmc, features=c('CD68', 'HLA-DR', 'HAM56'))
dev.off()
# pla-Endo
pdf('marker_pla_Endo.pdf', width=8,height=16)
FeaturePlot(pbmc, features=c('CD34' , 'VEGFR-2' ,'RAMP2' , 'ESAM' ,'VAMP5' ,  'PRCP','JAM2' ,'CD93' , 'VIM' , 'C8orf4','MGST2' ))
dev.off()
# deci-dsc
pdf('marker_deci_dsc.pdf')
FeaturePlot(pbmc, features=c('Prl8a2', 'Cryab', 'CD90',  'APOD', 'DKK1', 'C1R','CPXM1', 'HSPB1', 'HSPB6','MYL9', 'SERPING1', 'TIMP3','DCN', 'RBP1')
dev.off()
# deci-epi
pdf('marker_deci_epi.pdf')
FeaturePlot(pbmc, features=c('ANXA4', 'CD24', 'GPX3', 'TM4SF1'))
dev.off()
# deci-peri
pdf('marker_deci_peri.pdf')
FeaturePlot(pbmc, features=c('COL18A1', 'COL3A1', 'C11orf96', 'MGP', 'MYL9', 'CALD1'))
dev.off()
# deci-endo
pdf('marker_deci_endo.pdf')
FeaturePlot(pbmc, features=c('CD34', 'CD105', 'CD144','APP', 'CRIP2', 'CLDN5','GNG11', 'SEPW1', 'PCAT19','IGFBP7', 'EGFL7'))
dev.off()
# deci-NK
pdf('marker_deci_NK.pdf')
FeaturePlot(pbmc, features=c('CD56', 'CCL4', 'EVL', 'PRF1'))
dev.off()
# deci-M
pdf('marker_deci_M.pdf')
FeaturePlot(pbmc, features=c('CD11', 'CD14', 'CD68', 'HLA-DRA', 'CD74', 'APOE4', 'LYVE1'))
dev.off()
# deci-DC
pdf('marker_deci_DC.pdf')
FeaturePlot(pbmc, features=c('HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'CPVL', 'HLA-DMA', 'HLA-DMB', 'HLA-DPB1', 'HLA-LSP1', 'CST3', 'LYZ', 'SPl1', 'HLA-DMA'))
dev.off()
# deci-T
pdf('marker_deci_T.pdf')
FeaturePlot(pbmc, features=c('CCL5', 'CD3D', 'CD3E', 'CD3G', 'B2M', 'CD52', 'CXCR4'))
dev.off()

pdf('marker_FN1.pdf')
FeaturePlot(pbmc, features=c('FN1'))
dev.off()


##########**************follow zf******************#
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf('~/Downloads/HFZ_QC.pdf',width=12,height=7)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
> pbmc<-subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
> dim(pbmc)
[1] 13828 36797


pbmc.rename<-c('1'="NA", '2`="dS1", `3`="dNK1", `4`="EVT1", `5`="dS2", `6`="dS3",`7`="dM1", `8`="EVT2",
 `9`="dM2", `10`="dS4", `11`="dS5",`12`="dS6", `13`="dNK", `14`="dS7",`15`="VCT1", `16`="dS8",`17`="EVT3", `18`="EVT4", '19'="ENDO1",
 `20`="Fib1", '21'="VCT2",'22'="Fib2", '23'="Fib4",'24'="dS9", '25'="ElEP", '26'="Fib5",'27'="NK",
 `28`="Fib6",'29'="ds10",'30'="Fib7", '31'="CD4 T", '32'="Fib8",'33'="DC1",'34'="Fib9",'35'="DC2",'36'="DC3",'37'="DC4",'38'="PSc",'39'="Fib10",
 `40`="Ductual",'41'="VCT3",'42'="EVT5",'43'="EVT6",'44'="T", '45'="dS10",'46'="dNK2",'47'="ENDO2", '48'="Fib8",'49'="NA",'50'="NA")

1	Fib1
2	dS1
3	dNK1
4	EVT1
5	dS2
6	VCT1
7	dM1
8	EVT2
9	dM2
10	dS4
11	dS5
12	dS6
13	dNK
14	dS7
15	VCT2
16	dS8
17	EVT3
18	EVT4
19	ENDO1
20	EVT5
21	VCT3
22	Fib2
23	Fib3
24	Fib4
25	ElEP
26	Fib5
27	NK
28	Fib6
29	Fib7
30	Fib8
31	T
32	VCT4
33	EVT
34	Fib9
35	DC
36	VCT5
37	VCT6
38	PV
39	Fib9
40	Ductual
41	VCT7
42	EVT5
43	EVT6
44	NKT
45	dS10
46	dNK2
47	ENDO2
48	Fib10
49	Fib11
50	ElEp

# dM:decidual macrophage,DC:dendrtic cell，PV:perivascular cell,ElEp:Erythroid-like and erythroid precursor cells
current_cluster_ids <- c( 1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,
	35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)
new_cluster_ids <- c("1-Fib1", "2-dS1", "3-dNK1", "4-EVT1", "5-dS2", "6-VCT1","7-dM1","8-EVT2","9-dM2", "10-dS4", "11-dS5","12-dS6","13-dNK",
 "14-dS7","15-VCT2", "16-dS8","17-EVT3","18-EVT4", "19-ENDO1",
 "20-EVT5", "21-VCT3","22-Fib2", "23-Fib3","24-Fib4", "25-ElEp", "26-Fib5","27-NK","28-Fib6","29-Fib7","30-Fib8", "31-T","32-VCT4",
 "33-EVT","34-Fib9","35-DC","36-VCT5","37-VCT6","38-PV","39-Fib9"
 "40-Ductual","41-VCT7","42-EVT5","43-EVT6","44-NKT", "45-dS10","46-dNK2","47-ENDO2", "48-Fib10","49-Fib11","50-ElEp")

# Changing IDs to cell type
pbmc@meta.data$clust <- plyr::mapvalues(x = pbmc@meta.data$clust, 
                                from = current_cluster_ids, 
                                to = new_cluster_ids)
