# 将所有结果文件移动到一起
mv node9_fq/*.sprint nat_dev_organ_sprint/
# 将res文件提出，以保证原始文件不动
library(Seurat)
library(cowplot)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('/lustre/fyh/lab_other_work/singlecell/output/workspace')


deci1.data <- read.table(file = "1_decidua0117_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla1.data <- read.table(file = "1_placenta0423-2_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci2.data <- read.table(file = "2_decidua0417-2_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla2.data <- read.table(file = "2_placenta-2018_combined_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci3.data <- read.table(file = "3_decidua-2018_combined_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla3.data <- read.table(file = "3_Placenta20190402_L3_1000902_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci4.data <- read.table(file = "4_decidua20190215_combined_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla4.data <- read.table(file = "4_placenta-2019_combined_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci5.data <- read.table(file = "5_Decidua20190420_L3_1000903_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla5.data <- read.table(file = "5_placenta508-1_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci6.data <- read.table(file = "6_decidua-2019_combined_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla6.data <- read.table(file = "6_placenta508-2_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci7.data <- read.table(file = "7_decidua508_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
pla7.data <- read.table(file = "7_placenta514-2_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci8.data <- read.table(file = "8_decidua510_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)
deci9.data <- read.table(file = "9_decidua514-2_exon_tagged.dge.txt.gz", sep = "\t", header = T, row.names= 1)

######????????????????????????????????????????????????????????##################################
metadata <- read.table(file = "metadata.txt", sep = "\t", header = T)
######????????????????????????????????????????????????????????##################################

deci1 <- CreateSeuratObject(counts = deci1.data, project = "DECI-NM", min.cells = 3, min.features = 200, meta.data=metadata)
deci1 <- subset(deci1, subset = nFeature_RNA > 500)
deci1 <- NormalizeData(deci1, verbose = FALSE)
deci1 <- FindVariableFeatures(deci1, selection.method = "vst", nfeatures = 2000,)

deci2 <- CreateSeuratObject(counts = deci2.data, project = "DECI-NM", min.cells = 5,meta.data=metadata)
deci2 <- subset(deci2, subset = nFeature_RNA > 500)
deci2 <- NormalizeData(deci2, verbose = FALSE)
deci2 <- FindVariableFeatures(deci2, selection.method = "vst", nfeatures = 2000)

deci3 <- CreateSeuratObject(counts = deci3.data, project = "DECI-NM", min.cells = 5,meta.data=metadata)
deci3 <- subset(deci3, subset = nFeature_RNA > 500)
deci3 <- NormalizeData(deci3, verbose = FALSE)
deci3 <- FindVariableFeatures(deci3, selection.method = "vst", nfeatures = 2000)

deci4 <- CreateSeuratObject(counts = deci4.data, project = "DECI-PE", min.cells = 5,meta.data=metadata)
deci4 <- subset(deci4, subset = nFeature_RNA > 500)
deci4 <- NormalizeData(deci4, verbose = FALSE)
deci4 <- FindVariableFeatures(deci4, selection.method = "vst", nfeatures = 2000)

deci5 <- CreateSeuratObject(counts = deci5.data, project = "DECI-NM", min.cells = 5,meta.data=metadata)
deci5 <- subset(deci5, subset = nFeature_RNA > 500)
deci5 <- NormalizeData(deci5, verbose = FALSE)
deci5 <- FindVariableFeatures(deci5, selection.method = "vst", nfeatures = 2000)

deci6 <- CreateSeuratObject(counts = deci6.data, project = "DECI-PE", min.cells = 5,meta.data=metadata)
deci6 <- subset(deci6, subset = nFeature_RNA > 500)
deci6 <- NormalizeData(deci6, verbose = FALSE)
deci6 <- FindVariableFeatures(deci6, selection.method = "vst", nfeatures = 2000)

deci7new <- CreateSeuratObject(counts = deci7.data, project = "DECI-PE", min.cells = 5,meta.data=metadata,min.features=200)
deci7 <- subset(deci7, subset = nFeature_RNA > 500)
deci7 <- NormalizeData(deci7, verbose = FALSE)
deci7 <- FindVariableFeatures(deci7, selection.method = "vst", nfeatures = 2000)

deci8 <- CreateSeuratObject(counts = deci8.data, project = "DECI-PE", min.cells = 5,meta.data=metadata)
deci8 <- subset(deci8, subset = nFeature_RNA > 500)
deci8 <- NormalizeData(deci8, verbose = FALSE)
deci8 <- FindVariableFeatures(deci8, selection.method = "vst", nfeatures = 2000)

deci9 <- CreateSeuratObject(counts = deci9.data, project = "DECI-PE", min.cells = 5,meta.data=metadata)
deci9 <- subset(deci9, subset = nFeature_RNA > 500)
deci9 <- NormalizeData(deci9, verbose = FALSE)
deci9 <- FindVariableFeatures(deci9, selection.method = "vst", nfeatures = 2000)


pla1 <- CreateSeuratObject(counts = pla1.data, project = "PLA-NM", min.cells = 5,meta.data=metadata)
pla1 <- subset(pla1, subset = nFeature_RNA > 500)
pla1 <- NormalizeData(pla1, verbose = FALSE)
pla1 <- FindVariableFeatures(pla1, selection.method = "vst", nfeatures = 2000)

pla2 <- CreateSeuratObject(counts = pla2.data, project = "PLA-NM", min.cells = 5,meta.data=metadata)
pla2 <- subset(pla2, subset = nFeature_RNA > 500)
pla2 <- NormalizeData(pla2, verbose = FALSE)
pla2 <- FindVariableFeatures(pla2, selection.method = "vst", nfeatures = 2000)

pla3 <- CreateSeuratObject(counts = pla3.data, project = "PLA-NM", min.cells = 5,meta.data=metadata)
pla3 <- subset(pla3, subset = nFeature_RNA > 500)
pla3 <- NormalizeData(pla3, verbose = FALSE)
pla3 <- FindVariableFeatures(pla3, selection.method = "vst", nfeatures = 2000)

pla4 <- CreateSeuratObject(counts = pla4.data, project = "PLA-PE", min.cells = 5,meta.data=metadata)
pla4 <- subset(pla4, subset = nFeature_RNA > 500)
pla4 <- NormalizeData(pla4, verbose = FALSE)
pla4 <- FindVariableFeatures(pla4, selection.method = "vst", nfeatures = 2000)

pla5 <- CreateSeuratObject(counts = pla5.data, project = "PLA-PE", min.cells = 5,meta.data=metadata)
pla5 <- subset(pla5, subset = nFeature_RNA > 500)
pla5 <- NormalizeData(pla5, verbose = FALSE)
pla5 <- FindVariableFeatures(pla5, selection.method = "vst", nfeatures = 2000)

pla6 <- CreateSeuratObject(counts = pla6.data, project = "PLA-PE", min.cells = 5,meta.data=metadata)
pla6 <- subset(pla6, subset = nFeature_RNA > 500)
pla6 <- NormalizeData(pla6, verbose = FALSE)
pla6 <- FindVariableFeatures(pla6, selection.method = "vst", nfeatures = 2000)

pla7 <- CreateSeuratObject(counts = pla7.data, project = "PLA-NM", min.cells = 5,meta.data=metadata)
pla7 <- subset(pla7, subset = nFeature_RNA > 500)
pla7 <- NormalizeData(pla7, verbose = FALSE)
pla7 <- FindVariableFeatures(pla7, selection.method = "vst", nfeatures = 2000)


immune.anchors <- FindIntegrationAnchors(object.list = list(deci1,deci2,deci3, deci4,deci5,deci6,deci7,deci8,deci9,pla1,pla2,pla3,pla4,pla5,pla6,pla7), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#############

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

pbmcbig.anchors <- FindIntegrationAnchors(object.list = pbmc.big), dims = 1:20)
pbmcbig.combined <- IntegrateData(anchorset = pbmcbig.anchors, dims = 1:20)

# Visualization
p1 <- DimPlot(pbmc.big, reduction = "umap", group.by='orig.ident')
p2 <- DimPlot(pbmc.big, reduction = "umap", label = TRUE)
pdf(1_4.pdf)
plot_grid(p1,p2)
dev.off()





mito.genes <- grep(pattern = "^MT-", x = rownames(x = immune.combined), value = TRUE)
percent.mito <- Matrix::colSums(immune.combined@counts[mito.genes, ])/Matrix::colSums(immune.combined@counts)
spleen <- AddMetaData(object = immune.combined, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = immune.combined, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)



pbmc.big <- merge(deci1, y = c(deci2,deci3, deci4,deci5,deci6,deci7,deci8,deci9,pla1,pla2,pla3,pla4,pla5,pla6,pla7), add.cell.ids = c("deci1", "deci2","deci3", "deci4","deci5","deci6","deci7","deci8","deci9","pla1","pla2","pla3","pla4","pla5","pla6","pla7"), project = "PBMC15K")




pancreas <- CreateSeuratObject(pbmc.big, meta.data = metadata)
pancreas.list <- SplitObject(pancreas, split.by = "disease")


pbmc=EXP
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@counts[mito.genes, ])/Matrix::colSums(pbmc@counts)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(1000, -Inf), high.thresholds = c(2000, 0.2))

table(pbmc@ident)






