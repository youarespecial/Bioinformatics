pbmc=readRDS('pbmc1106.RDS')


VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=200
set.seed(1)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

tiff("ID.tiff", width = 8, height= 8, units = 'in',res = 400)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)+NoLegend()
dev.off()

################################################################


pbmc@meta.data$decidua=rep(0,ncol(pbmc))
pbmc@meta.data$decidua[which(pbmc@meta.data$batch %in% c('deci1-NM',
               'deci2-NM','deci7-NM','deci8-PE','deci9-PE',                                      
                'deci3-NM','deci6-PE','deci4-PE','deci5-NM'         
                                                     ))]=1
tiff("decidua.tiff", width = 8, height= 8, units = 'in',res = 400)
FeaturePlot(pbmc, features='decidua',order=TRUE)
dev.off()


pbmc@meta.data$placenta=rep(1,ncol(pbmc))
#pbmc@meta.data$placenta[which(pbmc@meta.data$batch %in% c('pla1-NM',
               # 'pla2-NM','pla3-NM','pla4-PE','pla5-PE',                                      
               #  'pla6-PE','pla7-PE'         
               #                                       ))]=0
pbmc@meta.data$placenta[which(pbmc@meta.data$batch %in% c('deci1-NM',
               'deci2-NM','deci7-NM','deci8-PE','deci9-PE',                                      
                'deci3-NM','deci6-PE','deci4-PE','deci5-NM'
                                                     ))]=0
tiff("placenta.tiff", width = 8, height= 8, units = 'in',res = 400)
FeaturePlot(pbmc, features='placenta',order=TRUE)
dev.off()
###########################     step4_giveCellName.R      ############

###########################     step4_giveCellName.R      ############


tiff("Cluster.tiff", width = 8, height= 8, units = 'in',res = 400)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
dev.off()

saveRDS(object = pbmc@meta.data,file='META.RDS')

####################



pbmc@meta.data$tag=rep('placenta',ncol(pbmc))
pbmc@meta.data$tag[which(pbmc@meta.data$batch %in% c('deci1-NM',
               'deci2-NM','deci7-NM','deci8-PE','deci9-PE',                                      
                'deci3-NM','deci6-PE','deci4-PE','deci5-NM'          
                                                     ))]='decidua'

TAB=table(pbmc@meta.data$clust, pbmc@meta.data$tag)


TAB=table(pbmc@meta.data$clust, pbmc@meta.data$batch)

.writeTable(t(TAB),PATH = 'TABLE.txt')
