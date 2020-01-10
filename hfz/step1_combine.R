library(Seurat)
library(cowplot)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
.set_python('/usr/bin/python')
setwd('/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/deci_NM')
#  rm(list = ls())

# *******************************deci NM************************************** #
D1=read.table(gzfile("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/decidua0117_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D2=read.table(gzfile("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/decidua0417-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D3=read.table(gzfile("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/decidua508_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D6=read.table(gzfile("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/decidua-2018_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D9=read.table(gzfile("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/Decidua20190420_L3_1000903_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
colnames(D1)=paste0('deci1','_',colnames(D1))
colnames(D2)=paste0('deci2','_',colnames(D2))
colnames(D3)=paste0('deci7','_',colnames(D3))
colnames(D6)=paste0('deci3','_',colnames(D6))
colnames(D9)=paste0('deci5','_',colnames(D9))
DD1=.simple_combine(D1,D2)$combine
DD2=.simple_combine(D3,D6)$combine
DD3=.simple_combine(DD1,DD2)$combine
DATA=.simple_combine(DD3,D9)$combine
# *******************************deci NM************************************** #


saveRDS(DATA, file='./DATA.RDS')
.get_batch<-function(x){
        y=unlist(strsplit(x,'_'))[1]
        return(y)
        } 

BATCH=apply(as.matrix(colnames(DATA)),1,.get_batch)
saveRDS(BATCH, file='./BATCH.RDS')

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')
##################################

#setwd('/users/zha8dh/tianlab/HFZ')
# source('./BEER.R')
mybeer=readRDS('./mybeer.RDS')


PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'

#############
pbmc <- mybeer$seurat
dim(mybeer$seurat)
# 16087 20000
################################################
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
dim(pbmc)
# 16087 12671
saveRDS(pbmc, file='pbmc.RDS')
