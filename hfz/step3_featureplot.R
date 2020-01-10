library(Seurat)
library(cowplot)
source('/home/wzk/tools/BEER-master/BEER.R')
setwd("/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata")
#.set_python(‘/usr/bin/python’)
.set_python('/usr/bin/python')

print(getwd())

pbmc=readRDS('pbmc.enh.RDS')


tiff("./Marker/dSC2.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('DKK1','IGFBP1','PRL','CYP11A1','TAGLN','IGFBP2','IGFBP6','ACTA2'),ncol=3,pt.size=0.01)
dev.off()
tiff("./Marker/NK.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('NKG7','NCAM1','GZMA','GZMB','MKI67','CS1','TOP2A','PRF1','XCL1','XCL2'),ncol=3,pt.size=0.01)
dev.off()
tiff("./Marker/T.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CCL5','CD3D','CD69','CXCR4','CD3G','LTB','IL7R','GZMA','GZMK'),ncol=3,pt.size=0.01)
dev.off()
tiff("./Marker/DC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD74','CD14','CLEC9A','CD52','CD1C','CD83','CD86'),ncol=3,pt.size=0.01)
dev.off()


tiff("./Marker/CTB.cytotrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('KRT7','PEG10','PARP1','PAGE4'),ncol=2,pt.size=0.01)
dev.off()
tiff("./Marker/STB.syncytiotrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CGA','HSD3B1','CYP19A1','ERVW-1'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/EVT.ExtravillousTrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('HLA-G','DIO2','LAIR2','HSD3B1'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/fib1_fib2.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('COL1A1','COL1A2','COL3A1','ACTA2'),ncol=2,pt.size=0.01)
dev.off()
tiff("./Marker/VEC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD34','PLVAP','CDH5','CD105','CD144','PCDH17','ICAM1'),ncol=3,pt.size=0.01)
dev.off()
tiff("./Marker/EEC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('PAEP','SLC18A2','PAEP'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/LEC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('LYVE1','STMN1','CD9','FABP5'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/M.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD74','MS4A4A','CD14','STAB1','SEPP1','MS4A7','CD74'),ncol=3,pt.size=0.01)
dev.off()

tiff("./Marker/fetal.endothelial.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD34','CDH5','PLVAP','CD133'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/HB.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD68','AIF1','CD14','CD68'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/SMC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('PI15','NDUFA4L2','MYH11','ACTA2','RGS5'),ncol=2,pt.size=0.01)
dev.off()
tiff("./Marker/pFB.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('COL1A1','COL3A1','COL1A2','ACTA2','DLK1'),ncol=2,pt.size=0.01)
dev.off()













