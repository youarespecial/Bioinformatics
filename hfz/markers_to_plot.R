# TB
FeaturePlot(pbmc, features = c('KRT7', 'ERVW-1', 'CYP19A1', 'CGA', 'HLA-G', 'DIO2', 'PAGE4', 'PEG10', 'PARP1'), 
            reduction = 'umap')
# HC
FeaturePlot(pbmc, features = c('CD163', 'CD14', 'CSF1R', 'CD68', 'CCL4', 'AIF1', 'FCGRT'), 
            reduction = 'umap')
# NK
FeaturePlot(pbmc, features = c('NKG7', 'NCAM1', 'GZMA', 'GZMB', 'MKI67', 'TOP2A', 'XCL1', 'XCL2', 'PRF1'), 
            reduction = 'umap')
# T
FeaturePlot(pbmc, features = c('CD3D', 'CD3G', 'CD69', 'LTB', 'CXCR4', 'CCL5', 'IL7R', 'GZMA', 'GZMK'), 
            reduction = 'umap')
# T reg
FeaturePlot(pbmc, features=c('CD4', 'IL2RA', 'FOXP3'), 
            reduction = 'umap')
# FB
FeaturePlot(pbmc, features = c('COL1A1', 'COL1A2', 'COL3A1', 'ACTA2', 'DLK1'), 
            reduction = 'umap')
# EB
FeaturePlot(pbmc, features = c('HBE1', 'HBG1', 'HBA1', 'HBM', 'HBG2', 'HBA2', 'HBZ', 'GYPA', 'GYPB', 'GYPC'), 
            reduction = 'umap')
# SC
# FeaturePlot(pbmc, features = c('DKK1', 'IGFBP1', 'PRL', 'APOA1', 'CHI3L2', 'SERPINA3', 'IL1B', 'PROK1'))
FeaturePlot(pbmc, features = c('DKK1', 'TAGLN', 'PRL', 'IGFBP1','IGFBP6','IGFBP2','CYP11A1'), 
            reduction = 'umap')
# VEC
FeaturePlot(pbmc, features = c('CD34', 'PLVAP', 'CDH5', 'PCDH17', 'ICAM1'), 
            reduction = 'umap')
# LEC
FeaturePlot(pbmc, features = c('LYVE1', 'STMN1', 'CD9', 'FABP5', 'PDPN', 'FLT4', 'PROX1'), 
            reduction = 'umap')
# deci FBs
# FeaturePlot(pbmc, features = c('MME', 'ITGB1', 'CD47', 'CD81', 'LRP1', 'IL1R1'))
# EEC
FeaturePlot(pbmc, features=c('PAEP', 'SLC18A2'), 
            reduction = 'umap')
# SMC
# FeaturePlot(pbmc, features = c('PI15', 'NDUFA4L2', 'MYH11', 'ACTA2', 'RGS5'))
FeaturePlot(pbmc, features = c('MCAM', 'PDGFRB', 'CNN1', 'MYH11', 'CSPG4', 'ANPEP'), 
            reduction = 'umap')
# M
FeaturePlot(pbmc, features = c('CD74', 'CD14', 'MS4A4A', 'STAB1', 'SEPP1', 'MS4A7'), 
            reduction = 'umap')
# DC adding: ID2, IRF8
FeaturePlot(pbmc, features = c('CD74', 'CD14', 'CLEC9A', 'CD52', 'CD1C', 'CD83', 'CD86', 'ID2', 'IRF8'), 
            reduction = 'umap')
