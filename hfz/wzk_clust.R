dec_DC_M <- c(55,5,99,200,162,153,30,171,143)
dec_ds <- c()
Endo <- c(195,23,91,101,61,57,17,25,52,103)
dec_NK <- c(93,12,82,50,167,163,80)
dec_T <- c(186,48,110,2,154,84,179,199,152,150,165,29,51,190,175,178,79,169,54,194,39,45)
dec_T_NK <- c(90,157)
dec_Epi <- c()
dec_PV <- c(76,129,125,135,189,60,177)
pla_CTB <- c(192,87,34,133,123,27,7,71,68,58,92,37,120,130,33,151,43,182)
pla_STB <- c(89,132,106,18)
pla_EVT <- c(145,141,187,1,193,97,94,149,161,19,188,56,47,147,116,155)
pla_MSC <- c()


for (i in 1:200){
    if (i %in% dec_DC_M){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'DC_M'
    }
    else if (i %in% dec_ds){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'dec_ds'
    }
    else if (i %in% Endo){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'Endo'
    }
    else if (i %in% dec_NK){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'NK'
    }
    else if (i %in% dec_T){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'T'
    }
    else if (i %in% dec_T_NK){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'T_NK'
    }
    else if (i %in% dec_Epi){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'dec_Epi'
    }
    else if (i %in% dec_PV){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'dec_PV'
    }
    else if (i %in% pla_CTB){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'CTB'
    }
    else if (i %in% pla_STB){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'STB'
    }
    else if (i %in% pla_EVT){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'EVT'
    }
    else if (i %in% pla_MSC){
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'MSC'
    }
    else{
        pbmc@meta.data$annoClust[pbmc@meta.data$clust == as.character(i)] <- 'UnKnow'
    }
}

DimPlot(pbmc, reduction = 'umap', group.by = 'annoClust', pt.size = 0.3, label = TRUE)
