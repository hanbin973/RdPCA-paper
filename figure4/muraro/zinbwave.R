library(Seurat)
library(Matrix)
library(zinbwave)
library(dplyr)

exp_data <- read.csv('../../data/muraro/matrix.csv', row.names=1)
exp_data_meta <- read.csv('../../data/muraro/meta.csv', row.names=1)[1:2285,]
colnames(exp_data) <- rownames(exp_data_meta)
size_factor <- exp_data_meta$nCount_RNA/mean(exp_data_meta$nCount_RNA)
pbmc <- CreateSeuratObject(exp_data)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method='vst', nfeatures=2000)

donor <- read.csv('/data01/hanbin973/batch/batch_labels.csv')$X0
pbmc@meta.data$donor <- donor[1:2285]
pbmc@meta.data$logsf <- log(size_factor)

pbmc@meta.data$logsf <- log(pbmc@meta.data$nCount_RNA/mean(pbmc@meta.data$nCount_RNA))
sceobj <- as.SingleCellExperiment(pbmc)
filter1 <- rowSums(assay(sceobj)>5)>5
sceobj <- sceobj[filter1,]
assay(sceobj) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sceobj)
vars <- sort(vars, decreasing = TRUE)
sceobj <- sceobj[names(vars)[1:2000],]
assay(sceobj) <- apply(assay(sceobj), c(1,2), function(x){
        as.integer(x)
})
sceobj <- zinbsurf(sceobj, K=30, X='~factor(donor)+logsf')
write.csv(reducedDim(sceobj), 'zinbwave.csv')

