library(Seurat)
library(Matrix)
library(zinbwave)
library(dplyr)

exp_data <- readMM('../../data/macparland/matrix.mtx')
exp_data <- apply(exp_data, c(1,2), function(x){
	as.integer(x)
})
exp_data_meta <- read.csv('../../data/macparland/meta.csv', row.names=1)
colnames(exp_data) <- rownames(exp_data_meta)
rownames(exp_data) <- read.csv('../../data/macparland/features.tsv')$x
size_factor <- exp_data_meta$nCount_RNA/mean(exp_data_meta$nCount_RNA)
pbmc <- CreateSeuratObject(exp_data)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method='vst', nfeatures=2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

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
sceobj <- zinbsurf(sceobj, K=20, X='~factor(orig.ident)+logsf')

write.csv(reducedDim(sceobj), 'zinbwave.csv')

