library(Seurat)
library(Matrix)
library(zinbwave)
library(dplyr)

expr <- readMM('../../data/pbmc_10k_v3/matrix.mtx')
rownames(expr) <- read.csv('../../data/pbmc_10k_v3/features.tsv')$x
colnames(expr) <- read.csv('../../data/pbmc_10k_v3/barcodes.tsv')$x
pbmc <- CreateSeuratObject(expr)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method='vst', nfeatures=3000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

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
sceobj <- zinbsurf(sceobj, K=20, X='~G2M.Score+S.Score+logsf')
write.csv(reducedDim(sceobj), 'zinbwave.csv')

