library(Seurat)
library(Matrix)
library(zinbwave)

expr <- Read10X('../../data/b_cells/')
pbmc <- CreateSeuratObject(expr)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method='vst', nfeatures=1000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pbmc@meta.data$logsf <- log(pbmc@meta.data$nCount_RNA/mean(pbmc@meta.data$nCount_RNA))
sceobj <- as.SingleCellExperiment(pbmc)
sceobj_vf <- sceobj[VariableFeatures(pbmc),]
sceobj_vf <- zinbwave(sceobj_vf, K=10, X='~G2M.Score+S.Score+logsf')

write.csv(reducedDim(sceobj_vf), 'zinbwave.csv')


