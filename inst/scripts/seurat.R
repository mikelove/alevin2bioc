library(Seurat)
library(dplyr)

# beginning from 'sce' defined in vignette

# create seurat object for whitelisted cells
cts <- assays(sce)[["counts"]]
pbmc <- CreateSeuratObject(cts)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filtering some outlier cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20 & nCount_RNA < 20000)
pbmc

# Normalize and runPCA
pbmc <- SCTransform(pbmc) %>% RunPCA()

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:25)
pbmc <- FindClusters(pbmc, resolution = 0.2)

# plotting umap embeddings
pbmc <- RunUMAP(pbmc, dims = 1:25)
DimPlot(pbmc, label = T)

# Find marker genes & extract top 10 genes
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
top10$symbol <- mcols(sce)[top10$gene,"SYMBOL"]
top10$cluster.name <- rep(c("CD14+ Monocytes","CD4 T","B","CD8 T","NK"), each=10)
df <- as.data.frame(top10)
write.csv(df, file="top10.csv")
saveRDS(Idents(pbmc), file="idents.rds")

# Assigning Cluster names based on marker genes
# "CD14+ Monocytes" "CD14" 0 (CST3 and LYZ)
# "CD4 T" "IL7R" 1
# "B" "MS4A1" 2
# "CD8 T" "CD8A/B" 3 (CD8A comes up in top 20)
# "Natural Killer" "GNLY" 4
new.cluster.ids <- c("CD14+ Monocytes", "CD4 T", "B", "CD8 T", "NK")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 
