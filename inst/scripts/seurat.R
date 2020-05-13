library(Seurat)
library(dplyr)

# beginning from 'sce' defined in vignette

# create seurat object for whitelisted cells
cts <- assays(sce)[["counts"]]
pbmc <- CreateSeuratObject(cts)
pbmc

mt.genes <- rownames(sce)[as.logical(seqnames(sce) == "chrM")]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt.genes)
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

# Assigning Cluster names based on marker genes
# "CD14+ Monocytes" "CST3" "ENSG00000101439.9"
# "CD4 T" "IL7R" "ENSG00000168685.15"
# "B" "MS4A1" "ENSG00000156738.17"
# "CD8 T" "KLRB1" "ENSG00000111796.4" (we use only for this particular set of cells)
# "NK" "GNLY" "ENSG00000115523.16"

markers.vec <- c("CD14+ Monocytes"="CST3", "CD4 T"="IL7R", "B"="MS4A1", "CD8 T"="KLRB1", "NK"="GNLY")
clusters <- character(5)
for (i in seq_along(markers.vec)) {
  idx <- which(top10$symbol == markers.vec[i])
  stopifnot(length(idx) == 1)
  cluster.idx <- as.numeric(as.character(top10$cluster[idx])) + 1
  clusters[cluster.idx] <- names(markers.vec)[i]
}

top10$cluster.name <- rep(clusters, each=10)
df <- as.data.frame(top10)
write.csv(df, file="../extdata/top10.csv")

names(clusters) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, clusters)
saveRDS(Idents(pbmc), file="../extdata/idents.rds")

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 
