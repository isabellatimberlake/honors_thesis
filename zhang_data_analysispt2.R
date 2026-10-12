# load packages and set up environment
library(Seurat)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(data.table)
library(stringr)
library(patchwork) # for plotting
library(presto) # for differential expression


setwd("/gpfs1/pi/kmccreer/scRNAseq/issy_26-27/honors_thesis")

## Visualizing data! 
obj <- readRDS("merged_seurat_obj_allRNAsamples.rds")

# fixing sample name column
obj$sampleName <- sub("-Cell-[0-9]+$", "", colnames(obj))

obj@meta.data <- obj@meta.data |>
  mutate(sampleName = str_remove_all(sampleName, "_Cell_[0-9]+$"))

# plotting reads per sample
obj@meta.data |>
  ggplot(aes(x = sampleName, y = nCount_RNA, fill = sampleName)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(fill = "none") +
  labs(x = NULL, y = "UMI counts per cell (log10)", title = "Reads per Cell by Sample")

# find marker genes for each group
FindMarkers(
  object,
  slot = "data",
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  fc.results = NULL,
  densify = FALSE
)

# PCA plot
# scale data first
# log normalizes
objN <- NormalizeData(obj)
objNFVM <- FindVariableFeatures(objN)
obj_scale <- ScaleData(objNFVM)

obj_PCA <- RunPCA(obj_scale, verbose = FALSE,assay="RNA")

## investigating the PCAs
# elbow plot
ElbowPlot(obj_PCA, ndims = 40)
# appears to taper off at about 20-30
# now, calculate clusters
obj_PCA <- FindNeighbors(obj_PCA, dims = 1:30)
obj_PCA <- FindClusters(obj_PCA, resolution = 0.9)

# visualize
obj_PCA <- RunUMAP(obj_PCA,dims = 1:30)
# look at clusters vs sample type
DimPlot(obj_PCA, reduction = "umap", group.by = "sampleName",
        alpha=0.4, ncol=2)


# pseudoBULK ####

# check this first:
n_distinct(obj$sampleName)

# Summed raw counts per sample (pseudobulk, e.g. for DESeq2)
pb <- AggregateExpression(obj, group.by = "sampleName", normalization.method = "LogNormalize")$RNA



