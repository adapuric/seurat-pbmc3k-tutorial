# RUN THROUGH OF scRNAseq WORKFLOW USING SEURAT 
# using 10x Genomics 3k PBMC dataset (classic, small enough to run on a laptop)
# Peripheral Blood Mononuclear Cells (PBMC) data - 2,700 single cells 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial?utm_source=chatgpt.com

# First Step = Reading the Data 
library(dplyr)
library(patchwork)
library(Seurat)

# steps in the terminal: 
# raw file is in downloads, move it into the working directory from terminal
#     mv pbmc3k_filtered_gene_bc_matrices.tar ~/intro_workflow_seurat
# Extract the archive (in your terminal, in the same folder as the tar)
#     mkdir pbmc3k
#.    tar -xvf pbmc3k_filtered_gene_bc_matrices.tar -C pbmc3k
# confirm stucture ls pbmc3k/filtered_gene_bc_matrices/hg19 --> see barcodes/genes/matrix.gsvm 


# make sure we are in the right project directory 
setwd("~/intro_workflow_seurat")
# check 
getwd()

#Load PBMC Data  
pbmc.data <- Read10X(data.dir = "pbmc3k/filtered_gene_bc_matrices/hg19")

## Initialize the Seurat object with the raw (non-normalized data).
# min.cells =3, means keep only genes that are expressed in at least 3 cells
#
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# get info on the object 
pbmc

# Examine a few genes in the first 30 cells (see how many transcripts of each gene are in each cell): 
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# all the . are zeros --> use sparse matrix representation
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
# this outputs 709591472 bytes

sparse.size <- object.size(pbmc.data)
sparse.size
# this outputs 29905192 bytes

dense.size/sparse.size
# this outputs 23.7 bytes


# PREPROCESSING START
# cell selection/filtration based on QC, normalisation, scaling, detection of highly variable features 
# example: add columns to object with [[ - % mitochondrial genome - to find dying cells 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# where does Seurat store QC data? 
# automatically stores # of unique genes & total molecules when CreateSeuratObject()
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

#Visualise QC metrics as a Violin Plot 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# use FeatureScatter to see feature-feature relationship 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# how to adjust dataset, to remove the low quality cells? 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# says, keep only cells with # of genes detected btw 2.5k and 200, and with 
# % reads of mitochondrial genome < 5% 




# NORMALISATION 
# normalise data (so that reads are comparable, cells/samples sequenced at different depths)
#use LogNormalize 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor =10000)
# not necessary to specify all parameters, can also just do: 
# pbmc <- NormalizeData(pbmc) ]



# HIGHLY VARIABLE GENES (FEATURE) SELECTION
# select ~ 2k HVGs with FindVariableFeatures() function
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features, both with and without labels 
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# SCALING THE DATA 
# we don't want the top HVGs to dominate further analysis (PCA)
# standardize each gene --> mean = 0 & variance = 1 ACROSS CELLS
# ScaleData() function
# create vector that includes names of all genes (by default are row names)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


# PCA - Linear Dimensionality Reduction 
# input = previously determined highly variable features by default 
# if we want different subset of features, can specify 
# for 1st PC -  Seurat outputs list of genes w most + / - loadings
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
# results saved under slot 'pca' in pbmc object
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# visualisation
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#Seurat picks most “extreme” 500 cells (both sides of the PC axis)

#heatmap for 15 PC's 
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# CLUSTERING - we cluster ALL the cells, based on distances calculated w PCA 
# so different PC's determine how close / far apart we plot cells 
# 1st step - build a graph of how cells are related (based on 10 PC's in this case)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#2nd Step - partition the graph into communities (i.e. clusters)
pbmc <- FindClusters(pbmc, resolution = 0.5)
Idents(pbmc) <- "seurat_clusters"

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# run UMAP or TSNE (non linear dim reduction) - visualise formed clusters
pbmc <- RunUMAP(pbmc, dims = 1:10) # always with 10 PCs like in clustering 

# plot
# set `label = TRUE` or use LabelClusters function to label individual clusters
DimPlot(pbmc, reduction = "umap")
# save this object so no need to repeat previous steps 
saveRDS(pbmc, file = "../intro_workflow_seurat/pbmc_tutorial.rds")


# DIFFERENTIAL EXPRESSION 
# goal= find differentially expressed markers (genes/features) that define clusters
# find top 5 MARKERS in CLUSTER 2 (comparing this cluster against all others)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# seurat order the genes/markers by level of significance (p value)

# distinguishing CLUSTER 5 FROM CLUSTER 0 and 3
# good for specific comparisons 
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0 ,3))
# look at top 5 marker genes that distinguish most significantly cluster 5 from 0&3
head(cluster5.markers, n = 5)

# FINDALLMARKERS- find markers for every cluster compared to all remaining cells
# report only positive (upregulated) by setting only.pos = TRUE
# filter to keep only genes with log2FC > 1 and group markers by cluster  
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# and now group by cluster and filter to only have values with log2FC>1
pbmc.markers %>%  
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# use a different test for Differential Expression - "ROC" 
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)


# VISUALISING MARKERS 
# violin plot 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# plot the raw counts 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), slot = "counts", log = TRUE)

# features plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# heatmap for the top 20 markers 
pbmc.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10 

# top 10 was defined as a data frame during HVG selection
# FindMarkers also defined as data frame marker selection

DoHeatmap(pbmc, features = top10$gene) + NoLegend()


# assign cell type identity for each cluster
# use atlases for canonical / well known genes, then annotate clusters w ID's
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# make improved UMAP and save it to the disk 
library(ggplot2)
# save it into 'plot' 
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot

# save it 
saveRDS(pbmc, file = "../intro_workflow_seurat/pbmc3k_final.rds")

