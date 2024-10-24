## load environment and dataset 
library(Seurat)
library(ggplot2)
library(scCustomize)
library(readr)
library(pheatmap)
library(matrixStats)
library(spdep)
library(geojsonR)

## set directory 
data_dir <- "/project/shared/spatial_data_camp/datasets/DATASET2/XENIUM_COLORECTAL_CANCER"

data <- ReadXenium(data_dir, outs = c("matrix", "microns"), type=c("centroids", "segmentations"))

names(data$matrix)

## read in extra information on the cells from csv file 
cell_meta_data <- read.csv(file.path(data_dir, "cells.csv.gz"))
rownames(cell_meta_data) <- cell_meta_data$cell_id
head(cell_meta_data)

## creating seurat object from the data 
seurat <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]],
                             assay = "XENIUM_CANCER",
                             meta.data = cell_meta_data)
seurat

## creating a FOV 
coords <- CreateFOV(coords = list(centroids = CreateCentroids(data$centroids), 
                                  segmentation = CreateSegmentation(data$segmentations)),
                    type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "XENIUM_CANCER")
seurat[["COLON"]] <- coords


## visualise the seurat object 
ImageFeaturePlot(seurat, 'nCount_XENIUM_CANCER', axes = T) + scale_fill_viridis()


## subset the seurat object based on centroid coordinates of region of interest 

cell_meta_data$x_centroid
range(cell_meta_data$x_centroid) ## seeing range of x axis

# create filters for seurat object on both x and y axis, and then apply the filters to subset on this 

seurat$x_FILTER <- seurat$x_centroid <= 6000 & seurat$x_centroid >= 1000
seurat$y_FILTER <- seurat$y_centroid <= 4000 & seurat$y_centroid >= 2000

seurat_subset <- subset(seurat, x_FILTER & y_FILTER)

## visualise the subsetted area selected 
ImageFeaturePlot(seurat_subset, 'nCount_XENIUM_CANCER', axes = T) + scale_fill_viridis()

## cleaning up the data - removing unused codewords and negative probes - have to filter these datasets also 
length(colnames(seurat_subset))
length(colnames(data$matrix[["Negative Control Codeword"]]))
length(intersect(colnames(seurat_subset), colnames(data$matrix[["Negative Control Codeword"]])))

seurat_subset[["Negative.Control.Codeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]][, colnames(seurat_subset)])
seurat_subset[["Negative.Control.Probe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]][, colnames(seurat_subset)])
seurat_subset[["Unassigned.Codeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]][, colnames(seurat_subset)])

seurat_subset

## QC AND VISUALISATION
## visualising number of transcripts per cell in space 
ImageFeaturePlot(seurat_subset, 'nCount_XENIUM_CANCER', axes = T) + scale_fill_viridis()

## visualising the number of genes per cell in space 
ImageFeaturePlot(seurat_subset, "nFeature_XENIUM_CANCER") + scale_fill_viridis_c()


## Examining the number of genes (features) per cell using a density plot - analysing distribution across object
ggplot(seurat_subset[[]], aes(nFeature_XENIUM_CANCER)) + geom_density() +
geom_vline(xintercept = quantile(seurat_subset$nFeature_XENIUM_CANCER, c(0.01, 0.1, 0.5, 0.9, 0.99)), lty=2)
quantile(seurat_subset$nFeature_XENIUM_CANCER, c(0.01, 0.1, 0.5, 0.9, 0.99)) 


## interrogate cell-to-nucleus area - provides insight on morphology (biological perturbations in disease vs control states)
ImageFeaturePlot(seurat_subset, "cell_area") + scale_fill_viridis_c()

seurat_subset$cell_nucleus_ratio <- seurat_subset$nucleus_area / seurat_subset$cell_area
ImageFeaturePlot(seurat_subset, "cell_nucleus_ratio") + scale_fill_viridis_c()

## plot cell area to see size of cells, and also plot transcript vs cell size 
ggplot(seurat_subset[[]], aes(cell_area)) + geom_density() +
  geom_vline(xintercept = quantile(seurat_subset$cell_area, c(0.01, 0.9, 0.9985)), lty=2)## seeing the quantiles to make sure we aren't removing too many cells

ggplot(seurat_subset[[]], aes(nCount_XENIUM_CANCER, cell_area)) + geom_point() 

## Filter large cells out - take 99th percentile as we do not want to remove too many cells 
seurat_subset[["SIZE_FILTER_LARGE"]] <- seurat_subset$cell_area < quantile(seurat_subset$cell_area, .9985)

## visualise the cells which have been flagged for filtering 
ImageDimPlot(seurat_subset, group.by = "SIZE_FILTER_LARGE")

## Filter very small cells - could be segmentation artefacts 
seurat_subset[["SIZE_FILTER_SMALL"]] <- seurat_subset$cell_area > quantile(seurat_subset$cell_area, .01)

ImageDimPlot(seurat_subset, group.by = "SIZE_FILTER_SMALL")

## sanity checking to see how cell size correlates with gene detection rate 

p1 <- VlnPlot(seurat_subset, "nFeature_XENIUM_CANCER", group.by = "SIZE_FILTER_SMALL", pt.size = .05, alpha = .5) + labs(title="Small Cell Filter")
p2 <- VlnPlot(seurat_subset, "nFeature_XENIUM_CANCER", group.by = "SIZE_FILTER_LARGE", pt.size = .1, alpha = .5)+ labs(title="Large Cell Filter")

p1 + p2


## overall transcript detection, setting the filter at > 15 transcripts per cell 
seurat_subset$TRANSCRIPT_FILTER <- seurat_subset$nCount_XENIUM_CANCER >= 15

## visualise the cells we will be losing by applying this filter
ImageDimPlot(seurat_subset, group.by="TRANSCRIPT_FILTER")


## visualise and filter out the negative controls to understand artifacts and noise in the data 
ImageFeaturePlot(seurat_subset, "nCount_Negative.Control.Codeword") + scale_fill_viridis_c()
ImageFeaturePlot(seurat_subset, "nCount_Negative.Control.Probe") + scale_fill_viridis_c()
ImageFeaturePlot(seurat_subset, "nCount_Unassigned.Codeword") + scale_fill_viridis_c()

seurat_subset$PROBE_FILTER <- seurat_subset$nCount_Unassigned.Codeword == 0 &
  seurat_subset$nCount_Negative.Control.Codeword == 0 &
  seurat_subset$nCount_Negative.Control.Probe == 0

ImageDimPlot(seurat_subset, group.by="PROBE_FILTER")

## Filter seurat object by all the filters we have created thus far
seurat_subset <- subset(seurat_subset, PROBE_FILTER & SIZE_FILTER_LARGE & SIZE_FILTER_SMALL & TRANSCRIPT_FILTER)

## examine the cleaned up object 
seurat_subset


## NORMALISATION 
## Using SCTransform, setting clip-range parameter to limit the range of the transformed values, to help downstream analyses
seurat_subset <- SCTransform(seurat_subset, assay = "XENIUM_CANCER", clip.range = c(-10, 10))

## Principal component analysis 
seurat_subset <- RunPCA(seurat_subset)

ElbowPlot(seurat_subset, 50)

## plotting the top genes, to identify the genes which are the most influential to the principal components
PC_Plotting(seurat_subset, dim_number = 1)

## Visualising the expression of a specific gene across cells in a PC 
FeaturePlot(seurat_subset, "SOX9", reduction = "pca") + scale_color_viridis_c()

## Can visualise how the principal components are distributed in space 
ImageFeaturePlot(seurat_subset, "PC_1") + scale_fill_viridis_c()


## Generate the UMAP - to visualise dimensionality for clustering 
seurat_subset <- RunUMAP(seurat_subset, dims = 1:20)
seurat_subset <- FindNeighbors(seurat_subset, reduction = "pca", dims = 1:20)
seurat_subset <- FindClusters(seurat_subset, resolution = 0.7)

## can change resolution by setting = i and pre-defining i, running in seq
## for ( i in seq(0.5:1, by = 0.25))
## seurat_subset <- FindClusters(seurat_subset, resolution = i)

## visualise the UMAP
DimPlot(seurat_subset, label=T, repel=T)

## visualise the clusters in space 
ImageDimPlot(seurat_subset, size= 0.5)

## Identifying marker genes which are DE for specific cell clusters 
markers <- FindMarkers(seurat_subset, ident.1 = 0, max.cells.per.ident = 500)
head(markers) ## gives us the highest DE genes i.e. our markers which we can visualise our clusters based on 

FeaturePlot(seurat_subset, "RRM2", label=T, repel=T)+ scale_color_viridis_c(direction=-1)
FeaturePlot(seurat_subset, "TYMS", label=T, repel=T)+  scale_color_viridis_c(direction=-1)
FeaturePlot(seurat_subset, "HMGB2", label=T, repel=T)+ scale_color_viridis_c(direction=-1)
FeaturePlot(seurat_subset, "C1QBP", label=T, repel=T)+ scale_color_viridis_c(direction=-1)
FeaturePlot(seurat_subset, "CA2", label=T, repel=T)+ scale_color_viridis_c(direction=-1)


## visualise the top markers for every cluster (not specifying a cluster of interest (ident.1 function))
markers <- FindAllMarkers(seurat_subset, max.cells.per.ident = 500)
head(markers)


## extracting the top 5 marker genes for each cluster using Extract_Top_Markers function 
top <- Extract_Top_Markers(markers, num_genes = 5, named_vector = FALSE, make_unique = TRUE)
top

## Plot the expression data and simultaneously clusters the genes and groups for enhanced visual interpretation
Clustered_DotPlot(seurat_subset, features = top, k=18)  ## k=18 determines the number of clusters for hierarchical clustering


## cell type identification (read in healthy control dataset)
ref <- readRDS("/project/shared/spatial_data_camp/datasets/SINGLE_CELL_REFERENCES/COLON_HC_5K_CELLS.RDS")
ref
DimPlot(ref)


## We want to evaluate how much structural information is lost in single-cell data when limiting ourselves to the targeted gene set. 
## Accurate cluster prediction is challenging if the current gene set does not adequately identify them. 
## To do this, we will quickly re-embedd the data using only the genes present in our spatial transcriptomics data 
## and keep the original cluster annotations derived from unbiased data.
ref <- SCTransform(ref, residual.features =rownames(seurat_subset))
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims=1:20)
DimPlot(ref, label=T, repel=T)

## Visualise the gene panel across our single cell reference clusters, to see how coverage is spread across cell types
ps <- AggregateExpression(ref, features = rownames(seurat_subset), normalization.method = "LogNormalize", assays="RNA", return.seurat = T)
ps <- ScaleData(ps, features=rownames(ps))
pheatmap(LayerData(ps, layer="scale.data"), show_rownames = F)






