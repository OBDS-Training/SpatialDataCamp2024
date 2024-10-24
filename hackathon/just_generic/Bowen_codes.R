library(Seurat)
library(ggplot2)
library(scCustomize)
library(readr)
library(pheatmap)
library(matrixStats)
library(spdep)
library(geojsonR)


# load data
data_dir <- "/project/shared/spatial_data_camp/datasets/DATASET3/XENIUM_5K_LYMPH_NODE"
data = ReadXenium(data_dir, outs = c("matrix"), type=c("centroids", "segmentations"))
cell_meta_data <- read.csv(file.path(data_dir, "cells.csv.gz"))
rownames(cell_meta_data) <- cell_meta_data$cell_id
head(cell_meta_data)


# create seurat object
seurat <- CreateSeuratObject(counts = data$`Gene Expression`,
                             assay = "XENIUM",
                             meta.data = cell_meta_data)
coords <- CreateFOV(coords = list(centroids = CreateCentroids(data1$centroids), 
                                  segmentation = CreateSegmentation(data1$segmentations)),
                    type = c("segmentation", "centroids"),
                    molecules = data1$microns,
                    assay = "XENIUM")
seurat[["LYMPH"]] <- coords  
seurat
seurat[["Negative.Control.Codeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
seurat[["Negative.Control.Probe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
seurat[["Unassigned.Codeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
seurat


# QC
ImageFeaturePlot(seurat, "nCount_XENIUM", max.cutoff = 'q99') + scale_fill_viridis_c()
ImageFeaturePlot(seurat, "nFeature_XENIUM", max.cutoff = 'q99') + scale_fill_viridis_c()
ggplot(seurat[[]], aes(nFeature_XENIUM)) + geom_density()
quantile(seurat$nFeature_XENIUM, c(0.01, 0.1, 0.5, 0.9, 0.99))


# backup naive seurat object
setwd("/project/shared/spatial_data_camp/HACKATHON/just_generic")
saveRDS(seurat, file = "seurat_BZ.RDS")


# clear env
rm(coords, cell_meta_data)
gc()


#continue with QC
ImageFeaturePlot(seurat, "cell_area", max.cutoff = 'q99') + scale_fill_viridis_c()
seurat$cell_nucleus_ratio <- seurat$nucleus_area / seurat$cell_area
ImageFeaturePlot(seurat, "cell_nucleus_ratio", min.cutoff = 'q1',  max.cutoff = 'q99') + scale_fill_viridis_c()
ggplot(seurat[[]], aes(cell_area)) + geom_density()
ggplot(seurat[[]], aes(nCount_XENIUM, cell_area)) + geom_point() 
ggplot(seurat[[]], aes(cell_nucleus_ratio)) + geom_density()


# remove large and small cells
seurat[["SIZE_FILTER_LARGE"]] <- seurat$cell_area < quantile(seurat$cell_area, .99)
ImageDimPlot(seurat, group.by="SIZE_FILTER_LARGE")
seurat[["SIZE_FILTER_SMALL"]] <- seurat$cell_area > quantile(seurat$cell_area, .01)
ImageDimPlot(seurat, group.by="SIZE_FILTER_SMALL")
p1 <- VlnPlot(seurat, "nFeature_XENIUM", group.by = "SIZE_FILTER_SMALL", pt.size = .1, alpha = .5) + labs(title="Small Cell Filter")
p2 <- VlnPlot(seurat, "nFeature_XENIUM", group.by = "SIZE_FILTER_LARGE", pt.size = .1, alpha = .5)+ labs(title="Large Cell Filter")


# filter cells by transcript number
p1 + p2
seurat[["SIZE_FILTER_SMALL"]] <- seurat$cell_area > quantile(seurat$cell_area, .1)
seurat$TRANSCRIPT_FILTER_20 <- seurat$nCount_XENIUM >= 20
ImageDimPlot(seurat, group.by="TRANSCRIPT_FILTER_20")


# filter by negative control
ImageFeaturePlot(seurat, "nCount_Negative.Control.Codeword") + scale_fill_viridis_c()
ImageFeaturePlot(seurat, "nCount_Negative.Control.Probe") + scale_fill_viridis_c()
ImageFeaturePlot(seurat, "nCount_Unassigned.Codeword") + scale_fill_viridis_c()
seurat$PROBE_FILTER <- seurat$nCount_Unassigned.Codeword == 0 &
  seurat$nCount_Negative.Control.Codeword == 0 &
  seurat$nCount_Negative.Control.Probe == 0
ImageDimPlot(seurat, group.by = "PROBE_FILTER")

saveRDS(seurat, file = "seurat_QC_BZ.RDS")

# apply all filters
seurat <- subset(seurat, PROBE_FILTER & SIZE_FILTER_LARGE & SIZE_FILTER_SMALL & TRANSCRIPT_FILTER)
saveRDS(seurat, file = "seurat_filter_BZ.RDS")


# transform data & runPCA
seurat <- SCTransform(seurat, assay = "XENIUM", clip.range = c(-10, 10))
saveRDS(seurat, file = "seurat_SCT_BZ.RDS")
seurat <- RunPCA(seurat)


# analyse PCA
ElbowPlot(seurat, 50)
ImageFeaturePlot(seurat, "PC_1", min.cutoff = 'q01') + scale_fill_viridis_c()


# run UMAP and clustering
seurat <- RunUMAP(seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.7)
DimPlot(seurat, label=T, repel=T)
ImageDimPlot(seurat, size=.5)
saveRDS(seurat, file = "seurat_CLUSTER_BZ.RDS")


# reading Tabula refdata
library(SeuratDisk)
Convert("/project/shared/spatial_data_camp/HACKATHON/just_generic/HCA_dataset/4de7ab62-475c-43f9-a1d4-bfd693e7df07/TabulaSapiens.h5ad", 
        dest = "h5seurat", overwrite = TRUE)
Tabula <- LoadH5Seurat("/project/shared/spatial_data_camp/HACKATHON/just_generic/HCA_dataset/4de7ab62-475c-43f9-a1d4-bfd693e7df07/TabulaSapiens.h5seurat",
                       misc = F)


# subset to lymph node
Tabula_Ly = subset(Tabula, organ_tissue == "Lymph_Node")
saveRDS(Tabula_Ly, file = "Tabula_lymph_node.RDS")


# normalise ref data 
Tabula_Ly <- SCTransform(Tabula_Ly, residual.features =rownames(seurat))
Tabula_Ly <- RunPCA(Tabula_Ly)
Tabula_Ly <- RunUMAP(Tabula_Ly, dims=1:20)
DimPlot(Tabula_Ly, label=T, repel=T, label.box = T,group.by = "cell_ontology_class")+NoLegend()


# transfer label using anchoring
anchors <- FindTransferAnchors(reference = Tabula_Ly, 
                               query = seurat, 
                               normalization.method = "SCT")
saveRDS(anchors, file = "Tabula_anchor.RDS")
seurat <- TransferData(anchorset = anchors, 
                       refdata = Tabula_Ly$cell_ontology_class, 
                       prediction.assay = TRUE,
                       weight.reduction = seurat[["pca"]], 
                       query = seurat, 
                       dims=1:30)


# check prediction results
DimPlot(seurat, group.by = "predicted.id", label = T, label.box = T, repel = T) + NoLegend()
ImageDimPlot(seurat, group.by = "predicted.id")
FeaturePlot(seurat, features = "predicted.id.score")
ImageFeaturePlot(seurat, features =  "predicted.id.score")

# export predictions

predictions = as.data.frame(seurat$predicted.id)
predictions$predictions.score = seurat$predicted.id.score
head(predictions)
names(predictions) = paste(c("predicted.id", "predicted.id.score"))
write.csv(predictions, file = "Tabula_predicted_labels.csv")



# spatially cluster the cells
coords <- GetTissueCoordinates(seurat, which = "centroids")
rownames(coords) <- coords$cell
neighbours <- FindNeighbors(as.matrix(coords[, c("x", "y")]), k.param = 20, return.neighbor=TRUE)



# find what cells are around B cells
cells <- colnames(subset(seurat, predicted.id == "b cell"))
adjacent <- TopNeighbors(neighbours, cells, n = 10)

Cells_by_Bcell = View(table(seurat$predicted.id[adjacent]))

Idents(seurat) <- "Other Cells"
seurat <- SetIdent(seurat, cells = adjacent, "Adjacent Cells")
seurat <- SetIdent(seurat, cells = cells, "Cells of Interest")


# play with segmentation

