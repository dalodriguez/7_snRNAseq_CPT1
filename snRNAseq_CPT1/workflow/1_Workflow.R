set.seed(1234)
library(Seurat)
library(hdf5r)
library(dplyr)

# CREATE SEURAT OBJECT
data <- Read10X_h5("Smpl2_filtered_feature_bc_matrix.h5")
seurat_obj <- CreateSeuratObject(counts = data, project = "snRNAseq")

# FILTERING AND QC
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-") # Mitochondrial genes
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]") # Hemoglobin genes
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]") # Ribosomal genes

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= 750 & nFeature_RNA <= 5000 &
           nCount_RNA <= 12500 &
           percent.mt < 0.5 &
           percent.hb < 0.05 &
           percent.ribo < 5
)

# REMOVING MITOCHONDRIAL AND HEMOGLOBIN GENES FROM SEURAT OBJECT
mito_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE)
hb_genes   <- grep("^Hb[ab]", rownames(seurat_obj), value = TRUE)

seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), c(mito_genes, hb_genes)))

# NORMALIZATION AND DIMENSIONALITY REDUCTION
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), seed.use = 1234)

# CLUSTERING & UMAP
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, seed.use = 1234)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, seed.use = 1234)

# SAVE SEURAT OBJECT
saveRDS(seurat_obj, "seurat_object.rds")

# CLUSTER MARKERS
markers <- FindAllMarkers(seurat_obj)