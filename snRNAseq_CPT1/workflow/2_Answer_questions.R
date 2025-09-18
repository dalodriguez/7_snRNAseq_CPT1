# QUESTION 1 ANSWER:
seurat <- readRDS('seurat_object.rds')
dim(seurat)

# QUESTION 2 ANSWER:
seurat <- readRDS('seurat_object.rds')
seurat <- FindVariableFeatures(seurat)
variable_genes <- VariableFeatures(seurat)
top3_genes <- head(variable_genes, 3)
top3_genes

# QUESTION 3 ANSWER
seurat <- readRDS('seurat_object.rds')
length(levels(seurat))

# QUESTION 4 ANSWER
seurat <- readRDS('seurat_object.rds')
data <- FetchData(seurat, vars = c("Slit2"))
data <- data[data$Slit2 > 0, ,drop=FALSE]

# QUESTION 5 ANSWER
seurat <- readRDS('seurat_object.rds')
cells_of_interest <- c("AAACCCAAGTTGCGCC-1", "AAACCCAGTCGCCACA-1")
expr_data <- as.data.frame(GetAssayData(seurat, slot = "data")[, cells_of_interest, drop = FALSE])
expr_data <- expr_data[expr_data$AAACCCAAGTTGCGCC-1 > 0, , drop = FALSE]
expr_data <- expr_data[expr_data$AAACCCAGTCGCCACA-1 > 0, , drop = FALSE]
dim(expr_data)
