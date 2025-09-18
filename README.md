#  README

## 1. Data sources

This task is based on publicly available sequencing data from the study "Tubule-Specific Compensatory Responses to Cpt1a Deletion in Aged Mice". The dataset is based on the GEO submission `GSE277335` and the workflow is based on a H5 count matrix file.

The data was originally sequenced using a Illumina NovaSeq 6000.

## 2. How to download

The original samples can be download from SRA using...

```bash
SRR_FILES=(GSM7863646)
prefetch ${SRR_FILES[@]}

```bash
for SRR in "${SRR_FILES[@]}"; do
  fasterq-dump Raw/${SRR}/${SRR}.sra \
    --split-files \
    --include-technical \
    -O Fastq \
    --threads $THREADS
  pigz -p $THREADS Fastq/${SRR}_1.fastq
  pigz -p $THREADS Fastq/${SRR}_2.fastq
  pigz -p $THREADS Fastq/${SRR}_3.fastq
done
```

## 3. Pre-processing

In the original paper, samples were processed using cellranger:

```bash
SRR_FILES=(GSM7863646)

cellranger mkref \
  --genome=GRCm39 \
  --fasta=Reference_ch10/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --genes=Reference_ch10/Mus_musculus.GRCm39.112.gtf

for SRR in "${SRR_FILES[@]}"; do
    cellranger count \
    --id=$SRR \
    --transcriptome=GRCm39 \
    --fastqs=Fastq \
    --sample=$SRR \
    --create-bam=false
done

cellranger aggr --id=Combined --csv=aggr.csv --normalize=mapped

```

## 4. How the workflow works

The workflow files are stored in `workflow/`. The files needed are `workflow/1_Workflow.R` and `workflow/2_Answer_Questions.R` 

### Step 1 - Converting to Seurat object and QC 

**Purpose:** Generating a seurat object and quality contorl
**Tools:** `Seurat`
**Inputs:** H5 file containing the count matrix
**Outputs:** seurat object
**Command:**

This step follow the original paper workflow: a single seurat object was created and
filtered to include only cells containing 750-5000 genes per cell, 12,500 maximum transcripts
per cell, and less than 0.5% mitochondrial gene content, 0.05% hemoglobin gene content, and
5% ribosomal subunit-associated gene content. All mitochondrial and hemoglobin genes were
deleted after inclusion filters.

```R
library(Seurat)

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
```
### Step 2 - Normalisation

**Purpose:** seurat log-normalisation
**Tools:** `Seurat`
**Inputs:** Need the filtered seurat object
**Outputs:** seurat object normalized
**Command:**

Cells were then log-transformed and scaled with determination of variable features and PCA analysis using default settings in Seurat.

```bash
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), seed.use = 1234)
```

### Step 3 - Clustering & UMAP

**Purpose:** Clustering and dimensionality reduction 
**Tools:** `Seurat`
**Inputs:** need the filtered and normalized seurat object
**Outputs:** seurat object
**Command:**

```bash
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, seed.use = 1234)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, seed.use = 1234)
```
### Step 4 - Quantification

**Purpose:** Generate markers per cluster
**Tools:** `Seurat`
**Inputs:** need the filtered and normalized seurat object
**Outputs:** seurat object
**Command:**

```bash
markers <- FindAllMarkers(seurat_obj)
```