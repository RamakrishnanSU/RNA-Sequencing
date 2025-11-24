



## 🧬 RNA-Seq Analysis of Breast Tumor Subtypes and Healthy Samples

This repository contains the code, analysis scripts, and resulting figures for an RNA Sequencing (RNA-Seq) project focused on identifying differentially expressed genes (DEGs) between different breast cancer subtypes (Triple-Negative Breast Cancer - **TNBC**, Non-TNBC, and HER2-positive) and healthy control samples.

The analysis pipeline was designed and executed on the IBU Cluster using dedicated bioinformatics containers (Apptainer/Singularity).

-----

### 📝 Project Overview

This study utilized publicly available RNA-Seq data from Eswaran et al. (2012) [cite\_start](GEO accession **GSE52194**) to explore subtype-specific gene expression patterns, which is critical for developing targeted therapies and understanding cancer heterogeneity[cite: 13, 18, 20, 24, 25, 30].

The dataset comprises **12 paired-end RNA-Seq libraries**, with three biological replicates for each condition:

  * **TNBC** (Triple-Negative Breast Cancer)
  * **Non-TNBC** (Non-Triple-Negative Breast Cancer)
  * **HER2** (HER2-positive Breast Cancer)
  * **Normal** (Healthy Control Samples) 

-----

### 🚀 Analysis Workflow

The entire bioinformatics pipeline was executed in a controlled environment using SLURM job arrays and Apptainer containers.

1.  **Quality Control (QC)**: Performed using **FastQC** to assess base quality, adapter contamination, and sequence duplication rates.
2.  **Read Alignment**: Reads were mapped to the $Homo\_sapiens.GRCh38.113$ reference genome using **HISAT2**.
      * *Result Highlight*: Normal samples showed high overall alignment rates $(>96\%)$ and high concordance $(>80\%)$, while tumor samples (TNBC, Non-TNBC, HER2) had lower concordance rates (30% to 50%), likely due to the structural rearrangements typical of cancer genomes.
3. **Read Counting**: Gene quantification was performed using **featureCounts** to determine the number of reads overlapping with annotated genes.
4. **Differential Expression Analysis (DEA)**: The raw count matrix was analyzed using the **DESeq2** R package to identify differentially expressed genes between the tumor subtypes and the Normal controls.

-----

### 📁 Repository Contents

| File/Directory | Description |
| :--- | :--- |
| `R_Code.txt` | The R script containing the DESeq2, PCA, and Heatmap generation code. |
| `Rna_Sequencing.sh` | SLURM submission scripts for HISAT2 indexing, alignment, sorting, indexing BAM files, and featureCounts. |
| `RNA-Seq Analysis Report.pdf` | The final project report detailing the materials, methods, results, and discussion. |
| `normalized_counts_matrix.csv` | Output file from the R script: Variance-Stabilized Transformed (VST) normalized gene expression counts. |
| `differential_expression_results.csv` | Output file from the R script: Full DESeq2 results comparing **TNBC** vs **Normal**. |
| `significant_genes.csv` | Output file from the R script: Subset of DEGs with an adjusted p-value $(padj < 0.05)$ |
| `PCA_plot.jpeg` | Output figure visualizing sample clustering based on the top two principal components.  |
| `heatmap_top_variable_genes.jpeg` | Output figure visualizing the expression of the top 50 most variable genes across all samples.  |

-----

### 📊 Key Findings

#### Sample Clustering

The Principal Component Analysis (PCA) plot (Figure 1 in the report) and the Heatmap (Figure 2 in the report) show that samples cluster distinctly based on their condition (HER2, NonTNBC, Normal, TNBC).

  * **PC1** accounts for **59%** of the variance, and **PC2** accounts for **16%**, demonstrating that gene expression profiles are highly condition-specific, making the data suitable for downstream differential analysis.

#### Differential Expression

Comparing the different subtypes, a total of **12,783 DEGs** were identified, with **4,541 upregulated** and **8,242 downregulated** genes.

  * **ENSG00000139618**: Showed **significantly higher expression in TNBC** compared to other subtypes, suggesting its potential role in the aggressive nature of TNBC.
  * **ENSG000001414510**: Exhibited **elevated expression in Normal and NonTNBC** samples but was lower in HER2 and TNBC, indicating its potential involvement in normal cellular processes and suppression in tumorigenic conditions.

#### Functional Enrichment

Overrepresentation analysis of the DEGs highlighted key biological processes, including:

  * **Nucleosome assembly** (GO:0006334)
  * **Cellular response to bacterial molecules** (GO:0071219, GO:0002237)
  * **Response to lipopolysaccharide** (GO:0032496)
  * **Cytoplasmic translation** (GO:0002181).

These terms suggest involvement in **immune responses**, **chromatin remodeling**, and **protein synthesis**. TNBC specifically showed upregulation in **cell cycle** and **DNA repair** genes, consistent with its aggressive phenotype.

-----

### 🛠️ Software and Tools

| Tool | Purpose | Container Path (from report) |
| :--- | :--- | :--- |
| **FastQC** | Initial quality control | `/containers/apptainer/fastqc-0.12.1.sif`  |
| **HISAT2** | Alignment to reference genome | `/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif`  |
| **Samtools** | SAM/BAM file sorting and indexing |Included in HISAT2 container |
| **featureCounts** | Read quantification (gene counting) | `subread_2.0.1--hed695b0_0.sif`  |
| **DESeq2** (R package) | Differential expression analysis | R environment |
| **pheatmap** (R package) | Heatmap generation | R environment |
| **ggplot2** (R package) | PCA plot generation | R environment |

-----

### 🧑‍💻 Code Snippet Example (DESeq2 Analysis)

The main differential expression analysis was performed using the following core steps in R (from `Rscript.txt`):

```r
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Proceed with DESeq2 analysis
dds <- DESeq(dds)

# Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Perform differential expression analysis (e.g., TNBC vs Normal)
res <- results(dds, contrast = c("condition", "TNBC", "Normal"))
```

-----

### 🔗 References

The analysis utilized data and methods from the following sources:

1.  Eswaran, J., et al. (2012). RNA sequencing of breast cancer samples. *Journal of Cancer Genomics, 5(3), 123-135*. 
2.  Love, M. I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-Seq data with DESeq2. *Genome Biology, 15(12), 550*. 

-----
