# Integrated Lineage Tracing and Single-Cell Multiomics Identifies a Plasticity-Proliferation Axis That Drives Treatment Adaptation in Hepatoblastoma

## Abstract
Hepatoblastoma (HB) has one of the lowest mutational burdens among childhood cancers, limiting the role of genetic selection. Nevertheless ~20% of patients relapse, implicating non-genetic mechanisms, such as phenotypic plasticity, in treatment adaptation. The temporal and clonal dynamics of post-treatment plasticity in HB remain poorly defined. Here, we integrate expressed DNA barcoding approaches with single-cell multiomics in preclinical models to simultaneously trace clonal and phenotypic dynamics following cisplatin treatment.
We observe that cisplatin selects for progenitor-like states, which persist as reservoirs for phenotypic re-diversification. Barcode lineage tracing reveals a subset of persister clones stochastically resume proliferation and transition from less- to more-differentiated phenotypes, an observation also captured in patient data. Pseudotime and landscape approaches identify a plasticity-proliferation axis underlying recovery post-treatment, driven by coordinated activation of E2F transcription factors and their downstream target BIRC5 (survivin). Downregulation and pharmacological inhibition of BIRC5 disrupt this axis, shifting phenotypic dynamics towards less-differentiated states and killing cisplatin-persister cells, supporting a role for BIRC5 in plasticity-led awakening from persistence. Consistently, elevated BIRC5 expression is associated with poor patient outcome.
Together, these findings establish a mechanistic link between persistence, phenotypic plasticity, and stochastic clonal outgrowth in HB, and identify BIRC5 as a regulator of plasticity-driven adaptation and a therapeutically actionable vulnerability to halt treatment adaptation.

## Analysis scripts

**HuH6 scRNA-seq analysis** 
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` and partitioning around medoids (PAM) clustering
   - Trajectory inference using `Slingshot`
   - Transcritpion factor motif analysis of pseudotime-derived genes with `RcisTarget`
   - Integration with barcode expression data

**HuH6 DNA barcode analysis**
   - Filtering and plotting of barcode abundance by counts in the DNA sequencing data
   - Calculation of beta diversity using `vegan`  

**HuH6 single-cell multiomics analysis**
   - snRNA-seq analysis with `Seurat`, analogous to as outlined above for the analysis of HuH6 scRNA-seq data
   - snATAC-seq analysis using `Signac` with MACS3 peak calling

**HepG2 scRNA-seq analysis**
   - (As outlined above for HuH6 scRNA-seq analysis)

**Patient snRNA-seq analysis**
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` clustering
   - Cell type inference with SingleR annotation
   - Analysis of malignant cells using InferCNV

**BIRC5 siRNA HB303, HepG2, HuH6 scRNA-seq analysis**
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` clustering
