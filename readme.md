# Integrated Lineage Tracing and Single-cell Multiomics Identifies Plasticity Drivers of Adaptation to Cisplatin Treatment in Hepatoblastoma

## Abstract
Background and Aims: Hepatoblastoma (HB), the most common paediatric liver cancer, has one of the lowest mutational burdens among childhood malignancies. Despite its genomic stability, approximately 20% of HB patients develop chemoresistance and relapse, suggesting a prominent role for non-genetic cancer mechanisms in treatment adaptation. The role and underlying mechanisms of plasticity in HB requires investigation after treatment at a more granular temporal scale to further elucidate its role in treatment adaptation. In this study, we investigate the evolutionary principles and molecular drivers of phenotypic plasticity during the early phases of HB adaptation to cisplatin, a standard of care chemotherapy in HB treatment.<br />
Approach and Results: By integrating expressed DNA barcoding with single-cell transcriptomics and epigenomics in dynamic preclinical models, we simultaneously trace the evolution of single clones and phenotypes of single cells. Before treatment, we identify phenotypic states in HB preclinical models that reflect the stages of liver development. We uncover progenitor-like states persist after treatment, but upon drug withdrawal clones stochastically resume proliferation and transition from less-to-more differentiated-like phenotypes, recapitulating the phenotypic heterogeneity observed at baseline. We found that this plasticity-mediated awakening of clones is driven by transcription factors in the MYC and E2F pathways, which regulate the expression of plasticity regulators, including BIRC5. Notably, chromatin at the BIRC5 locus becomes more accessible in post-treatment populations, suggesting an epigenetic memory that could prime cells for future cisplatin challenges. High BIRC5 expression is associated with poor prognosis in HB patients, underscoring the clinical relevance of plasticity-driven treatment adaptation.<br />
Conclusions: These findings reveal that phenotypic plasticity, regulated through transcriptional and epigenetic mechanisms, underpins clonal survival and repopulation following chemotherapy in HB. Targeting plasticity effectors, such as BIRC5, may offer new therapeutic avenues to prevent or overcome treatment resistance and improve patient outcomes.

## Analysis scripts

**HuH6 scRNA-seq analysis** 
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - Partitioning around medoids (PAM) clustering
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
