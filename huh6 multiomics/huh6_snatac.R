## HuH6 Signac Analysis

## /////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## [ Load Dependencies ] ----

library(Signac)
library(Seurat)
library(ChIPseeker)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(msigdbr)
library(TFBSTools)
library(RSQLite)
library(JASPAR2024)
library(motifmatchr)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(cowplot)
library(hues)
library(writexl)

plan("multicore", workers = 4)
options(future.globals.maxSize = 16000 * 1024^2)

## Colours for samples
require(RColorBrewer)
pal <- brewer.pal(12, "Paired")
huh6.cols <- c("Untreated_A" = pal[1],
               "Untreated_B" = pal[2],
               "Recovered_A" = pal[3],
               "Recovered_B" = pal[4])

## Colours for PAM clusters
require(paletteer)
pal <- paletteer_d("rcartocolor::Pastel")
pam.cols <- c("Hepatocytic 1" = pal[2], 
              "Hepatocytic 2" = pal[4], 
              "Hepatocytic 3" = pal[3], 
              "Late progenitor" = pal[5],
              "Early progenitor" = pal[1],
              "Stem-like" = pal[6])

## [ Preprocessing ] ----

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

## Untreated A
uta.counts <- Read10X_h5(filename = "signac_data/huh6_data/untreated_a_filtered_feature_bc_matrix.h5")
untreated.a <- CreateSeuratObject(
  counts = uta.counts$`Gene Expression`,
  assay = "RNA"
)

uta.fragpath <- "signac_data/huh6_data/untreated_a_atac_fragments.tsv.gz"
untreated.a[["ATAC"]] <- CreateChromatinAssay(
  counts = uta.counts$Peaks,
  sep = c(":", "-"),
  fragments = uta.fragpath,
  annotation = annotation
)
untreated.a

## Untreated B
utb.counts <- Read10X_h5(filename = "signac_data/huh6_data/untreated_b_filtered_feature_bc_matrix.h5")
untreated.b <- CreateSeuratObject(
  counts = utb.counts$`Gene Expression`,
  assay = "RNA"
)

utb.fragpath <- "signac_data/huh6_data/untreated_b_atac_fragments.tsv.gz"
untreated.b[["ATAC"]] <- CreateChromatinAssay(
  counts = utb.counts$Peaks,
  sep = c(":", "-"),
  fragments = utb.fragpath,
  annotation = annotation
)
untreated.b

## Recovered A
reca.counts <- Read10X_h5(filename = "signac_data/huh6_data/recovered_a_filtered_feature_bc_matrix.h5")
recovered.a <- CreateSeuratObject(
  counts = reca.counts$`Gene Expression`,
  assay = "RNA"
)

reca.fragpath <- "signac_data/huh6_data/recovered_a_atac_fragments.tsv.gz"
recovered.a[["ATAC"]] <- CreateChromatinAssay(
  counts = reca.counts$Peaks,
  sep = c(":", "-"),
  fragments = reca.fragpath,
  annotation = annotation
)
recovered.a

## Recovered B
recb.counts <- Read10X_h5(filename = "signac_data/huh6_data/recovered_b_filtered_feature_bc_matrix.h5")
recovered.b <- CreateSeuratObject(
  counts = recb.counts$`Gene Expression`,
  assay = "RNA"
)

recb.fragpath <- "signac_data/huh6_data/recovered_b_atac_fragments.tsv.gz"
recovered.b[["ATAC"]] <- CreateChromatinAssay(
  counts = recb.counts$Peaks,
  sep = c(":", "-"),
  fragments = recb.fragpath,
  annotation = annotation
)
recovered.b

## [ Quality Control ] ----

## Untreated A
DefaultAssay(untreated.a) <- "ATAC"
untreated.a <- NucleosomeSignal(untreated.a)
untreated.a <- TSSEnrichment(untreated.a)
saveRDS(untreated.a, "signac_data/huh6_data/untreated_a_qc.rds")

## Untreated B
DefaultAssay(untreated.b) <- "ATAC"
untreated.b <- NucleosomeSignal(untreated.b)
untreated.b <- TSSEnrichment(untreated.b)
saveRDS(untreated.b, "signac_data/huh6_data/untreated_b_qc.rds")

## Recovered A
DefaultAssay(recovered.a) <- "ATAC"
recovered.a <- NucleosomeSignal(recovered.a)
recovered.a <- TSSEnrichment(recovered.a)
saveRDS(recovered.a, "signac_data/huh6_data/recovered_a_qc.rds")

## Recovered B
DefaultAssay(recovered.b) <- "ATAC"
recovered.b <- NucleosomeSignal(recovered.b)
recovered.b <- TSSEnrichment(recovered.b)
saveRDS(recovered.b, "signac_data/huh6_data/recovered_b_qc.rds")

## /////////////////////////////////////////////////////////////////////////////
## QC plots ////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Untreated A
pdf("signac_plots/huh6_plots/untreated_a_tss_count.pdf", width = 5.8, height = 5.8)
png("signac_plots/huh6_plots/untreated_a_tss_count.png", units = "in", width = 5.8, height = 5.8, res = 1200)
DensityScatter(untreated.a, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
dev.off()

VlnPlot(object = untreated.a,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)
ggsave("signac_plots/huh6_plots/untreated_a_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/untreated_a_qc_plots.png", width = 8.3, height = 5.8)
dev.off()

## Untreated B
pdf("signac_plots/huh6_plots/untreated_b_tss_count.pdf", width = 5.8, height = 5.8)
png("signac_plots/huh6_plots/untreated_b_tss_count.png", units = "in", width = 5.8, height = 5.8, res = 1200)
DensityScatter(untreated.b, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
dev.off()

VlnPlot(object = untreated.b,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)
ggsave("signac_plots/huh6_plots/untreated_b_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/untreated_b_qc_plots.png", width = 8.3, height = 5.8)
dev.off()

## Recovered A
pdf("signac_plots/huh6_plots/recovered_a_tss_count.pdf", width = 5.8, height = 5.8)
png("signac_plots/huh6_plots/recovered_a_tss_count.png", units = "in", width = 5.8, height = 5.8, res = 1200)
DensityScatter(recovered.a, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
dev.off()

VlnPlot(object = recovered.a,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)
ggsave("signac_plots/huh6_plots/recovered_a_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/recovered_a_qc_plots.png", width = 8.3, height = 5.8)
dev.off()

## Recovered B
pdf("signac_plots/huh6_plots/recovered_b_tss_count.pdf", width = 5.8, height = 5.8)
png("signac_plots/huh6_plots/recovered_b_tss_count.png", units = "in", width = 5.8, height = 5.8, res = 1200)
DensityScatter(recovered.b, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
dev.off()

VlnPlot(object = recovered.b,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)
ggsave("signac_plots/huh6_plots/recovered_b_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/recovered_b_qc_plots.png", width = 8.3, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## QC filtering ////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

dim(untreated.a) ## 185031   6075
## Filter out low quality cells
untreated.a <- subset(x = untreated.a,
                      subset = nCount_ATAC < 100000 &
                        nCount_RNA < 50000 &
                        nCount_ATAC > 500 &
                        nCount_RNA > 500 &
                        nucleosome_signal < 2 &
                        TSS.enrichment > 1)
dim(untreated.a) ## 185031   5289
saveRDS(untreated.a, "signac_data/huh6_data/untreated_a_filt.rds")

dim(untreated.b) ## 227634   6365
## Filter out low quality cells
untreated.b <- subset(x = untreated.b,
                      subset = nCount_ATAC < 250000 &
                        nCount_RNA < 50000 &
                        nCount_ATAC > 1000 &
                        nCount_RNA > 1000 &
                        nucleosome_signal < 2 &
                        TSS.enrichment > 1)
dim(untreated.b) ## 227634   5308
saveRDS(untreated.b, "signac_data/huh6_data/untreated_b_filt.rds")

dim(recovered.a) ## 124638    562
## Filter out low quality cells
recovered.a <- subset(x = recovered.a,
                      subset = nCount_ATAC < 200000 &
                        nCount_RNA < 50000 &
                        nCount_ATAC > 500 &
                        nCount_RNA > 500 &
                        nucleosome_signal < 2 &
                        TSS.enrichment > 1)
dim(recovered.a) ## 124638    501
saveRDS(recovered.a, "signac_data/huh6_data/recovered_a_filt.rds")

dim(recovered.b) ## 184063   4746
## Filter out low quality cells
recovered.b <- subset(x = recovered.b,
                      subset = nCount_ATAC < 200000 &
                        nCount_RNA < 50000 &
                        nCount_ATAC > 1000 &
                        nCount_RNA > 1000 &
                        nucleosome_signal < 2 &
                        TSS.enrichment > 1)
dim(recovered.b) ## 184063   3736
saveRDS(recovered.b, "signac_data/huh6_data/recovered_b_filt.rds")

## [ Merge Objects ] ----

untreated.a <- readRDS("signac_data/huh6_data/untreated_a_filt.rds")
untreated.b <- readRDS("signac_data/huh6_data/untreated_b_filt.rds")
recovered.a <- readRDS("signac_data/huh6_data/recovered_a_filt.rds")
recovered.b <- readRDS("signac_data/huh6_data/recovered_b_filt.rds")

## /////////////////////////////////////////////////////////////////////////////
## Peak calling with MACS //////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

macs3.path = "/Users/echen/MACS3env/bin/macs3"

custom <- c("CD34", "PTPRC", "MCAM", ## Generally not expressed in normal liver, high in HCC; PTPRC = CD45 expressed in HPCs
            "FOXA1", "FOXA2", "GATA4", ## Early hepatic specification
            "EPCAM", "NCAM1", "CLDN3", "PROM1", "SOX17", ## HSC markers; EPCAM also HB marker
            "KRT8", "KRT18", "KRT19", ## HSC markers
            "CXCR4", ## Endoderm marker expressed in HSCs and HCC
            "AFP", "KRT7", "ICAM1", ## HB markers; AFP also fetal hepatocyte; KRT7 (CK7) also cholangiocyte marker
            "HNF4A", ## HSC-HB and HB-hepatocyte differentiation?
            "ONECUT2", ## HB migration
            "PROX1", "TBX3", ## HB proliferation and migration
            "SOX9", ## Hepatic progenitor and cholangiocyte
            "HNF1B", "SALL4", ## Cholangiocyte fate regulator
            "CEBPA", "ALB", "TTR", "RBPJ", "NR5A2", ## Hepatocytic cell fate
            "MAT1A", "NR1I2", ## Adult liver
            "MAT2A", ## Fetal liver and replaces MAT1A in HCC
            "TAF10", "TBP", ## Embryonic liver not adult
            "TAF4", ## Postnatal hepatocytes
            "TGFB1", "TGFB2", "TGFBR2") ## High TGFb = cholangiocyte differentiation

## Untreated A
uta.peaks <- CallPeaks(object = untreated.a, macs2.path = macs3.path)
## Remove peaks on nonstandard chromosomes and in genomic blacklist regions
uta.peaks <- keepStandardChromosomes(uta.peaks, pruning.mode = "coarse")
uta.peaks <- subsetByOverlaps(x = uta.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
## Remove sex chromosomes (X and Y)
chroms.to.keep <- which(seqnames(uta.peaks) %in% paste0("chr", 1:22))
uta.peaks <- uta.peaks[chroms.to.keep, ]
## Visualise coverage on genes of interest
pdf("signac_plots/huh6_plots/ut_a_custom_coverage_plots.pdf")
for (gene in custom) {
  print(CoveragePlot(object = untreated.a,
                     region = gene,
                     ranges = uta.peaks,
                     ranges.title = "MACS2"))
}
dev.off()
saveRDS(uta.peaks, "signac_data/huh6_data/untreated_a_peak_granges.rds")

## Untreated B
utb.peaks <- CallPeaks(object = untreated.b, macs2.path = macs3.path)
utb.peaks <- keepStandardChromosomes(utb.peaks, pruning.mode = "coarse")
utb.peaks <- subsetByOverlaps(x = utb.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
## Remove sex chromosomes (X and Y)
chroms.to.keep <- which(seqnames(utb.peaks) %in% paste0("chr", 1:22))
utb.peaks <- utb.peaks[chroms.to.keep, ]
## Visualise coverage on genes of interest
pdf("signac_plots/huh6_plots/ut_b_custom_coverage_plots.pdf")
for (gene in custom) {
  print(CoveragePlot(object = untreated.b,
                     region = gene,
                     ranges = utb.peaks,
                     ranges.title = "MACS2"))
}
dev.off()
saveRDS(utb.peaks, "signac_data/huh6_data/untreated_b_peak_granges.rds")

## Recovered A
reca.peaks <- CallPeaks(object = recovered.a, macs2.path = macs3.path)
reca.peaks <- keepStandardChromosomes(reca.peaks, pruning.mode = "coarse")
reca.peaks <- subsetByOverlaps(x = reca.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
## Remove sex chromosomes (X and Y)
chroms.to.keep <- which(seqnames(reca.peaks) %in% paste0("chr", 1:22))
reca.peaks <- reca.peaks[chroms.to.keep, ]
## Visualise coverage on genes of interest
pdf("signac_plots/huh6_plots/rec_a_custom_coverage_plots.pdf")
for (gene in custom) {
  print(CoveragePlot(object = recovered.a,
                     region = gene,
                     ranges = reca.peaks,
                     ranges.title = "MACS2"))
}
dev.off()
saveRDS(reca.peaks, "signac_data/huh6_data/recovered_a_peak_granges.rds")

## Recovered B
recb.peaks <- CallPeaks(object = recovered.b, macs2.path = macs3.path)
recb.peaks <- keepStandardChromosomes(recb.peaks, pruning.mode = "coarse")
recb.peaks <- subsetByOverlaps(x = recb.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
## Remove sex chromosomes (X and Y)
chroms.to.keep <- which(seqnames(recb.peaks) %in% paste0("chr", 1:22))
recb.peaks <- recb.peaks[chroms.to.keep, ]
## Visualise coverage on genes of interest
pdf("signac_plots/huh6_plots/rec_b_custom_coverage_plots.pdf")
for (gene in custom) {
  print(CoveragePlot(object = recovered.b,
                     region = gene,
                     ranges = recb.peaks,
                     ranges.title = "MACS2"))
}
dev.off()
saveRDS(recb.peaks, "signac_data/huh6_data/recovered_b_peak_granges.rds")

## /////////////////////////////////////////////////////////////////////////////
## Merge objects ///////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Create a common peak set
uta.peaks <- readRDS("signac_data/huh6_data/untreated_a_peak_granges.rds")
utb.peaks <- readRDS("signac_data/huh6_data/untreated_b_peak_granges.rds")
reca.peaks <- readRDS("signac_data/huh6_data/recovered_a_peak_granges.rds")
recb.peaks <- readRDS("signac_data/huh6_data/recovered_b_peak_granges.rds")

## Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(uta.peaks, utb.peaks, reca.peaks, recb.peaks))
## Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
min(peakwidths) ## 200
max(peakwidths) ## 3848
# combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
saveRDS(combined.peaks, "signac_data/huh6_data/combined_peaks.rds")

## Prepare metadata and fragment objects
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
## Untreated A
uta.meta <- untreated.a@meta.data
uta.frags <- Fragments(untreated.a)[[1]]
uta.counts <- FeatureMatrix(fragments = uta.frags,
                            features = combined.peaks,
                            cells = rownames(uta.meta))
uta.assay <- CreateChromatinAssay(uta.counts, fragments = uta.frags, annotation = annotation)
uta.seurat <- CreateSeuratObject(uta.assay, assay = "ATAC", meta.data = uta.meta)

## Untreated B
utb.meta <- untreated.b@meta.data
utb.frags <- Fragments(untreated.b)[[1]]
utb.counts <- FeatureMatrix(fragments = utb.frags,
                            features = combined.peaks,
                            cells = rownames(utb.meta))
utb.assay <- CreateChromatinAssay(utb.counts, fragments = utb.frags, annotation = annotation)
utb.seurat <- CreateSeuratObject(utb.assay, assay = "ATAC", meta.data = utb.meta)

## Recovered A
reca.meta <- recovered.a@meta.data
reca.frags <- Fragments(recovered.a)[[1]]
reca.counts <- FeatureMatrix(fragments = reca.frags,
                             features = combined.peaks,
                             cells = rownames(reca.meta))
reca.assay <- CreateChromatinAssay(reca.counts, fragments = reca.frags, annotation = annotation)
reca.seurat <- CreateSeuratObject(reca.assay, assay = "ATAC", meta.data = reca.meta)

## Recovered B
recb.meta <- recovered.b@meta.data
recb.frags <- Fragments(recovered.b)[[1]]
recb.counts <- FeatureMatrix(fragments = recb.frags,
                             features = combined.peaks,
                             cells = rownames(recb.meta))
recb.assay <- CreateChromatinAssay(recb.counts, fragments = recb.frags, annotation = annotation)
recb.seurat <- CreateSeuratObject(recb.assay, assay = "ATAC", meta.data = recb.meta)

## Add metadata column to identify dataset of origin
uta.seurat$Sample <- "Untreated_A"
utb.seurat$Sample <- "Untreated_B"
reca.seurat$Sample <- "Recovered_A"
recb.seurat$Sample <- "Recovered_B"

combined <- merge(x = uta.seurat,
                  y = c(utb.seurat, reca.seurat, recb.seurat),
                  add.cell.ids = c("untreated_a", "untreated_b",
                                   "recovered_a", "recovered_b"))
saveRDS(combined, "signac_data/huh6_data/combined_seurat.rds")

## Normalisation and dimensionality reduction
combined <- readRDS("signac_data/huh6_data/combined_seurat.rds")

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined)
combined <- RunSVD(combined)
## Don't use LSI dimension 1 for downstream as is correlated with sequencing depth
# FeatureScatter(combined, "nCount_ATAC", "LSI_1")
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = "Sample", pt.size = 0.1) +
  scale_colour_manual(values = huh6.cols) +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_umap.pdf", width = 5.8, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap.png", width = 5.8, height = 5.8)

combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
cluster.search <- function(seurat, from = 0.2, to = 0.8, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("ATAC_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}
combined <- cluster.search(combined)

DimPlot(object = combined, label = TRUE, group.by = "seurat_clusters.0.2") + NoLegend() +
  scale_colour_manual(values = iwanthue(length(unique(combined$seurat_clusters.0.2)))) +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_umap_clusters.pdf", width = 5.8, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap_clusters.png", width = 5.8, height = 5.8)

saveRDS(combined, "signac_data/huh6_data/combined_seurat_norm.rds")

## [ Gene Activity ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_norm.rds")

gene.activities <- GeneActivity(combined)

combined[["Activity"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(combined) <- "Activity"
combined <- NormalizeData(combined)
combined <- ScaleData(combined, features = rownames(combined))

saveRDS(combined, "signac_data/huh6_data/combined_seurat_activity.rds")

## Use scRNA-seq data to annotate all cells
reference <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/huh6_seurat_hvgs_pam_anno.rds")

transfer.anchors <- FindTransferAnchors(reference = reference,
                                        query = combined,
                                        features = VariableFeatures(object = reference),
                                        normalization.method = "LogNormalize",
                                        reference.assay = "RNA",
                                        query.assay = "Activity",
                                        reduction = "cca",
                                        dims = 2:50)

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$PAM_Name,
                                     weight.reduction = combined[["lsi"]], dims = 2:50)
saveRDS(celltype.predictions, "signac_data/huh6_data/pam_predictions.rds")

combined <- AddMetaData(combined, metadata = celltype.predictions[, c("predicted.id", "prediction.score.max")])

## Use snRNA-seq data to annotate cells with both sets of data
reference <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/multiomics/hb/snRNA-seq/data/huh6_seurat_hvgs_pam_anno.rds")

combined$CellID <- rownames(combined@meta.data)

rna.meta <- reference@meta.data[, c("CellID", "PAM_Name")]
rna.meta$CellID <- gsub("Recovered_A", "recovered_a", rna.meta$CellID)
rna.meta$CellID <- gsub("Recovered_B", "recovered_b", rna.meta$CellID)
rna.meta$CellID <- gsub("Untreated_A", "untreated_a", rna.meta$CellID)
rna.meta$CellID <- gsub("Untreated_B", "untreated_b", rna.meta$CellID)
rownames(rna.meta) <- rna.meta$CellID
rna.meta <- rna.meta[combined$CellID, ]
combined <- AddMetaData(combined, metadata = rna.meta)

combined$annotation_correct <- combined$predicted.id == combined$PAM_Name

umap1 <- DimPlot(object = combined, group.by = "predicted.id") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("Predicted PAM Annotation") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)

umap2 <- DimPlot(object = combined, group.by = "PAM_Name") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("PAM Annotation") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)

umap1 | umap2
ggsave("signac_plots/huh6_plots/combined_umap_pam_annotation.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap_pam_annotation.png", width = 8.3, height = 5.8)

predictions <- table(combined$PAM_Name, combined$predicted.id)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)

predictions.df <- predictions

p1 <- ggplot(predictions.df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") +
  xlab("Cluster annotation (RNA)") + ylab("Predicted cluster (ATAC)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 7),
        text = element_text(size = 15),
        aspect.ratio = 1)

correct <- length(which(combined$PAM_Name == combined$predicted.id))
incorrect <- length(which(combined$PAM_Name != combined$predicted.id))
data <- FetchData(combined, vars = c("prediction.score.max", "annotation_correct"))

p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) +
  theme_cowplot() +
  theme(text = element_text(size = 15),
        aspect.ratio = 1) +
  scale_fill_discrete(name = "Annotation Correct",
                      labels = c(paste0("FALSE (n = ", incorrect, ")"),
                                 paste0("TRUE (n = ", correct, ")"))) +
  scale_color_discrete(name = "Annotation Correct",
                       labels = c(paste0("FALSE (n = ", incorrect, ")"),
                                  paste0("TRUE (n = ", correct, ")"))) +
  xlab("Prediction Score")

p1 + p2
ggsave("signac_plots/huh6_plots/combined_pam_annotation_prediction.pdf", width = 11.7, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_pam_annotation_prediction.png", width = 11.7, height = 5.8)

saveRDS(combined, "signac_data/huh6_data/combined_seurat_pam_annotated.rds")

## [ Joint UMAP ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_pam_annotated.rds")
reference <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/multiomics/hb/snRNA-seq/data/huh6_seurat_hvgs_pam_anno.rds")

dim(combined) ## 19620 14834
dim(reference) ## 12473 13308

colnames(reference) <- gsub("Recovered_A", "recovered_a", colnames(reference))
colnames(reference) <- gsub("Recovered_B", "recovered_b", colnames(reference))
colnames(reference) <- gsub("Untreated_A", "untreated_a", colnames(reference))
colnames(reference) <- gsub("Untreated_B", "untreated_b", colnames(reference))

length(colnames(reference)[colnames(reference) %in% colnames(combined)]) ## 12515
cells.to.keep <- intersect(colnames(reference), colnames(combined))
combined <- subset(combined, cells = cells.to.keep)
reference <- subset(reference, cells = cells.to.keep)
dim(combined) ## 19620 12515
dim(reference) ## 12473 12515

combined[["RNA"]] <- reference[["RNA"]]

DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined)
combined <- RunPCA(combined)

combined <- FindMultiModalNeighbors(object = combined,
                                    reduction.list = list("pca", "lsi"), 
                                    dims.list = list(1:50, 2:50),
                                    modality.weight.name = "RNA.weight",
                                    verbose = TRUE)

combined <- RunUMAP(object = combined,
                    nn.name = "weighted.nn",
                    assay = "RNA",
                    verbose = TRUE)

DimPlot(object = combined, group.by = "Sample",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = huh6.cols) +
  ggtitle("Weighted Nearest Neighbour") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_umap_wnn.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap_wnn.png", width = 8.3, height = 5.8)

DimPlot(object = combined, group.by = "Sample",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = huh6.cols) +
  ggtitle("Weighted Nearest Neighbour") &
  theme(text = element_text(size = 15), aspect.ratio = 1,
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = NA))
ggsave("signac_plots/huh6_plots/combined_umap_wnn_bg.pdf", width = 8.3, height = 5.8, bg = "transparent")
ggsave("signac_plots/huh6_plots/combined_umap_wnn_bg.png", width = 8.3, height = 5.8, bg = "transparent")

DimPlot(object = combined, group.by = "PAM_Name",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("Weighted Nearest Neighbour") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_umap_wnn_pam_cluster.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap_wnn_pam_cluster.png", width = 8.3, height = 5.8)

DimPlot(object = combined, group.by = "predicted.id",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("Weighted Nearest Neighbour") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_umap_wnn_predicted_pam.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_umap_wnn_predicted_pam.png", width = 8.3, height = 5.8)

saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn.rds")

## [ Rerun UMAP ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")
DefaultAssay(combined) <- "ATAC"

## scATAC-seq data UMAPs after filtering for cells that are also in the snRNA-seq data
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
combined <- FindClusters(combined, resolution = 0.2)

DimPlot(object = combined, group.by = "Sample",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = huh6.cols) +
  ggtitle("Filtered scATAC-seq UMAP") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_filt_umap.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_filt_umap.png", width = 8.3, height = 5.8)

DimPlot(object = combined, group.by = "PAM_Name",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("Filtered scATAC-seq UMAP") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_filt_umap_pam_cluster.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_filt_umap_pam_cluster.png", width = 8.3, height = 5.8)

DimPlot(object = combined, group.by = "predicted.id",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = pam.cols) +
  ggtitle("Filtered scATAC-seq UMAP") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_filt_umap_predicted_pam.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_filt_umap_predicted_pam.png", width = 8.3, height = 5.8)

DimPlot(object = combined, group.by = "seurat_clusters.0.2",
        repel = TRUE, reduction = "umap") +
  scale_colour_manual(values = iwanthue(length(unique(combined$seurat_clusters.0.2)))) +
  ggtitle("Filtered scATAC-seq UMAP") +
  theme(text = element_text(size = 15),
        aspect.ratio = 1)
ggsave("signac_plots/huh6_plots/combined_filt_umap_seurat_cluster.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/combined_filt_umap_seurat_cluster.png", width = 8.3, height = 5.8)

saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_umap.rds")

## [ Differential Accessibility ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")

## /////////////////////////////////////////////////////////////////////////////
## Between conditions //////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Find differentially accessible peaks between experimental conditions
DefaultAssay(combined) <- "ATAC"
combined$Condition <- ifelse(combined$Sample %in% c("Untreated_A", "Untreated_B"), "Untreated", "Recovered")
Idents(combined) <- "Condition"
saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn.rds")

da.peaks <- FindMarkers(object = combined,
                        ident.1 = "Recovered",
                        ident.2 = "Untreated",
                        test.use = "wilcox")
head(da.peaks)
##                                  p_val avg_log2FC pct.1 pct.2     p_val_adj
## chr19-49363047-49364075   0.000000e+00   2.329838 0.358 0.081  0.000000e+00
## chr20-401476-402633      1.297009e-294   1.978032 0.399 0.119 2.881383e-289
## chr19-3284554-3286213    1.612156e-279   1.537275 0.511 0.210 3.581502e-274
## chr11-77299451-77300388  3.001028e-273   3.937448 0.177 0.012 6.666963e-268
## chr9-137149893-137151576 3.351422e-243   1.124679 0.636 0.341 7.445384e-238
## chr17-74416842-74418025  6.255167e-239   2.458304 0.252 0.051 1.389623e-233

saveRDS(da.peaks, "signac_data/huh6_data/recovered_vs_untreated_markers.rds")

peaks.filtered <- dplyr::filter(da.peaks, abs(avg_log2FC) > 2, p_val_adj < 0.05)
peaks.filtered <- peaks.filtered[order(peaks.filtered$avg_log2FC, decreasing = TRUE), ]
head(peaks.filtered)
##                                  p_val avg_log2FC pct.1 pct.2    p_val_adj
## chr1-61139300-61139672   3.952965e-32   4.563557 0.019 0.001 8.781750e-27
## chr2-142944126-142944537 2.198304e-18   4.461684 0.010 0.000 4.883665e-13
## chr1-61265690-61266058   6.626718e-22   4.238477 0.013 0.001 1.472165e-16
## chr8-101635744-101636024 2.708882e-40   4.183626 0.025 0.001 6.017944e-35
## chr6-21005765-21006103   1.136114e-31   4.180352 0.020 0.001 2.523946e-26
## chr8-16825982-16826328   3.533592e-37   4.179760 0.023 0.001 7.850087e-32

fc <- FoldChange(combined, ident.1 = "Recovered", ident.2 = "Untreated")
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)
##                          avg_log2FC pct.1 pct.2
## chr3-186860810-186861040   8.057358 0.005 0.000
## chr3-98391321-98391523     6.252327 0.001 0.000
## chr22-43707647-43707927    5.107043 0.009 0.000
## chr17-52673769-52674082    4.889194 0.007 0.000
## chr13-66657175-66657394    4.612757 0.006 0.000
## chr1-61139300-61139672     4.563557 0.019 0.001

open.rec <- rownames(peaks.filtered[peaks.filtered$avg_log2FC > 3, ])
open.ut <- rownames(peaks.filtered[peaks.filtered$avg_log2FC < -3, ])

closest.genes.rec <- ClosestFeature(combined, regions = open.rec)
closest.genes.ut <- ClosestFeature(combined, regions = open.ut)

head(closest.genes.rec)
##            tx_id    gene_name         gene_id   gene_biotype type           closest_region             query_region distance
## ENST00000485903   ENST00000485903         NFIA ENSG00000162599 protein_coding  gap   chr1-61088681-61277519   chr1-61139300-61139672        0
## ENST00000621320   ENST00000621320         KYNU ENSG00000115919 protein_coding  gap chr2-142927742-142954809 chr2-142944126-142944537        0
## ENST00000485903.1 ENST00000485903         NFIA ENSG00000162599 protein_coding  gap   chr1-61088681-61277519   chr1-61265690-61266058        0
## ENST00000517674   ENST00000517674        GRHL2 ENSG00000083307 protein_coding  gap chr8-101632366-101636896 chr8-101635744-101636024        0
## ENST00000378610   ENST00000378610       CDKAL1 ENSG00000145996 protein_coding  gap   chr6-21000373-21065047   chr6-21005765-21006103        0
## ENST00000521411   ENST00000521411 RP11-13N12.1 ENSG00000253496        lincRNA  gap   chr8-16783283-16844104   chr8-16825982-16826328        0
head(closest.genes.ut)
##            tx_id     gene_name         gene_id   gene_biotype type            closest_region              query_region distance
## ENSE00002197186 ENST00000529689       CCDC90B ENSG00000137500 protein_coding exon   chr11-83285873-83286407   chr11-83318644-83319407    32236
## ENST00000369796 ENST00000369796        STRIP1 ENSG00000143093 protein_coding  cds  chr1-110045015-110045078  chr1-110045018-110045495        0
## ENST00000635529 ENST00000635529       CSNK1G1 ENSG00000169118 protein_coding  utr   chr15-64355988-64356132   chr15-64355426-64356679        0
## ENST00000410059 ENST00000410059         DPP10 ENSG00000175497 protein_coding  utr  chr2-115842346-115845752  chr2-116756477-116756995   910724
## ENSE00002206814 ENST00000306533 RP11-136I14.5 ENSG00000255689        lincRNA exon chr11-115600221-115600339 chr11-115610711-115611096    10371
## ENST00000437570 ENST00000437570         RUFY1 ENSG00000176783 protein_coding  utr  chr5-179559692-179560038  chr5-179558918-179560016        0

peaks.filtered$query_region <- rownames(peaks.filtered)

rec.results <- merge(closest.genes.rec, peaks.filtered, by = "query_region")
ut.results <- merge(closest.genes.ut, peaks.filtered, by = "query_region")

write.csv(rec.results, "signac_data/huh6_data/recovery_da_results.csv", row.names = FALSE)
write.csv(ut.results, "signac_data/huh6_data/untreated_da_results.csv", row.names = FALSE)

## /////////////////////////////////////////////////////////////////////////////
## Between clusters ////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Find differentially accessible peaks between PAM clusters
Idents(combined) <- combined$PAM_Name

pam.markers <- FindAllMarkers(combined,
                              assay = "ATAC",
                              slot = "data",
                              only.pos = TRUE,
                              densify = TRUE,
                              min.pct = 0.05,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")

table(pam.markers$cluster)
## Hepatocytic_1 Progenitor-cholangiocytic    Progenitor-hepatocytic          Progenitor TGFb+          Hepatoblast-like 
##        14932                         5                     21360                     35488                      5019 

## Couldn't run ClosestFeature() with rownames(pam.markers)
## Some regions were outside the valid range
open.pam <- rownames(pam.markers)
split.regions <- strsplit(open.pam, "-")
chrom <- sapply(split.regions, `[`, 1)
start <- as.numeric(sapply(split.regions, `[`, 2))
end <- as.numeric(sapply(split.regions, `[`, 3))

summary(start)
summary(end)
invalid.start <- which(start > 2^31 | start < -2^31)
invalid.end <- which(end > 2^31 | end < -2^31)
open.pam[invalid.start]
open.pam[invalid.end]

valid.indices <- which(start <= 2^31 & start >= -2^31 & end <= 2^31 & end >= -2^31)
start <- start[valid.indices]
end <- end[valid.indices]
chrom <- chrom[valid.indices]

open.pam.gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
closest.genes <- ClosestFeature(combined, regions = open.pam.gr)

pam.markers$query_region <- pam.markers$gene
pam.markers$gene <- NULL

markers.df <- merge(pam.markers, closest.genes, by = "query_region")
pam.marker.list <- lapply(unique(markers.df$cluster), function(x){
  tmp.df <- markers.df[markers.df$cluster==x & markers.df$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})
names(pam.marker.list) <- unique(markers.df$cluster)
saveRDS(pam.marker.list, "signac_data/huh6_data/pam_markers.rds")

## GSEA on PAM cluster markers
wilcox.GO.BP.list <- lapply(pam.marker.list, function(x){
  enrichGO(x$gene_name,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
wilcox.C2.list <- lapply(pam.marker.list, function(x){
  enricher(x$gene_name,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
wilcox.H.list <- lapply(pam.marker.list, function(x){
  enricher(x$gene_name,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(wilcox.GO.BP.list, "signac_data/huh6_data/pam_markers_BP_ontology.rds")
saveRDS(wilcox.C2.list, "signac_data/huh6_data/pam_markers_C2_geneset.rds")
saveRDS(wilcox.H.list, "signac_data/huh6_data/pam_markers_Hallmarks.rds")

wilcox.lists <- list("wilcox.GO.BP.list" = wilcox.GO.BP.list,
                     "wilcox.C2.list" = wilcox.C2.list,
                     "wilcox.H.list" = wilcox.H.list,
                     "pam.list" = pam.marker.list)

wilcox.names <- sub("list", "sheets", names(wilcox.lists))

for(i in 1:length(wilcox.lists)){
  wilcox.list <- wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(wilcox.names[[i]], list.sheets)
}

write_xlsx(wilcox.GO.BP.sheets, "signac_data/huh6_data/pam_markers_BP_ontology.xlsx")
write_xlsx(wilcox.C2.sheets, "signac_data/huh6_data/pam_markers_C2_geneset.xlsx")
write_xlsx(wilcox.H.sheets, "signac_data/huh6_data/pam_markers_Hallmarks.xlsx")
write_xlsx(pam.sheets, "signac_data/huh6_data/pam_markers.xlsx")

## /////////////////////////////////////////////////////////////////////////////
## Pseudotime target genes /////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

combined <- readRDS("signac_data/huh6_data/combined_seurat_pam_annotated.rds")
DefaultAssay(combined) <- "ATAC"
peak.ranges <- StringToGRanges(rownames(combined), sep = c("-", "-"))
gene.anno <- Annotation(combined)
target.genes <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/trajectory_genes/test_genes.rds")
target.anno <- gene.anno[gene.anno$gene_name %in% target.genes]

overlaps <- findOverlaps(peak.ranges, target.anno)
target.peaks <- unique(rownames(combined)[queryHits(overlaps)])

combined$Condition <- ifelse(combined$Sample %in% c("Untreated_A", "Untreated_B"), "Untreated", "Recovered")
table(combined$Sample, combined$Condition)
##             Recovered Untreated
## Recovered_A       501         0
## Recovered_B      3736         0
## Untreated_A         0      5289
## Untreated_B         0      5308
Idents(combined) <- "Condition"

da.peaks <- FindMarkers(object = combined,
                        ident.1 = "Recovered",
                        ident.2 = "Untreated",
                        features = target.peaks,
                        test.use = "wilcox")
head(da.peaks)
saveRDS(da.peaks, "signac_data/huh6_data/recovered_vs_untreated_pseudotime_genes.rds")

da.ranges <- StringToGRanges(rownames(da.peaks), sep = c("-", "-"))
da.overlaps <- findOverlaps(da.ranges, target.anno)
peak.genes <- data.frame(peak = rownames(da.peaks)[queryHits(da.overlaps)],
                         gene = target.anno$gene_name[subjectHits(da.overlaps)])
multi.genes <- peak.genes %>%
  group_by(peak) %>%
  summarise(genes = paste(unique(gene), collapse = "; ")) %>%
  ungroup()
da.anno <- da.peaks %>%
  mutate(peak = rownames(.)) %>%
  left_join(multi.genes, by = "peak")
saveRDS(da.anno, "signac_data/huh6_data/recovered_vs_untreated_pseudotime_genes_anno.rds")

da.filt <- dplyr::filter(da.anno, avg_log2FC > 0, p_val_adj < 0.05)
da.filt <- da.filt[order(da.filt$avg_log2FC, decreasing = TRUE), ]
write.csv(da.filt, "signac_data/huh6_data/recovered_pseudotime_genes_filtered.csv", row.names = FALSE)

test.genes <- unique(da.filt$genes)
saveRDS(test.genes, "signac_data/huh6_data/pseudotime_targets.rds")

trrust.df <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/trajectory_genes/trrust/trrust_ref.rds")
trrust.filt <- trrust.df[trrust.df$Target %in% test.genes, ]
write.csv(trrust.filt, "signac_data/huh6_data/recovered_pseudotime_genes_filtered_annotated.csv", row.names = FALSE)

## Volcano plots
da.anno <- readRDS("signac_data/huh6_data/recovered_vs_untreated_pseudotime_genes_anno.rds")
da.anno$neg_log10_pval <- -log10(da.anno$p_val_adj)
## Define colour groups
da.anno$colour <- "NS"
da.anno$colour[da.anno$avg_log2FC > 0.5 & da.anno$p_val_adj < 0.05] <- "Up"
da.anno$colour[da.anno$avg_log2FC < -0.5 & da.anno$p_val_adj < 0.05] <- "Down"
## Label only a subset of peaks
da.anno$label <- ifelse(da.anno$avg_log2FC > 0.5 & da.anno$p_val_adj < 0.05, da.anno$genes, "")

ggplot(da.anno, aes(x = avg_log2FC, y = neg_log10_pval, colour = colour)) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.5, colour = "black", max.overlaps = 100) +
  scale_colour_manual(values = c("Up" = "#B1283AFF",
                                 "Down" = "#006A8EFF",
                                 "NS" = "#A8A6A7FF"),
                      labels = c("Up" = "Log2FC > 0.5 and p-value < 0.05",
                                 "Down" = "Log2FC < 0.5 and p-value < 0.05",
                                 "NS" = "Not significant"),
                      name = NULL) +
  theme_bw() +
  labs(x = "Log2 Fold Change",
       y = "-Log10 P-value",
       title = "Differential accessibility of filtered plasticity targets",
       colour = "Gene Group") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  theme(legend.position = "bottom",
        text = element_text(size = 12))
ggsave("signac_plots/huh6_plots/differential_accesibility_volcano.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/differential_accesibility_volcano.png", width = 8.3, height = 5.8)

## Filtering based on plot
da.filt <- dplyr::filter(da.anno, avg_log2FC > 0.5, p_val_adj < 0.05)
da.filt <- da.filt[order(da.filt$avg_log2FC, decreasing = TRUE), ]
test.genes <- unique(da.filt$genes)
saveRDS(test.genes, "signac_data/huh6_data/pseudotime_targets_filtered.rds")

## BIRC5 track plots
combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")
DefaultAssay(combined) <- "ATAC"
combined$Condition <- ifelse(combined$Sample %in% c("Untreated_A", "Untreated_B"), "Untreated", "Recovered")
table(combined$Sample, combined$Condition)
Idents(combined) <- "Condition"
CoveragePlot(object = combined,
             region = "BIRC5",
             extend.upstream = 10000,
             extend.downstream = 10000) &
  scale_fill_manual(values = c("Untreated" = "#4477AA",
                               "Recovered" = "#66CCEE")) &
  labs(title = "BIRC5") &
  theme(text = element_text(size = 15))
ggsave("signac_plots/huh6_plots/coverage_plots_recovered_birc5_condition.pdf", height = 5.8, width = 8.3)
ggsave("signac_plots/huh6_plots/coverage_plots_recovered_birc5_condition.png", height = 5.8, width = 8.3)

Idents(combined) <- "PAM_Name"
CoveragePlot(object = combined,
             region = "BIRC5",
             extend.upstream = 10000,
             extend.downstream = 10000,
             split.assays = FALSE) &
  scale_fill_manual(values = pam.cols) &
  labs(title = "BIRC5") &
  theme(text = element_text(size = 15))
ggsave("signac_plots/huh6_plots/coverage_plots_recovered_birc5_cluster.pdf", height = 5.8, width = 8.3)
ggsave("signac_plots/huh6_plots/coverage_plots_recovered_birc5_cluster.png", height = 5.8, width = 8.3)

## /////////////////////////////////////////////////////////////////////////////
## HB target genes /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

combined <- readRDS("signac_data/huh6_data/combined_seurat_pam_annotated.rds")
DefaultAssay(combined) <- "ATAC"
peak.ranges <- StringToGRanges(rownames(combined), sep = c("-", "-"))
gene.anno <- Annotation(combined)
hb.genes <- c("CTNNB1", "LEF1", "HNF4A", "TCF4", "AFP")
target.anno <- gene.anno[gene.anno$gene_name %in% hb.genes]

overlaps <- findOverlaps(peak.ranges, target.anno)
target.peaks <- unique(rownames(combined)[queryHits(overlaps)])

combined$Condition <- ifelse(combined$Sample %in% c("Untreated_A", "Untreated_B"), "Untreated", "Recovered")
table(combined$Sample, combined$Condition)
##             Recovered Untreated
## Recovered_A       501         0
## Recovered_B      3736         0
## Untreated_A         0      5289
## Untreated_B         0      5308
Idents(combined) <- "Condition"

da.peaks <- FindMarkers(object = combined,
                        ident.1 = "Recovered",
                        ident.2 = "Untreated",
                        features = target.peaks,
                        test.use = "wilcox")
head(da.peaks)
saveRDS(da.peaks, "signac_data/huh6_data/recovered_vs_untreated_hb_genes.rds")

da.ranges <- StringToGRanges(rownames(da.peaks), sep = c("-", "-"))
da.overlaps <- findOverlaps(da.ranges, target.anno)
peak.genes <- data.frame(peak = rownames(da.peaks)[queryHits(da.overlaps)],
                         gene = target.anno$gene_name[subjectHits(da.overlaps)])
multi.genes <- peak.genes %>%
  group_by(peak) %>%
  summarise(genes = paste(unique(gene), collapse = "; ")) %>%
  ungroup()
da.anno <- da.peaks %>%
  mutate(peak = rownames(.)) %>%
  left_join(multi.genes, by = "peak")
saveRDS(da.anno, "signac_data/huh6_data/recovered_vs_untreated_hb_genes_anno.rds")

da.filt <- dplyr::filter(da.anno, p_val_adj < 0.05)
da.filt <- da.filt[order(da.filt$avg_log2FC, decreasing = TRUE), ]
write.csv(da.filt, "signac_data/huh6_data/recovered_hb_genes_filtered.csv", row.names = FALSE)
write.csv(da.filt[-grep("HNF4A|LEF1", da.filt$genes),],
          "signac_data/huh6_data/recovered_hb_genes_final.csv", row.names = FALSE)

## Volcano plots
da.anno <- readRDS("signac_data/huh6_data/recovered_vs_untreated_hb_genes_anno.rds")
da.anno$neg_log10_pval <- -log10(da.anno$p_val_adj)
## Define colour groups
da.anno$colour <- "NS"
da.anno$colour[da.anno$avg_log2FC > 0.5 & da.anno$p_val_adj < 0.05] <- "Up"
da.anno$colour[da.anno$avg_log2FC < -0.5 & da.anno$p_val_adj < 0.05] <- "Down"
## Label only a subset of peaks
da.anno$label <- ifelse(da.anno$avg_log2FC > 0.5 & da.anno$p_val_adj < 0.05, da.anno$genes, "")

ggplot(da.anno, aes(x = avg_log2FC, y = neg_log10_pval, colour = colour)) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.5, colour = "black", max.overlaps = 100) +
  scale_colour_manual(values = c("Up" = "#B1283AFF",
                                 "Down" = "#006A8EFF",
                                 "NS" = "#A8A6A7FF"),
                      labels = c("Up" = "Log2FC > 0.5 and p-value < 0.05",
                                 "Down" = "Log2FC < 0.5 and p-value < 0.05",
                                 "NS" = "Not significant"),
                      name = NULL) +
  theme_bw() +
  labs(x = "Log2 Fold Change",
       y = "-Log10 P-value",
       title = "Differential accessibility of hepatoblastoma markers",
       colour = "Gene Group") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  theme(legend.position = "bottom",
        text = element_text(size = 12))
ggsave("signac_plots/huh6_plots/differential_accesibility_volcano_hb_genes.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/differential_accesibility_volcano_hb_genes.png", width = 8.3, height = 5.8)

## [ Plot Genomic Regions ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")

## Plot genomic regions for cells grouped by condition
da.peaks <- readRDS("signac_data/huh6_data/recovered_vs_untreated_markers.rds")
CoveragePlot(object = combined,
  region = rownames(da.peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000)
ggsave("signac_plots/huh6_plots/coverage_plot_example.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/coverage_plot_example.png", width = 8.3, height = 5.8)

peaks.filtered <- dplyr::filter(da.peaks, abs(avg_log2FC) > 2, p_val_adj < 0.05)
peaks.filtered <- peaks.filtered[order(peaks.filtered$avg_log2FC, decreasing = TRUE), ]
fc <- FoldChange(combined, ident.1 = "Recovered", ident.2 = "Untreated")
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]

CoveragePlot(object = combined,
             region = rownames(fc)[1],
             extend.upstream = 40000,
             extend.downstream = 20000)
ggsave("signac_plots/huh6_plots/coverage_plot_example_fc.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/coverage_plot_example_fc.png", width = 8.3, height = 5.8)

## Plot DNA accessibility in regions for markers used to annotate snRNA-seq clusters
custom <- c("CD34", "PTPRC", "MCAM", ## Generally not expressed in normal liver, high in HCC; PTPRC = CD45 expressed in HPCs
            "FOXA1", "FOXA2", "GATA4", ## Early hepatic specification
            "EPCAM", "NCAM1", "CLDN3", "PROM1", "SOX17", ## HSC markers; EPCAM also HB marker
            "KRT8", "KRT18", "KRT19", ## HSC markers
            "CXCR4", ## Endoderm marker expressed in HSCs and HCC
            "AFP", "KRT7", "ICAM1", ## HB markers; AFP also fetal hepatocyte; KRT7 (CK7) also cholangiocyte marker
            "HNF4A", ## HSC-HB and HB-hepatocyte differentiation?
            "ONECUT2", ## HB migration
            "PROX1", "TBX3", ## HB proliferation and migration
            "SOX9", ## Hepatic progenitor and cholangiocyte
            "HNF1B", "SALL4", ## Cholangiocyte fate regulator
            "CEBPA", "ALB", "TTR", "RBPJ", "NR5A2", ## Hepatocytic cell fate
            "MAT1A", "NR1I2", ## Adult liver
            "MAT2A", ## Fetal liver and replaces MAT1A in HCC
            "TAF10", "TBP", ## Embryonic liver not adult
            "TAF4", ## Postnatal hepatocytes
            "TGFB1", "TGFB2", "TGFBR2") ## High TGFb = cholangiocyte differentiation

pdf("signac_plots/huh6_plots/coverage_plots_custom_markers.pdf")
for (gene in custom) {
  print(CoveragePlot(object = combined,
                     region = gene,
                     extend.upstream = 10000,
                     extend.downstream = 10000))
}
dev.off()

## Plot genomic regions for cells grouped by cluster
Idents(combined) <- combined$PAM_Name
combined <- SortIdents(combined)

## Custom markers
pdf("signac_plots/huh6_plots/coverage_plots_custom_markers_pam_cluster.pdf")
for (gene in custom) {
  print(CoveragePlot(object = combined,
                     region = gene,
                     extend.upstream = 10000,
                     extend.downstream = 10000))
}
dev.off()

## [ Link Peaks to Genes ] ----

## For each gene, we can find the set of peaks that may regulate the gene
combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")
DefaultAssay(combined) <- "ATAC"

## Compute the GC content for each peak to correct for bias
combined <- RegionStats(combined, genome = BSgenome.Hsapiens.UCSC.hg38)

## Link peaks to HNF1A and LEF1 genes
combined <- LinkPeaks(object = combined,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("HNF1A", "LEF1")) ## No significant links found

## Link peaks to gene expression of custom markers
custom <- c("CD34", "PTPRC", "MCAM", ## Generally not expressed in normal liver, high in HCC; PTPRC = CD45 expressed in HPCs
            "FOXA1", "FOXA2", "GATA4", ## Early hepatic specification
            "EPCAM", "NCAM1", "CLDN3", "PROM1", "SOX17", ## HSC markers; EPCAM also HB marker
            "KRT8", "KRT18", "KRT19", ## HSC markers
            "CXCR4", ## Endoderm marker expressed in HSCs and HCC
            "AFP", "KRT7", "ICAM1", ## HB markers; AFP also fetal hepatocyte; KRT7 (CK7) also cholangiocyte marker
            "HNF4A", ## HSC-HB and HB-hepatocyte differentiation?
            "ONECUT2", ## HB migration
            "PROX1", "TBX3", ## HB proliferation and migration
            "SOX9", ## Hepatic progenitor and cholangiocyte
            "HNF1B", "SALL4", ## Cholangiocyte fate regulator
            "CEBPA", "ALB", "TTR", "RBPJ", "NR5A2", ## Hepatocytic cell fate
            "MAT1A", "NR1I2", ## Adult liver
            "MAT2A", ## Fetal liver and replaces MAT1A in HCC
            "TAF10", "TBP", ## Embryonic liver not adult
            "TAF4", ## Postnatal hepatocytes
            "TGFB1", "TGFB2", "TGFBR2") ## High TGFb = cholangiocyte differentiation

combined <- LinkPeaks(object = combined,
                      peak.assay = "ATAC",
                      expression.assay = "SCT",
                      genes.use = custom)

custom.filt <- custom[custom %in% rownames(combined@assays$SCT)]

Idents(combined) <- combined$Condition
idents.plot <- unique(combined$Condition)
pdf("signac_plots/huh6_plots/coverage_plots_custom_markers_link_peaks.pdf")
for (gene in custom.filt) {
  print(CoveragePlot(object = combined,
                     region = gene,
                     features = gene,
                     expression.assay = "SCT",
                     idents = idents.plot,
                     extend.upstream = 10000,
                     extend.downstream = 10000) &
          scale_fill_manual(values = c("Untreated" = "#BBCCEE",
                                       "Recovered" = "#CCDDAA")))
}
dev.off()

## Group cells by PAM cluster
Idents(combined) <- combined$PAM_Name
idents.plot <- unique(combined$PAM_Name)
pdf("signac_plots/huh6_plots/coverage_plots_custom_markers_pam_cluster_link_peaks.pdf")
for (gene in custom.filt) {
  print(CoveragePlot(object = combined,
                     region = gene,
                     features = gene,
                     expression.assay = "SCT",
                     idents = idents.plot,
                     extend.upstream = 10000,
                     extend.downstream = 10000) &
          scale_fill_manual(values = c(pam.cols)))
}
dev.off()

## [ Motif Analysis ] ----

combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")

## Load JASPAR2024 database
jaspar <- JASPAR2024()
jaspar.sq <- dbConnect(SQLite(), db(jaspar))

## Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(jaspar.sq,
                    list(collection = "CORE",
                         species = "Homo sapiens"))

combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)

## Load in differential peaks between clusters
pam.marker.list <- readRDS("signac_data/huh6_data/pam_markers.rds")
## Get top differentially accessible peaks (already filtered this list for adjusted p-value)
top.markers <- list()
for (name in names(pam.marker.list)) {
  filtered <- pam.marker.list[[name]][pam.marker.list[[name]]$p_val_adj < 0.05 &
                                        pam.marker.list[[name]]$avg_log2FC > 1, ]$query_region ## Only extract region
  top.markers[[name]] <- filtered
}

## Choose a set of background peaks
Idents(combined) <- combined$PAM_Name
open.peaks <- AccessiblePeaks(combined)

## Match the overall GC content in the peak set
meta.feature <- GetAssayData(combined, assay = "ATAC", layer = "meta.features")
peaks.matched <- list()
for (name in names(top.markers)) {
  matched <- MatchRegionStats(meta.feature = meta.feature[open.peaks, ],
                              query.feature = meta.feature[top.markers[[name]], ],
                              n = 50000)
  peaks.matched[[name]] <- matched
}

enriched.motifs <- list()
for (name in names(top.markers)) {
  enriched <- FindMotifs(object = combined,
                         features = top.markers[[name]],
                         background = peaks.matched[[name]])
  enriched.motifs[[name]] <- enriched
}
saveRDS(enriched.motifs, "signac_data/huh6_data/pam_markers_enriched_motifs.rds")

motifs.sheets <- lapply(enriched.motifs, as.data.frame)
write_xlsx(motifs.sheets, "signac_data/huh6_data/pam_markers_enriched_motifs.xlsx")

## [ Motif Footprinting ] ----

## Footprinting analysis on data before filtering for cells with snRNA-seq annotation
## Compare between untreated and recovered
combined <- readRDS("signac_data/huh6_data/combined_seurat_norm.rds")
up.tfs <- readRDS('/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/trajectory_genes/tfs_common_up.rds')

## Load JASPAR2024 database
jaspar <- JASPAR2024()
jaspar.sq <- dbConnect(SQLite(), db(jaspar))

## Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(jaspar.sq,
                    list(collection = "CORE",
                         species = "Homo sapiens"))

combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_motifs.rds")

## Footprinting for TF motifs enriched in pseudotime up genes
tf.names <- unlist(lapply(pfm, function(x) x@name))
filt.tfs <- up.tfs[up.tfs %in% tf.names]

Idents(combined) <- combined$Sample
combined <- Footprint(object = combined,
                      motif.name = filt.tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)

pdf("signac_plots/huh6_plots/motif_footprinting_sample.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf) &
          scale_colour_manual(values = huh6.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_sample_split.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "Sample") &
          scale_colour_manual(values = huh6.cols))
}
dev.off()

## Footprinting analysis on data after filtering for cells with snRNA-seq annotation
## Compare between PAM clusters annotated by snRNA-seq data
combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")

DefaultAssay(combined) <- "ATAC"
combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn_motifs.rds")

combined@assays$Activity <- NULL
combined@assays$RNA <- NULL
combined@assays$SCT <- NULL
Idents(combined) <- combined$PAM_Name
combined <- Footprint(object = combined,
                      motif.name = filt.tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn_footprint.rds")

pdf("signac_plots/huh6_plots/motif_footprinting_cluster.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf) &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_cluster_split.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "PAM_Name") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_cluster_sample.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "Sample") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

## Footprinting analysis on data before filtering for cells with snRNA-seq annotation
## Compare between clusters annotated using scRNA-seq data as a reference
combined <- readRDS("signac_data/huh6_data/combined_seurat_pam_annotated.rds")
up.tfs <- readRDS('/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/trajectory_genes/tfs_common_up.rds')

# Load JASPAR2024 database
jaspar <- JASPAR2024()
jaspar.sq <- dbConnect(SQLite(), db(jaspar))

## Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(jaspar.sq,
                    list(collection = "CORE",
                         species = "Homo sapiens"))

DefaultAssay(combined) <- "ATAC"
combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_pam_annotated_motifs.rds")

## Footprinting for TF motifs enriched in pseudotime up genes
tf.names <- unlist(lapply(pfm, function(x) x@name))
filt.tfs <- up.tfs[up.tfs %in% tf.names]

combined@assays$Activity <- NULL
Idents(combined) <- combined$predicted.id
combined <- Footprint(object = combined,
                      motif.name = filt.tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_pam_annotated_footprint.rds")

pdf("signac_plots/huh6_plots/motif_footprinting_predicted_cluster.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf) &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_predicted_cluster_split.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "predicted.id") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_predicted_cluster_sample.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "Sample") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

## /////////////////////////////////////////////////////////////////////////////
## Unvalidated collection //////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Footprinting for CEBPZ and YBX1 - not in JASPAR CORE collection
## Data before filtering for cells with snRNA-seq annotation
## Compare between untreated and recovered
combined <- readRDS("signac_data/huh6_data/combined_seurat_norm.rds")
up.tfs <- readRDS('/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/trajectory_genes/tfs_common_up.rds')

## Load JASPAR2024 database
jaspar <- JASPAR2024()
jaspar.sq <- dbConnect(SQLite(), db(jaspar))

## Get a list of motif position frequency matrices from the JASPAR database (UNVALIDATED collection)
pfm <- getMatrixSet(jaspar.sq,
                    list(collection = "UNVALIDATED",
                         species = "Homo sapiens"))

combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_motifs_unvalidated.rds")

## Footprinting for TF motifs enriched in pseudotime up genes
tf.names <- unlist(lapply(pfm, function(x) x@name))
filt.tfs <- up.tfs[up.tfs %in% tf.names]

Idents(combined) <- combined$Sample
combined <- Footprint(object = combined,
                      motif.name = filt.tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_footprint_unvalidated.rds")

pdf("signac_plots/huh6_plots/motif_footprinting_predicted_sample_unvalidated.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf) &
          scale_colour_manual(values = huh6.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_predicted_sample_split_unvalidated.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "Sample") &
          scale_colour_manual(values = huh6.cols))
}
dev.off()

## Footprinting analysis on data after filtering for cells with snRNA-seq annotation
## Compare between PAM clusters annotated by snRNA-seq data
combined <- readRDS("signac_data/huh6_data/combined_seurat_filt_wnn.rds")

DefaultAssay(combined) <- "ATAC"
combined <- AddMotifs(object = combined,
                      genome = BSgenome.Hsapiens.UCSC.hg38,
                      pfm = pfm)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn_motifs_unvalidated.rds")

combined@assays$Activity <- NULL
combined@assays$RNA <- NULL
combined@assays$SCT <- NULL
Idents(combined) <- combined$PAM_Name
combined <- Footprint(object = combined,
                      motif.name = filt.tfs,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(combined, "signac_data/huh6_data/combined_seurat_filt_wnn_footprint_unvalidated.rds")

pdf("signac_plots/huh6_plots/motif_footprinting_cluster_unvalidated.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf) &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_cluster_split_unvalidated.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "PAM_Name") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

pdf("signac_plots/huh6_plots/motif_footprinting_cluster_sample_unvalidated.pdf")
for (tf in filt.tfs) {
  print(PlotFootprint(combined, features = tf, split.by = "Sample") &
          scale_colour_manual(values = pam.cols))
}
dev.off()

## [ Annotate Regions ] ----

pam.marker.list <- readRDS("signac_data/huh6_data/pam_markers.rds")

## Filter top differentially accessible peaks by log2 fold change
top.markers <- list()
for (name in names(pam.marker.list)) {
  filtered <- pam.marker.list[[name]][pam.marker.list[[name]]$p_val_adj < 0.05 &
                                        pam.marker.list[[name]]$avg_log2FC > 1, ]
  top.markers[[name]] <- filtered
}
saveRDS(top.markers, "signac_data/huh6_data/pam_markers_filt.rds")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annotated.list <- list()
for (name in names(top.markers)) {
  df <- top.markers[[name]]
  df <- df %>% 
    separate(query_region, into = c("chrom", "start", "end"), sep = "-") %>% 
    mutate(start = as.numeric(start), end = as.numeric(end))
  
  gr <- GRanges(seqnames = df$chrom,
                ranges = IRanges(start = df$start, end = df$end),
                gene_name = df$gene_name,
                cluster = df$cluster)
  
  peak.annotations <- annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, 
                                   annoDb = "org.Hs.eg.db")
  
  annotated.peaks <- as.data.frame(peak.annotations)
  annotated.peaks$chrom <- annotated.peaks$seqnames
  annotated.peaks$seqnames <- NULL
  annotated.df <- merge(df, annotated.peaks)
  
  annotated.list[[name]] <- annotated.df
}
saveRDS(annotated.list, "signac_data/huh6_data/pam_markers_annotated.rds")

annotated.sheets <- lapply(annotated.list, as.data.frame)
write_xlsx(annotated.sheets, "signac_data/huh6_data/pam_markers_annotated.xlsx")

head(unique(annotated.list$`Progenitor-hepatocytic`$annotation), 10)
annotated.list <- lapply(annotated.list, function(df) {
  df$genome_annotation <- gsub("\\s*\\(.*?\\)", "", df$annotation)
  return(df)
})
unique(annotated.list$`Progenitor-hepatocytic`$genome_annotation)

pal <- paletteer_d("ggthemes::Classic_10_Medium")
genome.cols <- c("Distal Intergenic" = pal[1],
                 "Intron" = pal[2],
                 "Exon" = pal[3],
                 "Promoter" = pal[4],
                 "3' UTR" = pal[5],
                 "5' UTR" = pal[6],
                 "Downstream" = pal[7])

clean.string <- function(string) {
  lower.string <- tolower(string)
  cleaned.string <- gsub(" ", "", lower.string)
  return(cleaned.string)
}

for (name in names(annotated.list)) {
  df <- annotated.list[[name]]
  freq <- as.data.frame(table(df$genome_annotation))
  colnames(freq) <- c("genome_annotation", "frequency")
  
  p <- ggplot(freq, aes(x = genome_annotation, y = frequency, fill = genome_annotation)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Frequency of Annotations in", name),
         x = "Annotation Type",
         y = "Frequency") +
    scale_fill_manual(values = genome.cols) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),
          legend.title = element_blank())
  p + coord_flip()
  file.name <- clean.string(name)
  ggsave(paste0("signac_plots/huh6_plots/pam_markers_", file.name, "_genome_annotation.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("signac_plots/huh6_plots/pam_markers_", file.name, "_genome_annotation.png"), width = 8.3, height = 5.8)
}

combined.df <- bind_rows(annotated.list, .id = "source")
freq.df <- combined.df %>%
  group_by(cluster, genome_annotation) %>%
  summarise(freq = n(), .groups = "drop")

ggplot(freq.df, aes(x = genome_annotation, y = freq, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = pam.cols) +
  labs(title = "Annotation Frequency by PAM Cluster",
       x = "Cluster",
       y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        legend.title = element_blank())
ggsave("signac_plots/huh6_plots/pam_markers_combined_genome_annotation.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/pam_markers_combined_genome_annotation.png", width = 8.3, height = 5.8)

prop.df <- combined.df %>%
  group_by(cluster, genome_annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count))

ggplot(prop.df, aes(x = genome_annotation, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = pam.cols) +
  labs(title = "Annotation Proportion by PAM Cluster",
       x = "Cluster",
       y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        legend.title = element_blank())
ggsave("signac_plots/huh6_plots/pam_markers_combined_genome_annotation_proportion.pdf", width = 8.3, height = 5.8)
ggsave("signac_plots/huh6_plots/pam_markers_combined_genome_annotation_proportion.png", width = 8.3, height = 5.8)
