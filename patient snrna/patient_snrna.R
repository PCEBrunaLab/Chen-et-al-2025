## HB patients: snRNA-seq alone

## /////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

renv::load()

library(Seurat)
library(scran)
library(scater)
library(scuttle)
library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(patchwork)
library(reshape2)
library(pheatmap)

## [ Plot themes ] ----

## Custom ggplot theme
umap.theme <- function() {
  require(ggplot2)
  return(theme_bw() +
           theme(aspect.ratio = 1, 
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line = element_line(colour = "#16161D", linewidth = 0.8),
                 axis.ticks = element_line(colour = "#16161D", linewidth = 0.8),
                 legend.title = element_blank()))
}

condition.cols <- c("12-00122_F1" = "#F7C9E3",
                    "12-00122_F2" = "#E58ABD",
                    "12-01026_F1" = "#FFE0B8",
                    "12-01026_F2" = "#FFAA55",
                    "16-00875_F1" = "#CDEBFB",
                    "16-00875_F2" = "#76C4E4")

## [ QC filtering ] ----

## Load in data
patient.sce <- readRDS("data/HB_patients_SCE-norm_relaxed.RDS")
length(unique(colnames(patient.sce))) ## 75897 cells

## Compute logcounts
patient.sce <- logNormCounts(patient.sce)

## ENSEMBL database
# library(EnsDb.Hsapiens.v86)
# library(biomaRt)
# ens.86.genes <- genes(EnsDb.Hsapiens.v86)
# linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype == "lincRNA"]
# ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
# ens.bm <- getBM(attributes = c("ensembl_gene_id",
#                                "chromosome_name",
#                                "start_position",
#                                "end_position",
#                                "external_gene_name",
#                                "hgnc_symbol",
#                                "description"),
#                 mart = ensembl,
#                 filters = "ensembl_gene_id",
#                 values = rownames(patient.sce))
# # Save ensembl.csv for reuse
# write.csv(ens.bm, "data/patient_ensembl_biomart.csv")
# gc()
ens.bm <- read.csv("data/patient_ensembl_biomart.csv", row.names = NULL)
ens.bm$X <- NULL

## Identify duplicated gene names
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)] 
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])
## Duplicated gene names: LINC00595

## Created filtered ens.bm and keep only useful information
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)
## This filters out around 1000 genes
nrow(ens.filt.bm) ## 35456
nrow(patient.sce) ## 36601
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]
## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]
## Around 10000 genes filtered out
nrow(ens.filt.bm) ## 25577
nrow(patient.sce) ## 36601
## Keep track of everything that didn't get annotated
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 9879
length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 9733 are novel transcripts or novel transcripts antisense to a gene
grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest are mostly either novel, mitochondrially encoded, or just empty description
## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## "ELFN2"  "GOLGA8M" "LINC00595" "LINC01115" "LINC03021" "LINC03023" "LINC03025" "RAET1E-AS1"  "SPATA13"  

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 25568

## Identify mitochondrial, ribosomal and linc genes
mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Add QC metrics to SCE object
patient.sce <- addPerCellQCMetrics(patient.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))
patient.sce <- patient.sce[ens.filt.bm$ensembl_gene_id,]
rownames(patient.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

## Assign rowData to SCE object with filtered ens.bm
rowData(patient.sce) <- ens.filt.bm
saveRDS(patient.sce, "data/patient_sce.rds")

## Extract meta data from SCE object
patient.meta <- as.data.frame(colData(patient.sce))
patient.meta$batch <- NULL
## Generate Seurat object from SCE object
patient.seurat <- CreateSeuratObject(counts = assay(patient.sce, "counts"),
                                     assay = "RNA",
                                     meta.data = patient.meta)
saveRDS(patient.seurat, "data/patient_seurat.rds")

patient.sce$Condition <- paste(patient.sce$Patient_ID, patient.sce$Fraction, sep = "_")
patient.sce$Condition_Rep <- paste(patient.sce$Condition, patient.sce$Replicate, sep = "_")
patient.seurat$Condition <- paste(patient.seurat$Patient_ID, patient.seurat$Fraction, sep = "_")
patient.seurat$Condition_Rep <- paste(patient.seurat$Condition, patient.seurat$Replicate, sep = "_")
## Plot QC metrics
dir.create("plots/qc", recursive = TRUE)
det.sce.gg <- plotColData(patient.sce, y = "detected", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
umi.sce.gg <- plotColData(patient.sce, y = "total", x = "Condition_Rep", colour_by = "Condition_Rep") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")
mt.sce.gg <- plotColData(patient.sce, y = "subsets_mt_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
rb.sce.gg <- plotColData(patient.sce, y = "subsets_ribo_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.8,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/qc/patient_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/patient_qc_plots.png", width = 8.3, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
counts.assay <- counts(patient.sce)
counts.assay@x <- counts.assay@x/rep.int(colSums(counts.assay), diff(counts.assay@p))
top.expr <- order(Matrix::rowSums(counts.assay), decreasing = TRUE)[20:1]
top.expr <- as.matrix(t(counts.assay[top.expr, ]))
top.expr <- as.data.frame(top.expr)
ord.idx <- colnames(top.expr)
top.expr <- melt(top.expr, variable.name = "Gene", value.name = "Prop")
top.expr$Gene <- factor(top.expr$Gene, levels = ord.idx)

ggplot(top.expr, mapping = aes(x = Gene, y = Prop * 100)) +
  geom_boxplot() +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dotted")) +
  coord_flip() +
  labs(y = "% total expression", x = "Gene", title = "Percentage expression of a single gene per total cell expression")
ggsave("plots/qc/patient_qc_top_genes.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/patient_qc_top_genes.png", width = 8.3, height = 5.8)

rm(counts.assay)
gc()

## Visualise MALAT1 expression and lincRNA per sample
malat1.sce.gg <- plotExpression(patient.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Condition_Rep", 
                                colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Log2 Expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
lincrna.sce.gg <- plotColData(patient.sce, y = "subsets_linc_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.8,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/qc/patient_qc_malat1_lincrna.pdf", width = 5.8, height = 8.3)
ggsave("plots/qc/patient_qc_malat1_lincrna.png", width = 5.8, height = 8.3)

## Cell cycle scoring
patient.seurat <- NormalizeData(patient.seurat)
patient.seurat <- CellCycleScoring(patient.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
patient.sce$Seurat.Phase <- patient.seurat$Phase
patient.sce$Seurat.S <- patient.seurat$S.Score
patient.sce$Seurat.G2M <- patient.seurat$G2M.Score
saveRDS(patient.sce, "data/patient_sce_cell_cycle.rds")
saveRDS(patient.seurat, "data/patient_seurat_cell_cycle.rds")

# patient.sce <- readRDS("data/patient_sce_cell_cycle.rds")
# patient.seurat <- readRDS("data/patient_seurat_cell_cycle.rds")
## Check gene and cell meta data
colnames(colData(patient.sce))
colnames(rowData(patient.sce))
table(patient.sce$Condition)
## 12-00122_F1 12-00122_F2 12-01026_F1 12-01026_F2 16-00875_F1 16-00875_F2 
## 11785          73       23936       13347       19621        7135 
table(patient.sce$Condition_Rep)
## 12-00122_F1_i  12-00122_F1_ii   12-00122_F2_i   12-01026_F1_i  12-01026_F1_ii 12-01026_F1_iii   12-01026_F2_i  12-01026_F2_ii   16-00875_F1_i 
## 1401           10384              73            5241           11744            6951            6773            6574           10575 
## 16-00875_F1_ii   16-00875_F2_i  16-00875_F2_ii 
## 9046            5773            1362 

det.sce.gg <- plotColData(patient.sce, y = "detected", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dashed")
  # geom_hline(yintercept = 200, linetype = "dotted")
mt.sce.gg <- plotColData(patient.sce, y = "subsets_mt_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 15, linetype = "dotted")
rb.sce.gg <- plotColData(patient.sce, y = "subsets_ribo_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed")
det.sce.gg + mt.sce.gg + rb.sce.gg &
  theme(aspect.ratio = 0.7,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/qc/patient_qc_plots_filtering.pdf", width = 11.7, height = 5.8)
ggsave("plots/qc/patient_qc_plots_filtering.png", width = 11.7, height = 5.8)

## No need for expression filtering
detected.filt <- colnames(patient.sce)[patient.sce$detected > 1000] ## No cells removed
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(patient.sce)[Matrix::rowSums(counts(patient.sce)) > 5] ## 2375 genes removed
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(patient.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)
patient.filt.sce <- patient.sce[genes.filt, detected.filt]
dim(patient.sce) - dim(patient.filt.sce)
## We lose 2375 genes and 0 cells

## Gene expression of each cell has to have less than 15% mitochondrial
mito.filt <- colnames(patient.filt.sce)[patient.filt.sce$subsets_mt_percent < 15]
## No ribosomal filtering - too harsh for snRNA-seq
#ribo.filt <- colnames(patient.filt.sce)[patient.filt.sce$subsets_ribo_percent > 5]
#qc.filt <- intersect(mito.filt, ribo.filt)

patient.filt.sce <- patient.filt.sce[, mito.filt]
dim(patient.sce) - dim(patient.filt.sce)
## We lose 195 cells in total

## Percentage of the library that we filter out
round(100 - (100 * (table(patient.filt.sce$Condition_Rep) / table(patient.sce$Condition_Rep))), 1)
## 12-00122_F1_i  12-00122_F1_ii   12-00122_F2_i   12-01026_F1_i  12-01026_F1_ii 12-01026_F1_iii   12-01026_F2_i  12-01026_F2_ii   16-00875_F1_i 
## 0.0             0.0             0.0             0.2             0.0             0.0             1.3             0.3             0.5 
## 16-00875_F1_ii   16-00875_F2_i  16-00875_F2_ii 
## 0.2             0.0             0.0 

table(patient.filt.sce$Condition_Rep)
## 12-00122_F1_i  12-00122_F1_ii   12-00122_F2_i   12-01026_F1_i  12-01026_F1_ii 12-01026_F1_iii   12-01026_F2_i  12-01026_F2_ii   16-00875_F1_i 
## 1401           10384              73            5228           11742            6951            6686            6551           10524 
## 16-00875_F1_ii   16-00875_F2_i  16-00875_F2_ii 
## 9027            5773            1362 

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(patient.seurat)[Matrix::rowSums(GetAssayData(patient.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(patient.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

patient.filt.seurat <- subset(patient.seurat, 
                              cells = WhichCells(patient.seurat, 
                                                 expression = subsets_mt_percent < 15 &
                                                   detected > 1000),
                              features = genes.filt)

## The dimensions are equal
dim(patient.filt.seurat) == dim(patient.filt.sce)

saveRDS(patient.filt.sce, "data/patient_sce_filtered.rds")
saveRDS(patient.filt.seurat, "data/patient_seurat_filtered.rds")

## [ Filter HVGs ] ----

patient.filt.seurat <- readRDS("data/patient_seurat_filtered.rds")

## Yura's script to filter HVGs
patient.filt.seurat <- NormalizeData(patient.filt.seurat)
patient.filt.seurat <- FindVariableFeatures(patient.filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(patient.filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(patient.filt.seurat), 20), repel = TRUE)
ggsave("plots/filter_hvgs/patient_filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/patient_filter_hvgs.png", width = 8.3, height = 5.8)

## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(patient.filt.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(patient.filt.seurat), value = TRUE),
              grep("[\\.]", rownames(patient.filt.seurat), value = TRUE),
              grep("^LINC", rownames(patient.filt.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
VariableFeatures(patient.filt.seurat) <- VariableFeatures(patient.filt.seurat)[which(!VariableFeatures(patient.filt.seurat) %in% rm.genes)]
## Scale data and score cell cycle
patient.filt.seurat <- ScaleData(patient.filt.seurat, features = VariableFeatures(patient.filt.seurat))

tmp.seurat <- CellCycleScoring(object = patient.filt.seurat, 
                               g2m.features = cc.genes$g2m.genes, 
                               s.features = cc.genes$s.genes)
tmp.seurat$Cycle.Score <- tmp.seurat$S.Score - tmp.seurat$G2M.Score

seurat.cycle.melt <- ggsankey::make_long(
  data.frame("Sample" = patient.filt.seurat$Sample, 
             "SCTransform" = tmp.seurat$Phase,
             "Normal" = patient.filt.seurat$Phase),
  SCTransform, Normal)

library(ggsankey)
ggplot(seurat.cycle.melt, 
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node),
           label = node)) +
  geom_sankey(flow.alpha = 0.7) +
  geom_sankey_label() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = element_blank(), title = "Seurat Cell Cycle Prediction") +
  scale_fill_manual(values = hues::iwanthue(3))
## Not a huge difference
rm(tmp.seurat)
gc()

DefaultAssay(patient.filt.seurat) <- "RNA"
patient.filt.seurat$Seurat.Cycle.Score <- patient.filt.seurat$S.Score - patient.filt.seurat$G2M.Score
saveRDS(patient.filt.seurat, "data/patient_seurat_hvgs.rds")

# patient.seurat <- readRDS("data/patient_seurat_hvgs.rds")
patient.seurat <- RunPCA(patient.seurat, verbose = FALSE, npcs = 100, 
                         features = VariableFeatures(patient.seurat))
ElbowPlot(object = patient.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- harmony::RunHarmony(patient.seurat, group.by.vars = "Condition_Rep",
                                      theta = 0.5, lambda = 1, sigma = 0.05,
                                      assay.use = "RNA", reduction = "pca",
                                      dims.use = 1:30, reduction.save = "Harmony",
                                      max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition_Rep", order = TRUE) +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("plots/filter_hvgs/umap_harmony_t05_l1_s005.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/umap_harmony_t05_l1_s005.png", width = 8.3, height = 5.8)

# patient.seurat <- RunUMAP(patient.seurat, reduction = "pca", dims = 1:30)
# patient.seurat <- FindNeighbors(patient.seurat, reduction = "pca", dims = 1:30)

cluster.search <- function(seurat, from = 0.2, to = 0.8, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("RNA_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}
harmony.seurat <- cluster.search(harmony.seurat)
saveRDS(harmony.seurat, "data/patient_seurat_hvgs_harmony.rds")

patient.seurat <- readRDS("data/patient_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
patient.seurat <- readRDS("data/patient_seurat_hvgs.rds")
obj.list <- SplitObject(patient.seurat, split.by = "Condition")

## Rerun variable feature selection
obj.list <- lapply(obj.list, FindVariableFeatures)
## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(patient.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(patient.seurat), value = TRUE),
              grep("[\\.]", rownames(patient.seurat), value = TRUE),
              grep("^LINC", rownames(patient.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(patient.seurat, obj.list, obj)

## Scale data
output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
## Run PCA - can't use 100 PCs for object with 73 cells so use 50
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 50,
         features = VariableFeatures(obj))
})

## Run UMAP and compute nearest neighbours 
output <- lapply(output, RunUMAP, reduction = "pca", dims = 1:30)
output <- lapply(output, FindNeighbors, reduction = "pca", dims = 1:30)
## Clustering 
cluster.search <- function(seurat, from = 0.2, to = 0.8, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("RNA_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}
output <- lapply(output, cluster.search)
saveRDS(output, "data/patient_seurat_split_list.rds")

patient.list <- readRDS("data/patient_seurat_split_list.rds")
patient.seurat <- readRDS("data/patient_seurat_hvgs_harmony.rds")
## Look through UMAPs and choose best cluster resolution
# DimPlot(patient.seurat, group.by = "seurat_clusters.0.4", order = TRUE) +
#   umap.theme()
Idents(patient.list$`12-00122_F1`) <- patient.list$`12-00122_F1`$seurat_clusters.0.2
Idents(patient.list$`12-00122_F2`) <- patient.list$`12-00122_F2`$seurat_clusters.0.2
Idents(patient.list$`12-01026_F1`) <- patient.list$`12-01026_F1`$seurat_clusters.0.2
Idents(patient.list$`12-01026_F2`) <- patient.list$`12-01026_F2`$seurat_clusters.0.4
Idents(patient.list$`16-00875_F1`) <- patient.list$`16-00875_F1`$seurat_clusters.0.2
Idents(patient.list$`16-00875_F2`) <- patient.list$`16-00875_F2`$seurat_clusters.0.2
Idents(patient.seurat) <- patient.seurat$seurat_clusters.0.4

patient.seurat$Patient_Num <- recode(as.character(patient.seurat$Patient_ID),
                                     "12-00122" = "P1",
                                     "12-01026" = "P2",
                                     "16-00875" = "P3")
## Get cluster markers - MAST takes longer to run
p1f1.markers <- FindAllMarkers(patient.list$`12-00122_F1`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
p1f2.markers <- FindAllMarkers(patient.list$`12-00122_F2`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE) ## No markers found
p2f1.markers <- FindAllMarkers(patient.list$`12-01026_F1`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
p2f2.markers <- FindAllMarkers(patient.list$`12-01026_F2`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
p3f1.markers <- FindAllMarkers(patient.list$`16-00875_F1`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
p3f2.markers <- FindAllMarkers(patient.list$`16-00875_F2`, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
all.markers <- FindAllMarkers(patient.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

dir.create("data/jaccard/", recursive = TRUE)
saveRDS(p1f1.markers, "data/jaccard/12-00122_F1_markers.rds")
saveRDS(p1f2.markers, "data/jaccard/12-00122_F2_markers.rds")
saveRDS(p2f1.markers, "data/jaccard/12-01026_F1_markers.rds")
saveRDS(p2f2.markers, "data/jaccard/12-01026_F2_markers.rds")
saveRDS(p3f1.markers, "data/jaccard/16-00875_F1_markers.rds")
saveRDS(p3f2.markers, "data/jaccard/16-00875_F2_markers.rds")
saveRDS(all.markers, "data/jaccard/patient_all_markers.rds")

## Filter significant and higher expressed markers
p1f1.markers <- p1f1.markers[p1f1.markers$p_val_adj < 0.05, ]
p1f2.markers <- p1f2.markers[p1f2.markers$p_val_adj < 0.05, ]
p2f1.markers <- p2f1.markers[p2f1.markers$p_val_adj < 0.05, ]
p2f2.markers <- p2f2.markers[p2f2.markers$p_val_adj < 0.05, ]
p3f1.markers <- p3f1.markers[p3f1.markers$p_val_adj < 0.05, ]
p3f2.markers <- p3f2.markers[p3f2.markers$p_val_adj < 0.05, ]
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]

p1f1.markers <- p1f1.markers[p1f1.markers$avg_log2FC > 0.5, ]
p1f2.markers <- p1f2.markers[p1f2.markers$avg_log2FC > 0.5, ]
p2f1.markers <- p2f1.markers[p2f1.markers$avg_log2FC > 0.5, ]
p2f2.markers <- p2f2.markers[p2f2.markers$avg_log2FC > 0.5, ]
p3f1.markers <- p3f1.markers[p3f1.markers$avg_log2FC > 0.5, ]
p3f2.markers <- p3f2.markers[p3f2.markers$avg_log2FC > 0.5, ]
all.markers <- all.markers[all.markers$avg_log2FC > 0.5, ]

## Split into lists of genes per cluster
p1f1.markers <- split(p1f1.markers, p1f1.markers$cluster)
p1f2.markers <- split(p1f2.markers, p1f2.markers$cluster)
p2f1.markers <- split(p2f1.markers, p2f1.markers$cluster)
p2f2.markers <- split(p2f2.markers, p2f2.markers$cluster)
p3f1.markers <- split(p3f1.markers, p3f1.markers$cluster)
p3f2.markers <- split(p3f2.markers, p3f2.markers$cluster)
all.markers <- split(all.markers, all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(p1f1.markers, nrow) ## 0 = 1009, 1 = 601, 2 = 1559, 3 = 443, 4 = 1117, 5 = 1158, 6 = 635, 7 = 1009, 8 = 1416
lapply(p1f2.markers, nrow) ## No markers
lapply(p2f1.markers, nrow) ## 0 = 214, 1 = 538, 2 = 631, 3 = 502, 4 = 304, 5 = 523, 6 = 292, 7 = 1336, 8 = 632, 9 = 1244
lapply(p2f2.markers, nrow) ## 0 = 226, 1 = 118, 2 = 339, 3 = 394, 4 = 598, 5 = 438, 6 = 643, 7 = 572, 8 = 690
lapply(p3f1.markers, nrow) ## 0 = 652, 1 = 763, 2 = 1197, 3 = 617, 4 = 1120, 5 = 1178, 6 = 1621, 7 = 1603, 8 = 1349, 9 = 303, 10 = 1165, 11 = 471
lapply(p3f2.markers, nrow) ## 0 = 551, 1 = 332, 2 = 399, 3 = 249, 4 = 1516, 5 = 442, 6 = 1082, 7 = 162, 8 = 1098, 9 = 957, 10 = 778
lapply(all.markers, nrow) ## 0 = 1127, 1 = 931, 2 = 1494, 3 = 1287, 4 = 1390, 5 = 1607, 6 = 1337, 7 = 803, 8 = 976, 9 = 1547, 10 = 1493,
## 11 = 897, 12 = 1275, 13 = 1228, 14 = 808, 15 = 1507, 16 = 1800, 17 = 760, 18 = 1201, 19 = 1415, 20 = 1166, 21 = 1104

p1f1.markers <- lapply(p1f1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 400)
})
# p1f2.markers <- lapply(p1f2.markers, \(x) {
#   x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
#   head(x, 100)
# })
p2f1.markers <- lapply(p2f1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
p2f2.markers <- lapply(p2f2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
p3f1.markers <- lapply(p3f1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
p3f2.markers <- lapply(p3f2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
all.markers <- lapply(all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})

## Give your clusters unique names
names(p1f1.markers) <- paste0("P1_F1_", names(p1f1.markers))
# names(p1f2.markers) <- paste0("P1_F2_", names(p1f2.markers))
names(p2f1.markers) <- paste0("P2_F1_", names(p2f1.markers))
names(p2f2.markers) <- paste0("P2_F2_", names(p2f2.markers))
names(p3f1.markers) <- paste0("P3_F1_", names(p3f1.markers))
names(p3f2.markers) <- paste0("P3_F2_", names(p3f2.markers))
names(all.markers) <- paste0("All_", names(all.markers))

## Make a big list of all the cluster markers
patient.markers.list <- c(p1f1.markers,
                          p2f1.markers, p2f2.markers,
                          p3f1.markers, p3f2.markers,
                          all.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
patient.jaccard.mat <- matrix(data = NA, nrow = length(patient.markers.list),
                              ncol = length(patient.markers.list),
                              dimnames = list(names(patient.markers.list),
                                              names(patient.markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(patient.jaccard.mat)) {
  for (j in colnames(patient.jaccard.mat)) {
    patient.jaccard.mat[i,j] <- jaccard(patient.markers.list[[i]]$gene, patient.markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Sample" = gsub("_\\d+", "", colnames(patient.jaccard.mat)),
                      row.names = colnames(patient.jaccard.mat))
anno.col <- list("Sample" = c("All" = "black",
                              "P1_F1" = as.character(condition.cols[1]),
                              "P1_F2" = as.character(condition.cols[2]),
                              "P2_F1" = as.character(condition.cols[3]),
                              "P2_F2" = as.character(condition.cols[4]),
                              "P3_F1" = as.character(condition.cols[5]),
                              "P3_F2" = as.character(condition.cols[6])))
dev.off()
gt <- pheatmap(patient.jaccard.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_row = anno.df,
               annotation_col = anno.df,
               annotation_colors = anno.col)$gtable
ggsave("plots/jaccard/jaccard_heatmap.png", plot = gt, height = 16.5, width = 16.5)
ggsave("plots/jaccard/jaccard_heatmap.pdf", plot = gt, height = 16.5, width = 16.5)

## Run the parameter search for PAM
patient.jaccard.dist <- as.dist(1 - patient.jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
for (k in 2:30) {  # Assuming you want to check from 2 to 30 clusters
  pam.fit <- pam(patient.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:30
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 11
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:31))
## Draw the geom_vline at 10 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 10, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
patient.k11.pam <- pam(patient.jaccard.dist, k = 11, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- patient.k11.pam
## Assign the clusters back to your Seurat objects
patient.clusters.df <- anno.df$PAM.Cluster
names(patient.clusters.df) <- rownames(anno.df)
p1f1.set <- paste0("P1_F1_", as.character(patient.list$`12-00122_F1`$seurat_clusters.0.2))
p1f2.set <- paste0("P1_F2_", as.character(patient.list$`12-00122_F2`$seurat_clusters.0.2))
p2f1.set <- paste0("P2_F1_", as.character(patient.list$`12-01026_F1`$seurat_clusters.0.2))
p2f2.set <- paste0("P2_F2_", as.character(patient.list$`12-01026_F2`$seurat_clusters.0.4))
p3f1.set <- paste0("P3_F1_", as.character(patient.list$`16-00875_F1`$seurat_clusters.0.2))
p3f2.set <- paste0("P3_F2_", as.character(patient.list$`16-00875_F2`$seurat_clusters.0.2))
all.set <- paste0("All_", as.character(patient.seurat$seurat_clusters.0.4))

patient.list$`12-00122_F1`$PAM.Cluster <- factor(as.character(patient.clusters.df[p1f1.set]))
patient.list$`12-00122_F2`$PAM.Cluster <- factor(as.character(patient.clusters.df[p1f2.set]))
patient.list$`12-01026_F1`$PAM.Cluster <- factor(as.character(patient.clusters.df[p2f1.set]))
patient.list$`12-01026_F2`$PAM.Cluster <- factor(as.character(patient.clusters.df[p2f2.set]))
patient.list$`16-00875_F1`$PAM.Cluster <- factor(as.character(patient.clusters.df[p3f1.set]))
patient.list$`16-00875_F2`$PAM.Cluster <- factor(as.character(patient.clusters.df[p3f2.set]))
patient.seurat$PAM.Cluster <- factor(as.character(patient.clusters.df[all.set]))
saveRDS(patient.list, "data/patient_seurat_split_list_pam.rds")
saveRDS(patient.seurat, "data/patient_seurat_hvgs_pam.rds")

DimPlot(patient.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = condition.cols) +
  umap.theme() + labs(title = "Conditions")
ggsave("plots/umaps/umap_harmony_condition.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_condition.png", width = 8.3, height = 5.8)

patient.seurat$PAM.Cluster <- factor(patient.seurat$PAM.Cluster, levels = 1:11)
pam.cols <- hues::iwanthue(length(unique(patient.seurat$PAM.Cluster)))
names(pam.cols) <- 1:11
saveRDS(pam.cols, "data/pam_cluster_cols.rds")
DimPlot(patient.seurat, group.by = "PAM.Cluster", order = TRUE) +
  scale_colour_manual(values = pam.cols) +
  umap.theme() + labs(title = "PAM Cluster")
ggsave("plots/umaps/umap_harmony_pam.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_pam.png", width = 8.3, height = 5.8)

## Assign clusters from separated objects into combined
# merged <- scCustomize::Merge_Seurat_List(patient.list)
merged <- Reduce(function(x, y) merge(x, y), patient.list)
merged.meta <- merged@meta.data
length(rownames(merged.meta) %in% rownames(patient.seurat@meta.data)) ## 75702
merged.meta$CellID <- rownames(merged.meta)
patient.seurat$CellID <- rownames(patient.seurat@meta.data)
merged.ordered <- merged.meta[match(patient.seurat$CellID, merged.meta$CellID), ]
identical(rownames(merged.ordered), rownames(patient.seurat@meta.data)) ## TRUE
patient.seurat$Meta.Cluster <- merged.ordered$PAM.Cluster
patient.seurat$Meta.Cluster <- replace_na(patient.seurat$Meta.Cluster, "NA")

patient.seurat$Meta.Cluster <- factor(patient.seurat$Meta.Cluster, levels = c(1:11, "NA"))
meta.cols <- c(pam.cols, "NA" = "grey")
saveRDS(meta.cols, "data/meta_cluster_cols.rds")
DimPlot(patient.seurat, group.by = "Meta.Cluster", order = TRUE) +
  scale_colour_manual(values = meta.cols) +
  umap.theme() + labs(title = "Meta Cluster")
ggsave("plots/umaps/umap_harmony_meta.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_meta.png", width = 8.3, height = 5.8)

saveRDS(patient.seurat, "data/patient_seurat_hvgs_meta.rds")

## [ Cluster annotation ] ----

patient.seurat <- readRDS("data/patient_seurat_hvgs_meta.rds")
pam.cols <- readRDS("data/pam_cluster_cols.rds")
meta.cols <- readRDS("data/meta_cluster_cols.rds")

## HNF4A and LEF1 (Kluiver et al., 2023)
features = c("HNF4A", "LEF1")
Idents(patient.seurat) <- patient.seurat$PAM.Cluster
VlnPlot(patient.seurat, features = features,
        pt.size = 0, cols = pam.cols) &
  xlab("") &
  theme(legend.position = "none",
        text = element_text(size = 12),
        aspect.ratio = 0.5)
ggsave("plots/anno/hnf4a_lef1_violin_pam.pdf", width = 8.7, height = 5.8)
ggsave("plots/anno/hnf4a_lef1_violin_pam.png", width = 8.7, height = 5.8)

Idents(patient.seurat) <- patient.seurat$Meta.Cluster
VlnPlot(patient.seurat, features = features,
        pt.size = 0, cols = meta.cols) &
  xlab("") &
  theme(legend.position = "none",
        text = element_text(size = 12),
        aspect.ratio = 0.5)
ggsave("plots/anno/hnf4a_lef1_violin_meta.pdf", width = 8.7, height = 5.8)
ggsave("plots/anno/hnf4a_lef1_violin_meta.png", width = 8.7, height = 5.8)

## Use signature scores to annotate clusters
## HB signatures from literature
sig.df <- as.data.frame(readxl::read_xlsx(path = "data/hb_sigs.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")
sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})
names(sigs) <- unique(sig.df$Signature)
sigs
sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)
patient.seurat <- AddModuleScore(patient.seurat, features = sigs, assay = "RNA", seed = 12345, 
                                 name = sig.names)
colnames(patient.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(patient.seurat@meta.data))
patient.seurat[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = patient.seurat, vars = sig.names)))
DoHeatmap(patient.seurat, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "PAM.Cluster", label = FALSE, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/anno/hb_sigs_heatmap_pam.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/hb_sigs_heatmap_pam.png", width = 8.3, height = 8.3)
DoHeatmap(patient.seurat, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "Meta.Cluster", label = FALSE, group.colors = meta.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"), aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Meta"))
ggsave("plots/anno/hb_sigs_heatmap_meta.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/hb_sigs_heatmap_meta.png", width = 8.3, height = 8.3)

## MSigDB Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
H.list <- lapply(unique(H.df$gs_name), function(x){
  H.df[H.df$gs_name==x, "gene_symbol"]
})
names(H.list) <- unique(H.df$gs_name)
H.list
H.names <- paste(names(H.list), ".Sig", sep = "")
H.names <- gsub("_", "-", H.names)
patient.seurat <- AddModuleScore(patient.seurat, features = H.list, assay = "RNA", seed = 12345, 
                                 name = H.names)
colnames(patient.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(patient.seurat@meta.data))
patient.seurat[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = patient.seurat, vars = H.names)))
DoHeatmap(patient.seurat, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/anno/msigdb_H_heatmap_pam.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/msigdb_H_heatmap_pam.png", width = 8.3, height = 8.3)
DoHeatmap(patient.seurat, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "Meta.Cluster", size = 4, group.colors = meta.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Meta"))
ggsave("plots/anno/msigdb_H_heatmap_meta.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/msigdb_H_heatmap_meta.png", width = 8.3, height = 8.3)

## MSigDB C2 gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
liver.match <- c("HEPATOCYTE", "HEPATOBLAST", "LIVER", "EMT", "EPITHELIAL_MESENCHYMAL",
                 "TGFB", "TGF_BETA", "WNT", "NOTCH", "HEDGEHOG", "SHH")
liver.gs <- C2.df[grepl(paste(liver.match, collapse = "|"), C2.df$gs_name), ]
liver.list <- lapply(unique(liver.gs$gs_name), function(x){
  liver.gs[liver.gs$gs_name==x, "gene_symbol"]
})
names(liver.list) <- unique(liver.gs$gs_name)
liver.list
liver.names <- paste(names(liver.list), ".Sig", sep = "")
liver.names <- gsub("_", "-", liver.names)
patient.seurat <- AddModuleScore(patient.seurat, features = liver.list, assay = "RNA", seed = 12345, 
                                 name = liver.names)
colnames(patient.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(patient.seurat@meta.data))
patient.seurat[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = patient.seurat, vars = liver.names)))
abbr_sig <- function(x) {
  core <- sub("\\.Sig$", "", x)
  parts <- strsplit(core, "-", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(parts)
  paste(
    c(parts[1], substr(parts[-1], 1, 1)),
    collapse = "-"
  )
}
liver.names.abbr <- sapply(liver.names, abbr_sig, USE.NAMES = TRUE)
DoHeatmap(patient.seurat, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  scale_y_discrete(labels = liver.names.abbr) +
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 2)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/anno/msigdb_C2_liver_heatmap_pam.pdf", width = 8.3, height = 8.3)
# ggsave("plots/anno/msigdb_C2_liver_heatmap_pam.png", width = 8.3, height = 8.3)
DoHeatmap(patient.seurat, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "seurat_clusters.0.4", size = 4, group.colors = meta.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  scale_y_discrete(labels = liver.names.abbr) +
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 2)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Meta"))
ggsave("plots/anno/msigdb_C2_liver_heatmap_meta.pdf", width = 8.3, height = 8.3)
# ggsave("plots/anno/msigdb_C2_liver_heatmap_meta.png", width = 8.3, height = 8.3)

saveRDS(patient.seurat, "data/patient_seurat_hvgs_meta_sigs.rds")

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk //////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

patient.seurat <- readRDS("data/patient_seurat_hvgs_meta.rds")
pam.cols <- readRDS("data/pam_cluster_cols.rds")
meta.cols <- readRDS("data/meta_cluster_cols.rds")

pseudo.pam <- AggregateExpression(patient.seurat, assays = "RNA", return.seurat = T,
                                  group.by = c("PAM.Cluster"))
saveRDS(pseudo.pam, "data/patient_pseudobulk_pam.rds")

pseudo.meta <- AggregateExpression(patient.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Meta.Cluster"))
saveRDS(pseudo.meta, "data/patient_pseudobulk_meta.rds")

## MSigDB Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
H.list <- lapply(unique(H.df$gs_name), function(x){
  H.df[H.df$gs_name==x, "gene_symbol"]
})
names(H.list) <- unique(H.df$gs_name)
H.list
H.names <- paste(names(H.list), ".Sig", sep = "")
H.names <- gsub("_", "-", H.names)

pseudo.pam <- AddModuleScore(pseudo.pam, features = H.list, assay = "RNA", seed = 12345, 
                             name = H.names)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.pam, vars = H.names)))
pseudo.pam$PAM.Cluster <- gsub("^g", "", pseudo.pam$PAM.Cluster)
pseudo.pam$PAM.Cluster <- factor(pseudo.pam$PAM.Cluster, levels = 1:11)
DoHeatmap(pseudo.pam, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_pam_msigdb_H_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_pam_msigdb_H_heatmap.png", width = 8.3, height = 8.3)

pam.h.mat <- as.matrix(pseudo.pam@assays$MSigDB_H@data)
pam.anno.df <- pseudo.pam@meta.data["PAM.Cluster"]
gt <- pheatmap(pam.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = pam.anno.df,
               annotation_colors = list(PAM.Cluster = pam.cols))$gtable
ggsave("plots/anno/pseudo_pam_msigdb_H_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_pam_msigdb_H_pheatmap.png", plot = gt)

pseudo.meta <- AddModuleScore(pseudo.meta, features = H.list, assay = "RNA", seed = 12345, 
                              name = H.names)
colnames(pseudo.meta@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.meta@meta.data))
pseudo.meta[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.meta, vars = H.names)))
pseudo.meta$Meta.Cluster <- gsub("^g", "", pseudo.meta$Meta.Cluster)
pseudo.meta$Meta.Cluster <- factor(pseudo.meta$Meta.Cluster, levels = c(1:11, "NA"))
DoHeatmap(pseudo.meta, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "Meta.Cluster", size = 4, draw.lines = FALSE, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_meta_msigdb_H_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_meta_msigdb_H_heatmap.png", width = 8.3, height = 8.3)

meta.h.mat <- as.matrix(pseudo.meta@assays$MSigDB_H@data)
meta.anno.df <- pseudo.meta@meta.data["Meta.Cluster"]
gt <- pheatmap(meta.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = meta.anno.df,
               annotation_colors = list(Meta.Cluster = meta.cols))$gtable
ggsave("plots/anno/pseudo_meta_msigdb_H_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_meta_msigdb_H_pheatmap.png", plot = gt)

## MSigDB C2 gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
liver.match <- c("HEPATOCYTE", "HEPATOBLAST", "LIVER", "EMT", "EPITHELIAL_MESENCHYMAL",
                 "TGFB", "TGF_BETA", "WNT", "NOTCH", "HEDGEHOG", "SHH")
liver.gs <- C2.df[grepl(paste(liver.match, collapse = "|"), C2.df$gs_name), ]
liver.list <- lapply(unique(liver.gs$gs_name), function(x){
  liver.gs[liver.gs$gs_name==x, "gene_symbol"]
})
names(liver.list) <- unique(liver.gs$gs_name)
liver.list
liver.names <- paste(names(liver.list), ".Sig", sep = "")
liver.names <- gsub("_", "-", liver.names)
abbr_sig <- function(x) {
  core <- sub("\\.Sig$", "", x)
  parts <- strsplit(core, "-", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(parts)
  paste(
    c(parts[1], substr(parts[-1], 1, 1)),
    collapse = "-"
  )
}
liver.names.abbr <- sapply(liver.names, abbr_sig, USE.NAMES = TRUE)

pseudo.pam <- AddModuleScore(pseudo.pam, features = liver.list, assay = "RNA", seed = 12345, 
                             name = liver.names)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.pam, vars = liver.names)))
DoHeatmap(pseudo.pam, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE, group.colors = pam.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  scale_y_discrete(labels = liver.names.abbr) +
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 2)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_pam_msigdb_C2_liver_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_pam_msigdb_C2_liver_heatmap.png", width = 8.3, height = 8.3)
pam.c2.mat <- as.matrix(pseudo.pam@assays$MSigDB_C2@data)
pam.reactome.mat <- pam.c2.mat[which(grepl("REACTOME", rownames(pam.c2.mat))), ]
gt <- pheatmap(pam.reactome.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = pam.anno.df,
               annotation_colors = list(PAM.Cluster = pam.cols))$gtable
ggsave("plots/anno/pseudo_pam_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_pam_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 8.3, height = 8.3)
pam.cairo.mat <- pam.c2.mat[which(grepl("CAIRO", rownames(pam.c2.mat))), ]
gt <- pheatmap(pam.cairo.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = pam.anno.df,
               annotation_colors = list(PAM.Cluster = pam.cols))$gtable
ggsave("plots/anno/pseudo_pam_msigdb_C2_cairo_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_pam_msigdb_C2_cairo_pheatmap.png", plot = gt)

pseudo.meta <- AddModuleScore(pseudo.meta, features = liver.list, assay = "RNA", seed = 12345, 
                              name = liver.names)
colnames(pseudo.meta@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.meta@meta.data))
pseudo.meta[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.meta, vars = liver.names)))
DoHeatmap(pseudo.meta, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "Meta.Cluster", size = 4, draw.lines = FALSE, group.colors = meta.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  scale_y_discrete(labels = liver.names.abbr) +
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 2)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_meta_msigdb_C2_liver_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_meta_msigdb_C2_liver_heatmap.png", width = 8.3, height = 8.3)
meta.c2.mat <- as.matrix(pseudo.meta@assays$MSigDB_C2@data)
meta.reactome.mat <- meta.c2.mat[which(grepl("REACTOME", rownames(meta.c2.mat))), ]
gt <- pheatmap(meta.reactome.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = meta.anno.df,
               annotation_colors = list(Meta.Cluster = meta.cols))$gtable
ggsave("plots/anno/pseudo_meta_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_meta_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 8.3, height = 8.3)
meta.cairo.mat <- meta.c2.mat[which(grepl("CAIRO", rownames(meta.c2.mat))), ]
gt <- pheatmap(meta.cairo.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = meta.anno.df,
               annotation_colors = list(Meta.Cluster = meta.cols))$gtable
ggsave("plots/anno/pseudo_meta_msigdb_C2_cairo_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_meta_msigdb_C2_cairo_pheatmap.png", plot = gt)

saveRDS(pseudo.pam, "data/patient_pseudobulk_pam_sigs.rds")
saveRDS(pseudo.meta, "data/patient_pseudobulk_meta_sigs.rds")

## Liver differentiation markers
## (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/; https://pmc.ncbi.nlm.nih.gov/articles/PMC11114060/)
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
custom %in% rownames(pseudo.pam)
custom %in% rownames(pseudo.meta)

pam.custom.mat <- as.matrix(pseudo.pam@assays$RNA$scale.data[custom, ])
gt <- pheatmap(pam.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = pam.anno.df,
               annotation_colors = list(PAM.Cluster = pam.cols))$gtable
ggsave("plots/anno/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_pam_custom_pheatmap.png", plot = gt)

meta.custom.mat <- as.matrix(pseudo.meta@assays$RNA$scale.data[custom, ])
gt <- pheatmap(meta.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = meta.anno.df,
               annotation_colors = list(Meta.Cluster = meta.cols))$gtable
ggsave("plots/anno/pseudo_meta_custom_pheatmap.pdf", plot = gt)
ggsave("plots/anno/pseudo_meta_custom_pheatmap.png", plot = gt)

## [ Cell annotation ] ----

patient.seurat <- readRDS("data/patient_seurat_hvgs_meta.rds")
patient.sce <- as.SingleCellExperiment(patient.seurat)
hpca.se <- celldex::HumanPrimaryCellAtlasData()

library(SingleR)
pred.patient <- SingleR(test = patient.sce, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)
saveRDS(pred.patient, "data/patient_singler_prediction.rds")

table(pred.patient$labels)
## Astrocyte               B_cell         Chondrocytes                   DC Embryonic_stem_cells    Endothelial_cells     Epithelial_cells 
##      6669                   26                 4891                  133                    8                 1781                   15 
## Erythroblast          Fibroblasts                  GMP          Hepatocytes           HSC_-G-CSF            HSC_CD34+            iPS_cells 
##            3                 1541                    3                44920                    1                    1                  294 
## Keratinocytes           Macrophage             Monocyte                  MSC Neuroepithelial_cell              Neurons              NK_cell 
##             1                 1295                  573                  400                  322                 3634                   65 
## Osteoblasts            Platelets  Smooth_muscle_cells              T_cells    Tissue_stem_cells 
##        2376                   41                 6012                  145                  552 

pdf("plots/anno/patient_singler_heatmap.pdf", width = 8.3, height = 8.3)
png("plots/anno/patient_singler_heatmap.png", units = "in", res = 200, width = 8.3, height = 8.3)
plotScoreHeatmap(pred.patient)
dev.off()

plotDeltaDistribution(pred.patient, ncol = 3)
ggsave("plots/anno/patient_singler_delta.pdf", width = 5.8, height = 11.7)
ggsave("plots/anno/patient_singler_delta.png", width = 5.8, height = 11.7)

colnames(pred.patient) <- paste("SingleR", colnames(pred.patient), sep = ".")
colData(patient.sce) <- merge(colData(patient.sce), pred.patient, by = "row.names", all.x = TRUE)
identical(colnames(patient.seurat), rownames(pred.patient)) ## TRUE
patient.seurat <- AddMetaData(object = patient.seurat,
                              metadata = as.data.frame(pred.patient))
saveRDS(patient.seurat, "data/patient_seurat_singler.rds")

anno.cols <- c("#D62728", "#FF9896", "#1F77B4", "#FFBB78",
               "#FF7F0E", "#17BECF", "#9467BD", "#2CA02C",
               "#E377C2", "#AEC7E8", "#843C39", "#BCBD22",
               "#9EDAE5", "#393B79", "#C5B0D5", "#8C564B",
               "#7F7F7F", "#F7B6D2", "#7B4173", "#C7C7C7",
               "#DBDB8D", "#C49C94", "#98DF8A", "#637939",
               "#8C6D31", "#525252")
names(anno.cols) <- unique(patient.seurat$SingleR.labels)
saveRDS(anno.cols, "data/singler_cols.rds")

DimPlot(patient.seurat, group.by = "SingleR.labels", order = TRUE) +
  scale_colour_manual(values = anno.cols) +
  umap.theme() + labs(title = "SingleR annotation") +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave("plots/umaps/patient_umap_singler.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/patient_umap_singler.png", width = 8.3, height = 5.8)

pam.cols <- readRDS("data/pam_cluster_cols.rds")
meta.cols <- readRDS("data/meta_cluster_cols.rds")

pam.df <- as.data.frame.matrix(table(patient.seurat$PAM.Cluster, patient.seurat$SingleR.labels))
pam.pheatmap.cols <- list(PAM = pam.cols,
                          CellType = anno.cols)
pam.anno.row <- data.frame(PAM = factor(rownames(pam.df)))
rownames(pam.anno.row) <- rownames(pam.df)
pam.anno.col <- data.frame(CellType = factor(colnames(pam.df)))
rownames(pam.anno.col) <- colnames(pam.df)
pam.log.df <- log1p(pam.df)
gt <- pheatmap(pam.log.df,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               main = "log(Cell Number + 1)",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_row = pam.anno.row,
               annotation_col = pam.anno.col,
               annotation_colors = pam.pheatmap.cols)
ggsave(plot = gt, "plots/anno/patient_heatmap_pam_cell.pdf", width = 11.7, height = 16.5)
ggsave(plot = gt, "plots/anno/patient_heatmap_pam_cell.png", width = 11.7, height = 16.5)

meta.df <- as.data.frame.matrix(table(patient.seurat$Meta.Cluster, patient.seurat$SingleR.labels))
meta.pheatmap.cols <- list(Meta = meta.cols,
                           CellType = anno.cols)
meta.anno.row <- data.frame(Meta = factor(rownames(meta.df)))
rownames(meta.anno.row) <- rownames(meta.df)
meta.anno.col <- data.frame(CellType = factor(colnames(meta.df)))
rownames(meta.anno.col) <- colnames(meta.df)
meta.log.df <- log1p(meta.df)
gt <- pheatmap(meta.log.df,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               main = "log(Cell Number + 1)",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_row = meta.anno.row,
               annotation_col = meta.anno.col,
               annotation_colors = meta.pheatmap.cols)
ggsave(plot = gt, "plots/anno/patient_heatmap_meta_cell.pdf", width = 11.7, height = 16.5)
ggsave(plot = gt, "plots/anno/patient_heatmap_meta_cell.png", width = 11.7, height = 16.5)

## [ InferCNV ] ----

## Gene reference
patient.seurat <- readRDS("data/patient_seurat_singler.rds")
ens.bm <- read.csv("data/patient_ensembl_biomart.csv", row.names = NULL)
ens.ref <- ens.bm[ens.bm$external_gene_name %in% rownames(patient.seurat), ]
ens.ref <- ens.ref[!duplicated(ens.ref$external_gene_name), ]
gene.order <- ens.ref[c("external_gene_name", "chromosome_name", "start_position", "end_position")]
saveRDS(gene.order, "data/gene_order_ref.rds")
write.table(gene.order, "data/gene_order_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

anno.df <- patient.seurat@meta.data[c("CellID", "SingleR.labels")]
immune.cells <- c("B_cell", "NK_cell", "T_cells", "Monocyte", "Macrophage", "DC")
anno.df$SingleR.labels[anno.df$SingleR.labels %in% immune.cells] <- "Immune_ref"
saveRDS(anno.df, "data/cell_anno_ref.rds")
write.table(anno.df, "data/cell_anno_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

## Run InferCNV object on HPC cluster
# ## Create InferCNV object
# library(infercnv)
# infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = patient.seurat@assays$RNA$counts,
#                                      annotations_file = "data/cell_anno_ref.txt",
#                                      delim = "\t",
#                                      gene_order_file = "data/gene_order_ref.txt",
#                                      ref_group_names = "Immune_ref")
# saveRDS(infercnv.obj, "data/infercnv_input.rds")
# options(scipen = 999)
# ## Run InferCNV
# infercnv.obj <- infercnv::run(infercnv.obj,
#                               cutoff = 0.1, ## Use 1 for smart-seq, 0.1 for 10x-genomics
#                               out_dir = "infercnv",  
#                               cluster_by_groups = F,
#                               denoise = T,
#                               HMM = F) ## Was failing hspike modelling
# saveRDS(infercnv.obj, "data/infercnv_output.rds")

infercnv.obj <- readRDS("data/infercnv/15_tumor_subclusters.leiden.infercnv_obj")

expr.mat <- infercnv.obj@expr.data
nrow(expr.mat) ## 8761
ncol(expr.mat) ## 75702

## Calculate CNV scores
cnv.scores <- colMeans(abs(expr.mat))
cnv.scores.df <- data.frame(cell = names(cnv.scores),
                            cnv_score = as.numeric(cnv.scores),
                            stringsAsFactors = FALSE)
write.table(cnv.scores.df,
            file = "data/infercnv/cnv_scores_per_cell.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

gene.anno <- readRDS("data/gene_order_ref.rds")
gene.anno$gene_idx <- match(gene.anno$external_gene_name, rownames(expr.mat))
gene.anno <- gene.anno[!is.na(gene.anno$gene_idx), ]
gene.anno$chr <- paste0("chr", gene.anno$chr)

## Calculate scores for specific chromosomes
calculate_score <- function(chr_name, expr_matrix, gene_info) {
  ## Handle whole chromosomes (e.g., "chr7")
  if (!grepl("[pq]$", chr_name)) {
    gene_indices <- gene_info$gene_idx[gene_info$chr == chr_name]
  } else {
    ## Handle chromosome arms (e.g., "chr17q")
    chr_base <- gsub("[pq]$", "", chr_name)
    arm <- substr(chr_name, nchar(chr_name), nchar(chr_name))
    ## Get centromere position (median of gene positions)
    chr_genes <- gene_info[gene_info$chr == chr_base, ]
    if (nrow(chr_genes) == 0) {
      return(NULL)
    }
    centromere_pos <- median(chr_genes$start)
    if (arm == "p") {
      gene_indices <- chr_genes$gene_idx[chr_genes$start < centromere_pos]
    } else {
      gene_indices <- chr_genes$gene_idx[chr_genes$start >= centromere_pos]
    }
  }
  if (length(gene_indices) < 10) {
    message("  WARNING: Skipping ", chr_name, " (only ", length(gene_indices), " genes)")
    return(NULL)
  }
  ## Extract expression for this chromosome/arm
  chr_expr <- expr_matrix[gene_indices, , drop = FALSE]
  ## Calculate score per cell: mean absolute deviation
  chr_scores <- colMeans(abs(chr_expr))
  message(
    "  ", chr_name, ": ", length(gene_indices), " genes, mean score = ",
    round(mean(chr_scores), 4)
  )
  return(chr_scores)
}
chr.arms <- c("chr1q", "chr2", "chr8", "chr20", ## Frequent gains
              "chr1p", "chr4", "chr11q", "chr18q", "chr5p", "chr16q") ## Frequent losses
## Tomlinson and Kappler, (2012); doi:10.1002/pbc.24213
## Barros et al., (2021); doi: 10.3389/fonc.2021.741526
## 1p also gain? Wu et al., (2013); doi: 10.1007/s12072-012-9350-y
chr.scores.list <- list()
for (arm in chr.arms) {
  scores <- calculate_score(arm, expr.mat, gene.anno)
  if (!is.null(scores)) {
    chr.scores.list[[arm]] <- scores
  }
}

## Add CNV scores
patient.seurat <- readRDS("data/patient_seurat_singler.rds")
patient.seurat$infercnv_score <- NA
common.cells <- intersect(colnames(patient.seurat), cnv.scores.df$cell) ## 75702
patient.seurat@meta.data[common.cells, "infercnv_score"] <- cnv.scores.df$cnv_score[match(common.cells, cnv.scores.df$cell)]

## Add chromosome-specific scores
for (arm in names(chr.scores.list)) {
  col_name <- paste0(arm, "_score")
  patient.seurat@meta.data[[col_name]] <- NA
  chr.common.cells <- intersect(colnames(patient.seurat), names(chr.scores.list[[arm]]))
  patient.seurat@meta.data[chr.common.cells, col_name] <- chr.scores.list[[arm]][chr.common.cells]
}

saveRDS(patient.seurat, "data/patient_seurat_infercnv.rds")

## CNV score histograms
ggplot(patient.seurat@meta.data, aes(x = infercnv_score)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001,
                 fill = "#1F77B4", colour = "#1F77B4", alpha = 0.7) +
  geom_density(colour = "black", size = 1.2) +
  labs(title = "InferCNV Score Distribution", x = "InferCNV score") +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("plots/infercnv/infercnv_score_histogram.pdf", width = 8.3, height = 5.8)
ggsave("plots/infercnv/infercnv_score_histogram.png", width = 8.3, height = 5.8)

chr.score <- colnames(patient.seurat@meta.data)
chr.score <- chr.score[grep("^chr.*_score$", chr.score)]
make_title <- function(x) {
  x |> 
    gsub("_", " ", x = _) |> 
    tools::toTitleCase()
}
for (score in chr.score) {
  name <- score
  p <- ggplot(patient.seurat@meta.data, aes(x = .data[[score]])) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01,
                   fill = "#1F77B4", colour = "#1F77B4", alpha = 0.7) +
    geom_density(colour = "black", size = 1.2) +
    labs(title = paste(make_title(score), "Distribution"), x = "InferCNV score") +
    theme_bw() +
    theme(text = element_text(size = 15))
  ggsave(plot = p, paste0("plots/infercnv/", score, "_histogram.pdf"), width = 8.3, height = 5.8)
  ggsave(plot = p, paste0("plots/infercnv/", score, "_histogram.png"), width = 8.3, height = 5.8)
}

## CNV score UMAPs
FeaturePlot(patient.seurat, features = "infercnv_score") +
  scale_colour_gradientn(colours = pals::coolwarm(100)) +
  labs(title = "InferCNV Score", colour = "") +
  umap.theme() +
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))
ggsave("plots/infercnv/infercnv_score_umap.pdf", width = 8.3, height = 5.8)
ggsave("plots/infercnv/infercnv_score_umap.png", width = 8.3, height = 5.8)

for (score in chr.score) {
  FeaturePlot(patient.seurat, features = score) +
    scale_colour_gradientn(colours = pals::coolwarm(100)) +
    labs(title = paste(make_title(score)), colour = "") +
    umap.theme() +
    theme(text = element_text(size = 15),
          legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/infercnv/", score, "_umap.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/infercnv/", score, "_umap.png"), width = 8.3, height = 5.8)
}

## [ Meta-clustering ] ----

## Meta-clustering using Jaccard similarity and PAM clustering had limited success
## in finding similarity between clusters in different samples
## Try manual meta-clustering by comparing the expression similarity between clusters in different samples

patient.list <- readRDS("data/patient_seurat_split_list_pam.rds")
patient.seurat <- readRDS("data/patient_seurat_singler.rds")

## Assign inital clusters from separated objects into combined
patient.list$`12-00122_F1`$initial_clusters <- paste0("P1_F1_", patient.list$`12-00122_F1`$seurat_clusters.0.2)
patient.list$`12-00122_F2`$initial_clusters <- paste0("P1_F2_", patient.list$`12-00122_F2`$seurat_clusters.0.2)
patient.list$`12-01026_F1`$initial_clusters <- paste0("P2_F1_", patient.list$`12-01026_F1`$seurat_clusters.0.2)
patient.list$`12-01026_F2`$initial_clusters <- paste0("P2_F2_", patient.list$`12-01026_F2`$seurat_clusters.0.4)
patient.list$`16-00875_F1`$initial_clusters <- paste0("P3_F1_", patient.list$`16-00875_F1`$seurat_clusters.0.2)
patient.list$`16-00875_F2`$initial_clusters <- paste0("P3_F2_", patient.list$`16-00875_F2`$seurat_clusters.0.2)

merged <- Reduce(function(x, y) merge(x, y), patient.list)
merged.meta <- merged@meta.data
length(rownames(merged.meta) %in% rownames(patient.seurat@meta.data)) ## 75702
merged.meta$CellID <- rownames(merged.meta)
merged.ordered <- merged.meta[match(patient.seurat$CellID, merged.meta$CellID), ]
identical(rownames(merged.ordered), rownames(patient.seurat@meta.data)) ## TRUE
patient.seurat$Init.Cluster <- merged.ordered$initial_clusters

clusters <- unique(patient.seurat$Init.Cluster)
parts <- do.call(rbind, strsplit(clusters, "_"))
p.num <- as.numeric(sub("P", "", parts[,1]))
f.num    <- as.numeric(sub("F", "", parts[,2]))
c.num <- as.numeric(parts[,3])
order.num <- order(p.num, f.num, c.num)
order.levels <- clusters[order.num]
patient.seurat$Init.Cluster <- factor(patient.seurat$Init.Cluster, levels = order.levels)
saveRDS(patient.seurat, "data/patient_seurat_hvgs_init.rds")

init.cols <- hues::iwanthue(length(unique(patient.seurat$Init.Cluster)))
names(init.cols) <- order.levels
saveRDS(init.cols, "data/initial_cluster_cols.rds")

for (patient in names(patient.list)) {
  DimPlot(patient.list[[patient]], group.by = "initial_clusters", order = TRUE) +
    scale_colour_manual(values = init.cols) +
    umap.theme() + labs(title = paste(patient, "Initial Clusters"))
  ggsave(paste0("plots/umaps/umap_init_", patient,".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_init_", patient,".png"), width = 8.3, height = 5.8)
}

DimPlot(patient.seurat, group.by = "Init.Cluster", order = TRUE) +
  scale_colour_manual(values = init.cols) +
  umap.theme() + labs(title = "Initial Clusters")
ggsave("plots/umaps/umap_harmony_init.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_init.png", width = 8.3, height = 5.8)

pam.df <- patient.seurat@meta.data %>%
  count(Condition, PAM.Cluster) %>%
  group_by(Condition) %>%                 
  mutate(prop = n / sum(n))  
pam.cols <- readRDS("data/pam_cluster_cols.rds")
ggplot(pam.df, aes(x = Condition, y = prop, fill = PAM.Cluster)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Sample") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = pam.cols) +
  guides(fill = guide_legend(title = "PAM clusters")) +
  theme_bw() +
  theme(text = element_text(size = 15), aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/anno/proportion_pam.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/proportion_pam.png", width = 8.3, height = 5.8)

meta.df <- patient.seurat@meta.data %>%
  count(Condition, Meta.Cluster) %>%
  group_by(Condition) %>%                 
  mutate(prop = n / sum(n))  
meta.cols <- readRDS("data/meta_cluster_cols.rds")
ggplot(meta.df, aes(x = Condition, y = prop, fill = Meta.Cluster)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Sample") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = meta.cols) +
  guides(fill = guide_legend(title = "Meta clusters")) +
  theme_bw() +
  theme(text = element_text(size = 15), aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/anno/proportion_meta.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/proportion_meta.png", width = 8.3, height = 5.8)

init.df <- patient.seurat@meta.data %>%
  count(Condition, Init.Cluster) %>%
  group_by(Condition) %>%                 
  mutate(prop = n / sum(n))  
ggplot(init.df, aes(x = Condition, y = prop, fill = Init.Cluster)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Sample") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = init.cols) +
  guides(fill = guide_legend(title = "Initial clusters", ncol = 4)) +
  theme_bw() +
  theme(text = element_text(size = 15), aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/anno/proportion_init.pdf", width = 11.7, height = 5.8)
ggsave("plots/anno/proportion_init.png", width = 11.7, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk //////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

pseudo.init <- AggregateExpression(patient.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Init.Cluster"))
saveRDS(pseudo.init, "data/patient_pseudobulk_init.rds")

## MSigDB Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
H.list <- lapply(unique(H.df$gs_name), function(x){
  H.df[H.df$gs_name==x, "gene_symbol"]
})
names(H.list) <- unique(H.df$gs_name)
H.list
H.names <- paste(names(H.list), ".Sig", sep = "")
H.names <- gsub("_", "-", H.names)

pseudo.init <- AddModuleScore(pseudo.init, features = H.list, assay = "RNA", seed = 12345, 
                              name = H.names)
colnames(pseudo.init@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.init@meta.data))
pseudo.init[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.init, vars = H.names)))
init.levels <- gsub("_", "-", levels(patient.seurat$Init.Cluster))
pseudo.init$Init.Cluster <- factor(pseudo.init$Init.Cluster, levels = init.levels)
names(init.cols) <- init.levels
DoHeatmap(pseudo.init, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "Init.Cluster", size = 2, draw.lines = FALSE, group.colors = init.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_init_msigdb_H_heatmap.pdf", width = 11.7, height = 11.7)
ggsave("plots/anno/pseudo_init_msigdb_H_heatmap.png", width = 11.7, height = 11.7)

init.h.mat <- as.matrix(pseudo.init@assays$MSigDB_H@data)
init.anno.df <- pseudo.init@meta.data["Init.Cluster"]
gt <- pheatmap(init.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/anno/pseudo_init_msigdb_H_pheatmap.pdf", plot = gt, width = 11.7, height = 23.4)
ggsave("plots/anno/pseudo_init_msigdb_H_pheatmap.png", plot = gt, width = 11.7, height = 23.4)

## MSigDB C2 gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
liver.match <- c("HEPATOCYTE", "HEPATOBLAST", "LIVER", "EMT", "EPITHELIAL_MESENCHYMAL",
                 "TGFB", "TGF_BETA", "WNT", "NOTCH", "HEDGEHOG", "SHH")
liver.gs <- C2.df[grepl(paste(liver.match, collapse = "|"), C2.df$gs_name), ]
liver.list <- lapply(unique(liver.gs$gs_name), function(x){
  liver.gs[liver.gs$gs_name==x, "gene_symbol"]
})
names(liver.list) <- unique(liver.gs$gs_name)
liver.list
liver.names <- paste(names(liver.list), ".Sig", sep = "")
liver.names <- gsub("_", "-", liver.names)
abbr_sig <- function(x) {
  core <- sub("\\.Sig$", "", x)
  parts <- strsplit(core, "-", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(parts)
  paste(
    c(parts[1], substr(parts[-1], 1, 1)),
    collapse = "-"
  )
}
liver.names.abbr <- sapply(liver.names, abbr_sig, USE.NAMES = TRUE)

pseudo.init <- AddModuleScore(pseudo.init, features = liver.list, assay = "RNA", seed = 12345, 
                              name = liver.names)
colnames(pseudo.init@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.init@meta.data))
pseudo.init[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.init, vars = liver.names)))
DoHeatmap(pseudo.init, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "Init.Cluster", size = 2, draw.lines = FALSE, group.colors = init.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  scale_y_discrete(labels = liver.names.abbr) +
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 2)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_init_msigdb_C2_liver_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_msigdb_C2_liver_heatmap.png", width = 8.3, height = 8.3)
init.c2.mat <- as.matrix(pseudo.init@assays$MSigDB_C2@data)
init.reactome.mat <- init.c2.mat[which(grepl("REACTOME", rownames(init.c2.mat))), ]
gt <- pheatmap(init.reactome.mat,
               border_color = NA,
               cellwidth = 5, cellheight = 5,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/anno/pseudo_init_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 11.7, height = 23.4)
ggsave("plots/anno/pseudo_init_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 11.7, height = 23.4)
init.cairo.mat <- init.c2.mat[which(grepl("CAIRO", rownames(init.c2.mat))), ]
gt <- pheatmap(init.cairo.mat,
               border_color = NA,
               cellwidth = 10, cellheight = 10,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/anno/pseudo_init_msigdb_C2_cairo_pheatmap.pdf", plot = gt, width = 16.5, height = 23.4)
ggsave("plots/anno/pseudo_init_msigdb_C2_cairo_pheatmap.png", plot = gt, width = 16.5, height = 23.4)

## HB signatures from literature
sig.df <- as.data.frame(readxl::read_xlsx(path = "data/hb_sigs.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")
sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})
names(sigs) <- unique(sig.df$Signature)
sigs
sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)
pseudo.init <- AddModuleScore(pseudo.init, features = sigs, assay = "RNA", seed = 12345, 
                              name = sig.names)
colnames(pseudo.init@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.init@meta.data))
pseudo.init[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.init, vars = sig.names)))
DoHeatmap(pseudo.init, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "Init.Cluster", size = 2, draw.lines = FALSE, group.colors = init.cols) +
  scale_fill_gradientn(colors = c("#F2F0EF", "#C23B22")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(0.8, "cm"),
        aspect.ratio = 1, axis.text.y = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/anno/pseudo_init_hb_sigs_heatmap.pdf", width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_hb_sigs_heatmap.png", width = 8.3, height = 8.3)

init.hb.mat <- as.matrix(pseudo.init@assays$HB_sigs@data)
gt <- pheatmap(init.hb.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/anno/pseudo_init_hb_sigs_pheatmap.pdf", plot = gt, width = 16.5, height = 23.4)
ggsave("plots/anno/pseudo_init_hb_sigs_pheatmap.png", plot = gt, width = 16.5, height = 23.4)

saveRDS(pseudo.init, "data/patient_pseudobulk_init_sigs.rds")

## Liver differentiation markers
## (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/; https://pmc.ncbi.nlm.nih.gov/articles/PMC11114060/)
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
custom %in% rownames(pseudo.init)

init.custom.mat <- as.matrix(pseudo.init@assays$RNA$scale.data[custom, ])
gt <- pheatmap(init.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/anno/pseudo_init_custom_pheatmap.pdf", plot = gt, width = 11.7, height = 23.4)
ggsave("plots/anno/pseudo_init_custom_pheatmap.png", plot = gt, width = 11.7, height = 23.4)

## Cell type annotation
init.df <- as.data.frame.matrix(table(patient.seurat$Init.Cluster, patient.seurat$SingleR.labels))
init.cols <- readRDS("data/initial_cluster_cols.rds")
anno.cols <- readRDS("data/singler_cols.rds")
init.pheatmap.cols <- list(Initial = init.cols,
                           CellType = anno.cols)
init.anno.row <- data.frame(Initial = factor(rownames(init.df)))
rownames(init.anno.row) <- rownames(init.df)
init.anno.col <- data.frame(CellType = factor(colnames(init.df)))
rownames(init.anno.col) <- colnames(init.df)
init.log.df <- log1p(init.df)
gt <- pheatmap(init.log.df,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               main = "log(Cell Number + 1)",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_row = init.anno.row,
               annotation_col = init.anno.col,
               annotation_colors = init.pheatmap.cols)
ggsave(plot = gt, "plots/anno/patient_heatmap_init_cell.pdf", width = 11.7, height = 23.4)
ggsave(plot = gt, "plots/anno/patient_heatmap_init_cell.png", width = 11.7, height = 23.4)

## Heatmaps by patient
# pseudo.init <- readRDS("data/patient_pseudobulk_init_sigs.rds")
# init.anno.df <- pseudo.init@meta.data["Init.Cluster"]
# init.cols <- readRDS("data/initial_cluster_cols.rds")
# names(init.cols) <- gsub("_", "-", names(init.cols))
## Then make init.custom.mat

sample.id <- sub("-.*$", "", colnames(init.custom.mat))
for (p in unique(sample.id)) {
  cols <- which(sample.id == p)
  mat.sub  <- init.custom.mat[, cols, drop = FALSE]
  anno.sub <- init.anno.df[colnames(mat.sub), , drop = FALSE]
  used.clusters <- unique(anno.sub$Init.Cluster)
  anno.cols.sub <- list(Init.Cluster = init.cols[used.clusters])
  gt <- pheatmap(mat.sub,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 annotation_col = anno.sub,
                 annotation_colors = anno.cols.sub,
                 main = p)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_init_custom_pheatmap_", p, ".pdf"), plot = gt)
  ggsave(filename = paste0("plots/anno/pseudo_init_custom_pheatmap_", p, ".png"), plot = gt)
}

## [ Patient UMAPs ] ----

## Split data
patient.seurat <- readRDS("data/patient_seurat_hvgs_init.rds")
obj.list <- SplitObject(patient.seurat, split.by = "Patient_ID")
obj.list <- lapply(obj.list, FindVariableFeatures)
rm.genes <- c(grep("^MT-", rownames(patient.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(patient.seurat), value = TRUE),
              grep("[\\.]", rownames(patient.seurat), value = TRUE),
              grep("^LINC", rownames(patient.seurat), value = TRUE),
              c("MALAT1"))
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(patient.seurat, obj.list, obj)

output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 100,
         features = VariableFeatures(obj))
})
output <- lapply(output, RunUMAP, reduction = "pca", dims = 1:30)
output <- lapply(output, FindNeighbors, reduction = "pca", dims = 1:30)

saveRDS(output, "data/patient_seurat_split_id_list.rds")

patient.list <- readRDS("data/patient_seurat_split_id_list.rds")
for (patient in names(patient.list)) {
  DimPlot(patient.list[[patient]], group.by = "Condition", order = TRUE) +
    scale_colour_manual(values = condition.cols) +
    umap.theme() + labs(title = paste(patient, " - Conditions"))
  ggsave(paste0("plots/umaps/umap_condition_", patient, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_condition_", patient, ".png"), width = 8.3, height = 5.8)
}
init.cols <- readRDS("data/initial_cluster_cols.rds")
for (patient in names(patient.list)) {
  DimPlot(patient.list[[patient]], group.by = "Init.Cluster", order = TRUE) +
    scale_colour_manual(values = init.cols) +
    umap.theme() + labs(title = paste(patient, " - Initial Cluster"))
  ggsave(paste0("plots/umaps/umap_init_", patient, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_init_", patient, ".png"), width = 8.3, height = 5.8)
}

## Patient 1
patient.list$`12-00122`$Final.Cluster <- recode(as.character(patient.list$`12-00122`$Init.Cluster),
                                                "P1_F1_0" = "Late hepatoblast",
                                                "P1_F1_1" = "Early hepatoblast",
                                                "P1_F1_2" = "TAMs",
                                                "P1_F1_3" = "Infiltrated",
                                                "P1_F1_4" = "Infiltrated",
                                                "P1_F1_5" = "Endothelial",
                                                "P1_F1_6" = "Infiltrated",
                                                "P1_F1_7" = "Stroma",
                                                "P1_F1_8" = "Stroma",
                                                "P1_F2_0" = "Late hepatoblast")
## Patient 2
patient.list$`12-01026`$Final.Cluster <- recode(as.character(patient.list$`12-01026`$Init.Cluster),
                                                "P2_F1_0" = "Late hepatoblast",
                                                "P2_F1_1" = "Hepatocytic",
                                                "P2_F1_2" = "Hepatocytic",
                                                "P2_F1_3" = "Late hepatoblast",
                                                "P2_F1_4" = "Late hepatoblast",
                                                "P2_F1_5" = "Stroma",
                                                "P2_F1_6" = "Infiltrated",
                                                "P2_F1_7" = "Endothelial",
                                                "P2_F1_8" = "Infiltrated",
                                                "P2_F1_9" = "Hepatocytic",
                                                "P2_F2_0" = "Hepatocyte 3",
                                                "P2_F2_1" = "Hepatocyte 1",
                                                "P2_F2_2" = "Mature hepatocyte",
                                                "P2_F2_3" = "Mature hepatocyte",
                                                "P2_F2_4" = "Hepatocyte 2",
                                                "P2_F2_5" = "Hepatocyte 1",
                                                "P2_F2_6" = "Late hepatoblast",
                                                "P2_F2_7" = "Infiltrated",
                                                "P2_F2_8" = "Endothelial")
## Patient 3
patient.list$`16-00875`$Final.Cluster <- recode(as.character(patient.list$`16-00875`$Init.Cluster),
                                                "P3_F1_0" = "TSC-like TGFb",
                                                "P3_F1_1" = "TSC-like TGFb",
                                                "P3_F1_2" = "TSC-like EMT",
                                                "P3_F1_3" = "TSC-like TGFb",
                                                "P3_F1_4" = "TSC-like EMT",
                                                "P3_F1_5" = "Endothelial",
                                                "P3_F1_6" = "MSC-like",
                                                "P3_F1_7" = "HSC-like",
                                                "P3_F1_8" = "Stroma", ## Unclear
                                                "P3_F1_9" = "Stroma", ## Could be Infiltrated or TSC-like EMT?
                                                "P3_F1_10" = "Hepatocyte",
                                                "P3_F1_11" = "MSC-like",
                                                "P3_F2_0" = "Astrocyte 1",
                                                "P3_F2_1" = "Astrocyte 2",
                                                "P3_F2_2" = "Astrocyte 1",
                                                "P3_F2_3" = "Astrocyte 1",
                                                "P3_F2_4" = "Astrocyte 2",
                                                "P3_F2_5" = "Astrocyte 1",
                                                "P3_F2_6" = "Astrocyte 2",
                                                "P3_F2_7" = "Astrocyte 1",
                                                "P3_F2_8" = "Infiltrated",
                                                "P3_F2_9" = "Endothelial",
                                                "P3_F2_10" = "TSC-like EMT")
saveRDS(patient.list, "data/patient_seurat_split_id_list_final.rds")

p1.cols <- c("Early hepatoblast" = "#EDB074",
             "Late hepatoblast" = "#E37A2B",
             "Infiltrated" = "#69B564",
             "TAMs" = "#77D195",
             "Stroma" = "#B5B064",
             "Endothelial" = "#92B564")
p2.cols <- c("Mature hepatocyte" = "#8A5146",
             "Hepatocyte 1" = "#ED7474",
             "Hepatocyte 2" = "#DE3535",
             "Hepatocyte 3" = "#AD2B2B",
             "Hepatocytic" = "#ED74B0",
             "Late hepatoblast" = "#E37A2B",
             "Infiltrated" = "#69B564",
             "Stroma" = "#B5B064",
             "Endothelial" = "#92B564")
p3.cols <- c("Hepatocyte" = "#ED7474",
             "Astrocyte 1" = "#AB9AED",
             "Astrocyte 2" = "#D49AED",
             "HSC-like" = "#82CFBA",
             "TSC-like EMT" = "#8297CF",
             "TSC-like TGFb" = "#82BDCF",
             "MSC-like" = "#2E5F75",
             "Infiltrated" = "#69B564",
             "Stroma" = "#B5B064",
             "Endothelial" = "#92B564")
final.cols <- c(p1.cols, p2.cols, p3.cols)
final.cols <- final.cols[!duplicated(names(final.cols))]
saveRDS(final.cols, "data/final_cluster_cols.rds")

all.clusters <-c("Early hepatoblast", "Late hepatoblast",
                 "Hepatocytic", "Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3",
                 "Hepatocyte", "Mature hepatocyte",
                 "HSC-like", "MSC-like", "TSC-like TGFb", "TSC-like EMT",
                 "Astrocyte 1", "Astrocyte 2",
                 "Infiltrated", "TAMs", "Stroma", "Endothelial")
for (patient in names(patient.list)) {
  seurat <- patient.list[[patient]]
  order <- all.clusters[all.clusters %in% unique(seurat$Final.Cluster)]
  seurat$Final.Cluster <- factor(seurat$Final.Cluster, levels = order)
  DimPlot(seurat, group.by = "Final.Cluster", order = TRUE) +
    scale_colour_manual(values = final.cols) +
    umap.theme() + labs(title = paste(patient, " - Final Cluster"))
  ggsave(paste0("plots/umaps/umap_final_", patient, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_final_", patient, ".png"), width = 8.3, height = 5.8)
}

meta.combined <- do.call(rbind, lapply(seq_along(patient.list), function(i) {
  df <- patient.list[[i]]@meta.data[, c("Condition", "SingleR.labels", "Final.Cluster")]
  df$Patient <- names(patient.list)[i]
  return(df)}))
type.df <- meta.combined %>%
  group_by(Condition, SingleR.labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
anno.cols <- readRDS("data/singler_cols.rds")
ggplot(type.df, aes(x = Condition, y = prop, fill = SingleR.labels)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = anno.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("plots/anno/patient_singler_bar.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/patient_singler_bar.png", width = 8.3, height = 5.8)

cluster.df <- meta.combined %>%
  group_by(Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
ggplot(cluster.df, aes(x = Condition, y = prop, fill = Final.Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = final.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("plots/anno/patient_final_cluster_bar.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/patient_final_cluster_bar.png", width = 8.3, height = 5.8)

cluster.filt <- meta.combined %>%
  filter(!Final.Cluster %in% c("TAMs", "Infiltrated", "Endothelial", "Stroma")) %>%
  group_by(Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
cluster.filt$Patient <- sub("_[^_]+$", "", cluster.filt$Condition)
patients <- unique(cluster.filt$Patient)
for (patient in patients) {
  df <- cluster.filt[cluster.filt$Patient == patient, ]
  ggplot(df, aes(x = Condition, y = prop, fill = Final.Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = final.cols) +
    scale_y_continuous(labels = scales::percent_format()) +
    ylab("Proportion") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          aspect.ratio = 1.5)
  ggsave(paste0("plots/anno/", patient,"_final_cluster_bar.pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/anno/", patient,"_final_cluster_bar.png"), width = 5.8, height = 5.8)
}

## /////////////////////////////////////////////////////////////////////////////
## InferCNV scores /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

patient.seurat <- readRDS("data/patient_seurat_infercnv.rds")
patient.list <- readRDS("data/patient_seurat_split_id_list_final.rds")

cnv.scores <- grep("^(infercnv|chr).*_score$", colnames(patient.seurat@meta.data), value = TRUE)
for (i in seq_along(patient.list)) {
  obj <- patient.list[[i]]
  common.cells <- intersect(colnames(obj), colnames(patient.seurat))
  meta.df <- patient.seurat@meta.data[common.cells, cnv.scores, drop = FALSE]
  obj@meta.data[common.cells, cnv.scores] <- meta.df
  patient.list[[i]] <- obj
}
saveRDS(patient.list, "data/patient_seurat_split_id_list_infercnv_final.rds")

make_title <- function(x) {
  x |> 
    gsub("_", " ", x = _) |> 
    tools::toTitleCase()
}
for (patient in names(patient.list)) {
  seurat <- patient.list[[patient]]
  for (score in cnv.scores) {
    FeaturePlot(seurat, features = score) +
      scale_colour_gradientn(colours = pals::coolwarm(100),
                             values = scales::rescale(c(min(seurat[[score]]), 1, max(seurat[[score]])))) +
      labs(title = paste(make_title(score)), colour = "") +
      umap.theme() +
      theme(text = element_text(size = 15),
            legend.key.size = unit(1, "cm"))
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_umap.pdf"), width = 8.3, height = 5.8)
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_umap.png"), width = 8.3, height = 5.8)
  }
}

init.cols <- readRDS("data/initial_cluster_cols.rds")
for (patient in names(patient.list)) {
  seurat <- patient.list[[patient]]
  for (score in cnv.scores) {
    Idents(seurat) <- seurat$Init.Cluster
    VlnPlot(seurat, features = score,
            pt.size = 0, cols = init.cols) &
      xlab("") &
      labs(title = paste(make_title(score))) &
      theme(legend.position = "none",
            text = element_text(size = 12),
            aspect.ratio = 0.5)
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_initial_violin.pdf"), width = 8.3, height = 5.8)
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_initial_violin.png"), width = 8.3, height = 5.8)
  }
}

final.cols <- readRDS("data/final_cluster_cols.rds")
all.clusters <-c("Early hepatoblast", "Late hepatoblast",
                 "Hepatocytic", "Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3",
                 "Hepatocyte", "Mature hepatocyte",
                 "HSC-like", "MSC-like", "TSC-like TGFb", "TSC-like EMT",
                 "Astrocyte 1", "Astrocyte 2",
                 "Infiltrated", "TAMs", "Stroma", "Endothelial")
for (patient in names(patient.list)) {
  seurat <- patient.list[[patient]]
  for (score in cnv.scores) {
    order <- all.clusters[all.clusters %in% unique(seurat$Final.Cluster)]
    seurat$Final.Cluster <- factor(seurat$Final.Cluster, levels = order)
    Idents(seurat) <- seurat$Final.Cluster
    VlnPlot(seurat, features = score,
            pt.size = 0, cols = final.cols) &
      xlab("") &
      labs(title = paste(make_title(score))) &
      theme(legend.position = "none",
            text = element_text(size = 12),
            aspect.ratio = 0.5)
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_final_violin.pdf"), width = 8.3, height = 5.8)
    ggsave(paste0("plots/anno/infercnv/", patient, "_", score, "_final_violin.png"), width = 8.3, height = 5.8)
  }
}

## [ Patient markers ] ----

patient.list <- readRDS("data/patient_seurat_split_id_list_final.rds")

library(writexl)
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
for (patient in names(patient.list)) {
  obj <- patient.list[[patient]]
  
  Idents(obj) <- obj$Final.Cluster
  obj.markers <- FindAllMarkers(obj, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
  saveRDS(obj.markers, paste0("data/markers/", patient, "_markers.rds"))
  table(obj.markers$cluster)
  
  obj.list <- lapply(unique(obj.markers$cluster), function(x){
    tmp.df <- obj.markers[obj.markers$cluster==x & obj.markers$p_val_adj < 0.05, ]
    tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
  })
  names(obj.list) <- unique(obj.markers$cluster)
  saveRDS(obj.list, paste0("data/markers/", patient, "_markers_list.rds"))
  
  obj.GO.BP.list <- lapply(obj.list, function(x){
    enrichGO(x$gene,
             OrgDb = "org.Hs.eg.db",
             keyType = "SYMBOL",
             ont = "BP", 
             readable = TRUE, 
             pAdjustMethod = "BH")
  })
  
  C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
  C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
  obj.C2.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = C2.df)
  })
  
  H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
  H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
  obj.H.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = H.df)
  })
  saveRDS(obj.GO.BP.list, paste0("data/markers/", patient, "_markers_BP_ontology.rds"))
  saveRDS(obj.C2.list, paste0("data/markers/", patient, "_markers_C2_geneset.rds"))
  saveRDS(obj.H.list, paste0("data/markers/", patient, "_markers_Hallmarks.rds"))
  
  obj.lists <- list("obj.GO.BP.list" = obj.GO.BP.list,
                    "obj.C2.list" = obj.C2.list,
                    "obj.H.list" = obj.H.list,
                    "obj.list" = obj.list)
  obj.names <- sub("list", "sheets", names(obj.lists))
  for(i in 1:length(obj.lists)){
    obj.list <- obj.lists[[i]]
    list.sheets <- lapply(obj.list, as.data.frame)
    assign(obj.names[[i]], list.sheets)
  }
  write_xlsx(obj.GO.BP.sheets, paste0("data/markers/", patient, "_markers_BP_ontology.xlsx"))
  write_xlsx(obj.C2.sheets, paste0("data/markers/", patient, "_markers_C2_geneset.xlsx"))
  write_xlsx(obj.H.sheets, paste0("data/markers/", patient, "_markers_Hallmarks.xlsx"))
  write_xlsx(obj.sheets, paste0("data/markers/", patient, "_markers.xlsx"))
  
  gc()
}

## Differential expression between clusters of interest
exclude <-c("Infiltrated", "TAMs", "Stroma", "Endothelial")
for (patient in names(patient.list)) {
  obj <- patient.list[[patient]]
  include <- unique(obj$Final.Cluster)[!(unique(obj$Final.Cluster) %in% exclude)]
  obj <- subset(obj, subset = Final.Cluster %in% include)
  
  Idents(obj) <- obj$Final.Cluster
  obj.markers <- FindAllMarkers(obj, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
  saveRDS(obj.markers, paste0("data/markers_filt/", patient, "_markers.rds"))
  table(obj.markers$cluster)
  
  obj.list <- lapply(unique(obj.markers$cluster), function(x){
    tmp.df <- obj.markers[obj.markers$cluster==x & obj.markers$p_val_adj < 0.05, ]
    tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
  })
  names(obj.list) <- unique(obj.markers$cluster)
  saveRDS(obj.list, paste0("data/markers_filt/", patient, "_markers_list.rds"))
  
  obj.GO.BP.list <- lapply(obj.list, function(x){
    enrichGO(x$gene,
             OrgDb = "org.Hs.eg.db",
             keyType = "SYMBOL",
             ont = "BP", 
             readable = TRUE, 
             pAdjustMethod = "BH")
  })
  
  C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
  C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
  obj.C2.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = C2.df)
  })
  
  H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
  H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
  obj.H.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = H.df)
  })
  saveRDS(obj.GO.BP.list, paste0("data/markers_filt/", patient, "_markers_BP_ontology.rds"))
  saveRDS(obj.C2.list, paste0("data/markers_filt/", patient, "_markers_C2_geneset.rds"))
  saveRDS(obj.H.list, paste0("data/markers_filt/", patient, "_markers_Hallmarks.rds"))
  
  obj.lists <- list("obj.GO.BP.list" = obj.GO.BP.list,
                    "obj.C2.list" = obj.C2.list,
                    "obj.H.list" = obj.H.list,
                    "obj.list" = obj.list)
  obj.names <- sub("list", "sheets", names(obj.lists))
  for(i in 1:length(obj.lists)){
    obj.list <- obj.lists[[i]]
    list.sheets <- lapply(obj.list, as.data.frame)
    assign(obj.names[[i]], list.sheets)
  }
  write_xlsx(obj.GO.BP.sheets, paste0("data/markers_filt/", patient, "_markers_BP_ontology.xlsx"))
  write_xlsx(obj.C2.sheets, paste0("data/markers_filt/", patient, "_markers_C2_geneset.xlsx"))
  write_xlsx(obj.H.sheets, paste0("data/markers_filt/", patient, "_markers_Hallmarks.xlsx"))
  write_xlsx(obj.sheets, paste0("data/markers_filt/", patient, "_markers.xlsx"))
  
  gc()
}

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk //////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

patient.list <- readRDS("data/patient_seurat_split_id_list_final.rds")
pseudo.list <- list()
for (patient in names(patient.list)) {
  obj <- patient.list[[patient]]
  pseudo <- AggregateExpression(obj, assays = "RNA", return.seurat = T,
                                group.by = c("Final.Cluster"))
  pseudo.list[[patient]] <- pseudo
}
saveRDS(pseudo.list, "data/patient_pseudobulk_split_list_final.rds")

final.cols <- readRDS("data/final_cluster_cols.rds")
all.clusters <-c("Early hepatoblast", "Late hepatoblast",
                 "Hepatocytic", "Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3",
                 "Hepatocyte", "Mature hepatocyte",
                 "HSC-like", "MSC-like", "TSC-like TGFb", "TSC-like EMT",
                 "Astrocyte 1", "Astrocyte 2",
                 "Infiltrated", "TAMs", "Stroma", "Endothelial")

## HB signatures from literature
sig.df <- as.data.frame(readxl::read_xlsx(path = "data/hb_sigs.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")
sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})
names(sigs) <- unique(sig.df$Signature)
sigs
sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)

for (patient in names(pseudo.list)) {
  pseudo <- pseudo.list[[patient]]
  pseudo <- AddModuleScore(pseudo, features = sigs, assay = "RNA", seed = 12345, 
                                name = sig.names)
  colnames(pseudo@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo@meta.data))
  pseudo[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo, vars = sig.names)))
  pseudo.list[[patient]] <- pseudo
  hb.mat <- as.matrix(pseudo@assays$HB_sigs@data)
  order <- all.clusters[all.clusters %in% unique(pseudo$Final.Cluster)]
  pseudo$Final.Cluster <- factor(pseudo$Final.Cluster, levels = order)
  anno.df <- pseudo@meta.data["Final.Cluster"]
  cols.sub <- list(Final.Cluster = final.cols[order])
  gt <- pheatmap(hb.mat,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 annotation_col = anno.df,
                 annotation_colors = cols.sub,
                 main = patient)$gtable
  ggsave(paste0("plots/anno/pseudo_final_hb_sigs_pheatmap_", patient, ".pdf"), plot = gt)
  ggsave(paste0("plots/anno/pseudo_final_hb_sigs_pheatmap_", patient, ".png"), plot = gt)
}

exclude <- c("TAMs", "Infiltrated", "Endothelial", "Stroma")
pseudo.filt <- lapply(pseudo.list, function(pseudo) {
  keep.cells <- setdiff(colnames(pseudo), exclude)
  subset(pseudo, cells = keep.cells)
})

for (patient in names(pseudo.filt)) {
  pseudo <- pseudo.filt[[patient]]
  hb.mat <- as.matrix(pseudo@assays$HB_sigs@data)
  order <- all.clusters[all.clusters %in% unique(pseudo$Final.Cluster)]
  pseudo$Final.Cluster <- factor(pseudo$Final.Cluster, levels = order)
  anno.df <- pseudo@meta.data["Final.Cluster"]
  cols.sub <- list(Final.Cluster = final.cols[order])
  gt <- pheatmap(hb.mat,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 annotation_col = anno.df,
                 annotation_colors = cols.sub,
                 main = patient)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_final_hb_sigs_pheatmap_", patient, "_tumour.pdf"), plot = gt)
  ggsave(filename = paste0("plots/anno/pseudo_final_hb_sigs_pheatmap_", patient, "_tumour.png"), plot = gt)
}

## Liver differentiation markers
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

for (patient in names(pseudo.list)) {
  pseudo <- pseudo.list[[patient]]
  custom.mat <- as.matrix(pseudo@assays$RNA$scale.data[custom, ])
  order <- all.clusters[all.clusters %in% unique(pseudo$Final.Cluster)]
  pseudo$Final.Cluster <- factor(pseudo$Final.Cluster, levels = order)
  anno.df <- pseudo@meta.data["Final.Cluster"]
  cols.sub <- list(Final.Cluster = final.cols[order])
  gt <- pheatmap(custom.mat,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 annotation_col = anno.df,
                 annotation_colors = cols.sub,
                 main = patient)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", patient, ".pdf"), plot = gt)
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", patient, ".png"), plot = gt)
}

# exclude <- c("TAMs", "Infiltrated", "Endothelial", "Stroma")
pseudo.filt <- lapply(pseudo.list, function(pseudo) {
  keep.cells <- setdiff(colnames(pseudo), exclude)
  subset(pseudo, cells = keep.cells)
})

for (patient in names(pseudo.filt)) {
  pseudo <- pseudo.filt[[patient]]
  custom.mat <- as.matrix(pseudo@assays$RNA$scale.data[custom, ])
  order <- all.clusters[all.clusters %in% unique(pseudo$Final.Cluster)]
  pseudo$Final.Cluster <- factor(pseudo$Final.Cluster, levels = order)
  anno.df <- pseudo@meta.data["Final.Cluster"]
  cols.sub <- list(Final.Cluster = final.cols[order])
  gt <- pheatmap(custom.mat,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 annotation_col = anno.df,
                 annotation_colors = cols.sub,
                 main = patient)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", patient, "_tumour.pdf"), plot = gt)
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", patient, "_tumour.png"), plot = gt)
}

## [ Relapse markers ] ----

patient.list <- readRDS("data/patient_seurat_split_id_list_final.rds")

## Filter for tumour cells
patient.filt <- lapply(
  patient.list,
  function(obj) {
    cells.keep <- rownames(obj@meta.data)[
      !obj$Final.Cluster %in% c("Infiltrated", "TAMs", "Stroma", "Endothelial")
    ]
    subset(obj, cells = cells.keep)
  }
)

library(writexl)
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
for (patient in names(patient.filt)) {
  obj <- patient.filt[[patient]]
  
  Idents(obj) <- obj$Fraction
  obj.markers <- FindAllMarkers(obj, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
  saveRDS(obj.markers, paste0("data/markers_relapse/", patient, "_markers.rds"))
  table(obj.markers$cluster)
  
  obj.list <- lapply(unique(obj.markers$cluster), function(x){
    tmp.df <- obj.markers[obj.markers$cluster==x & obj.markers$p_val_adj < 0.05, ]
    tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
  })
  names(obj.list) <- unique(obj.markers$cluster)
  saveRDS(obj.list, paste0("data/markers_relapse/", patient, "_markers_list.rds"))
  
  obj.GO.BP.list <- lapply(obj.list, function(x){
    enrichGO(x$gene,
             OrgDb = "org.Hs.eg.db",
             keyType = "SYMBOL",
             ont = "BP", 
             readable = TRUE, 
             pAdjustMethod = "BH")
  })
  
  C2.msigsdb <- msigdbr(species = "Homo sapiens", collection = "C2")
  C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])
  obj.C2.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = C2.df)
  })
  
  H.msigsdb <- msigdbr(species = "Homo sapiens", collection = "H")
  H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])
  obj.H.list <- lapply(obj.list, function(x){
    enricher(x$gene,
             pAdjustMethod = "BH", TERM2GENE = H.df)
  })
  saveRDS(obj.GO.BP.list, paste0("data/markers_relapse/", patient, "_markers_BP_ontology.rds"))
  saveRDS(obj.C2.list, paste0("data/markers_relapse/", patient, "_markers_C2_geneset.rds"))
  saveRDS(obj.H.list, paste0("data/markers_relapse/", patient, "_markers_Hallmarks.rds"))
  
  obj.lists <- list("obj.GO.BP.list" = obj.GO.BP.list,
                    "obj.C2.list" = obj.C2.list,
                    "obj.H.list" = obj.H.list,
                    "obj.list" = obj.list)
  obj.names <- sub("list", "sheets", names(obj.lists))
  for(i in 1:length(obj.lists)){
    obj.list <- obj.lists[[i]]
    list.sheets <- lapply(obj.list, as.data.frame)
    assign(obj.names[[i]], list.sheets)
  }
  write_xlsx(obj.GO.BP.sheets, paste0("data/markers_relapse/", patient, "_markers_BP_ontology.xlsx"))
  write_xlsx(obj.C2.sheets, paste0("data/markers_relapse/", patient, "_markers_C2_geneset.xlsx"))
  write_xlsx(obj.H.sheets, paste0("data/markers_relapse/", patient, "_markers_Hallmarks.xlsx"))
  write_xlsx(obj.sheets, paste0("data/markers_relapse/", patient, "_markers.xlsx"))
  
  gc()
}

## [ Patient 3 markers ] ----

## Check patient 3 for expression of canonical neuroendocrine markers
## e.g. synaptophysin (SYP), chromogranin (CHGA), serotonin (SLC6A4), and CD56 (NCAM)

pseudo.list <- readRDS("data/patient_pseudobulk_split_list_final.rds")
patient3.pseudo <- pseudo.list$`16-00875`

custom <- c("SYP", "CHGA", "SLC6A4", "NCAM1",
            "ASCL1", "NEUROD1", "POU2F3",
            "POU5F1", "SOX2", "NANOG",
            "CD34", "PROM1")
custom %in% rownames(patient3.pseudo)
custom.mat <- as.matrix(patient3.pseudo@assays$RNA$scale.data[custom, ])
all.clusters <-c("Early hepatoblast", "Late hepatoblast",
                 "Hepatocytic", "Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3",
                 "Hepatocyte", "Mature hepatocyte",
                 "HSC-like", "MSC-like", "TSC-like TGFb", "TSC-like EMT",
                 "Astrocyte 1", "Astrocyte 2",
                 "Infiltrated", "TAMs", "Stroma", "Endothelial")
order <- all.clusters[all.clusters %in% unique(patient3.pseudo$Final.Cluster)]
patient3.pseudo$Final.Cluster <- factor(patient3.pseudo$Final.Cluster, levels = order)
anno.df <- patient3.pseudo@meta.data["Final.Cluster"]
final.cols <- readRDS("data/final_cluster_cols.rds")
cols.sub <- list(Final.Cluster = final.cols[order])
gt <- pheatmap(custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = cols.sub,
               main = "Patient 3")$gtable
ggsave("plots/anno/pseudo_final_neuro_pheatmap_16-00875.pdf", plot = gt)
ggsave("plots/anno/pseudo_final_neuro_pheatmap_16-00875.png", plot = gt)

## Malignant cells only
exclude <- c("TAMs", "Infiltrated", "Endothelial", "Stroma")
keep.cells <- setdiff(colnames(patient3.pseudo), exclude)
pseudo.filt <- subset(patient3.pseudo, cells = keep.cells)

custom.mat <- as.matrix(pseudo.filt@assays$RNA$scale.data[custom, ])
order <- all.clusters[all.clusters %in% unique(pseudo.filt$Final.Cluster)]
pseudo.filt$Final.Cluster <- factor(pseudo.filt$Final.Cluster, levels = order)
anno.df <- pseudo.filt@meta.data["Final.Cluster"]
cols.sub <- list(Final.Cluster = final.cols[order])
gt <- pheatmap(custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = cols.sub,
               main = "Patient 3")$gtable
ggsave("plots/anno/pseudo_final_neuro_pheatmap_16-00875_tumour.pdf", plot = gt)
ggsave("plots/anno/pseudo_final_neuro_pheatmap_16-00875_tumour.png", plot = gt)

## [ Supplementary Figure 3 ] ----

## Plot Supplementary Fig. 3 without P3
pseudo.init <- readRDS("data/patient_pseudobulk_init_sigs.rds")
pseudo.init <- subset(pseudo.init, cells = colnames(pseudo.init)[!grepl("^P3", colnames(pseudo.init))])
colnames(pseudo.init)

## Liver differentiation markers
## (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/; https://pmc.ncbi.nlm.nih.gov/articles/PMC11114060/)
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
custom %in% rownames(pseudo.init)

init.custom.mat <- as.matrix(pseudo.init@assays$RNA$scale.data[custom, ])
init.anno.df <- pseudo.init@meta.data["Init.Cluster"]
init.cols <- readRDS("data/initial_cluster_cols.rds")
init.cols <- init.cols[!grepl("^P3", names(init.cols))]
names(init.cols) <- gsub("_", "-", names(init.cols))
gt <- pheatmap(init.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/suppfig3/pseudo_init_custom_pheatmap_p1_p2.pdf", plot = gt, width = 11.7, height = 11.7)
ggsave("plots/suppfig3/pseudo_init_custom_pheatmap_p1_p2.png", plot = gt, width = 11.7, height = 11.7)

init.hb.mat <- as.matrix(pseudo.init@assays$HB_sigs@data)
gt <- pheatmap(init.hb.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = init.anno.df,
               annotation_colors = list(Init.Cluster = init.cols))$gtable
ggsave("plots/suppfig3/pseudo_init_hb_sigs_pheatmap_p1_p2.pdf", plot = gt, width = 8.3, height = 11.7)
ggsave("plots/suppfig3/pseudo_init_hb_sigs_pheatmap_p1_p2.png", plot = gt, width = 8.3, height = 11.7)

patient.list <- readRDS("data/patient_seurat_split_id_list_final.rds")
patient.list$`16-00875` <- NULL
meta.combined <- do.call(rbind, lapply(seq_along(patient.list), function(i) {
  df <- patient.list[[i]]@meta.data[, c("Condition", "SingleR.labels", "Final.Cluster")]
  df$Patient <- names(patient.list)[i]
  return(df)}))
type.df <- meta.combined %>%
  group_by(Condition, SingleR.labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
anno.cols <- readRDS("data/singler_cols.rds")
ggplot(type.df, aes(x = Condition, y = prop, fill = SingleR.labels)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = anno.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("plots/suppfig3/patient_singler_bar_p1_p2.pdf", width = 8.3, height = 5.8)
ggsave("plots/suppfig3/patient_singler_bar_p1_p2.png", width = 8.3, height = 5.8)

cluster.df <- meta.combined %>%
  group_by(Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
final.cols <- readRDS("data/final_cluster_cols.rds")
ggplot(cluster.df, aes(x = Condition, y = prop, fill = Final.Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = final.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("plots/suppfig3/patient_final_cluster_bar_p1_p2.pdf", width = 8.3, height = 5.8)
ggsave("plots/suppfig3/patient_final_cluster_bar_p1_p2.png", width = 8.3, height = 5.8)
