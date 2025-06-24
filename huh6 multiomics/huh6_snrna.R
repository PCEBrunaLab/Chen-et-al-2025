## HuH6 snRNA-seq data analysis

## /////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(scran)
library(scater)
library(scuttle)
library(Seurat)
library(SeuratObject)
library(DropletUtils)
library(MatrixGenerics)
library(igraph)
library(BiocParallel)
library(biomaRt)
library(harmony)
library(cluster)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(celldex)
library(SingleR)
library(tidyverse)
library(ggplot2)
library(ggsankey)
library(patchwork)
library(pheatmap)
library(ggbeeswarm)
library(cowplot)
library(reshape2)
library(hues)
library(readxl)
library(writexl)

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

## Sample colours
huh6.cols <- c("Recovered_A" = "#117733",
               "Recovered_B" = "#44AA99",
               "Untreated_A" = "#88CCEE",
               "Untreated_B" = "#332288")

## QC palette
cool.warm.pal <- colorRampPalette(c("#1d4877", "#1b8a5a", "#fbb021", "#f68838", "#ee3e32"))

## Colours for PAM clusters
require(paletteer)
pal <- paletteer_d("rcartocolor::Pastel")
pam.cols <- c("Hepatocytic 1" = pal[2], 
              "Hepatocytic 2" = pal[4], 
              "Hepatocytic 3" = pal[3], 
              "Late progenitor" = pal[5],
              "Early progenitor" = pal[1])

## [ Data prep ] ----

## Load in data
hb.dirs <- list.dirs("/data/scratch/DMP/DUDMP/PAEDCANC/shamer/HB/multiomics", recursive = FALSE)
to.match <- c("fastq", "logs", "src")
hb.dirs <- hb.dirs[lapply(hb.dirs, function(x) length(grep(paste(to.match, collapse = "|"), x, value = FALSE))) == 0]
hb.data <- list()

for (sample in 1:length(hb.dirs)) {
  hb.data[sample] <- read10xCounts(paste0(hb.dirs[sample],"/outs/raw_feature_bc_matrix/")) 
}

hb.names <- tolower(list.dirs("/data/scratch/DMP/DUDMP/PAEDCANC/shamer/HB/multiomics", recursive = FALSE, full.names = FALSE))
hb.names <- hb.names[lapply(hb.names, function(x) length(grep(paste(to.match, collapse = "|"), x, value = FALSE))) == 0]
names(hb.data) <- hb.names

## Remove empty drops
for (sample in 1:length(hb.data)) {
  colData(hb.data[[sample]])$Sample <- gsub(pattern = ".*/multiomics/", replacement = "", colData(hb.data[[sample]])$Sample)
  colData(hb.data[[sample]])$Sample <- gsub(pattern = "/outs.*", replacement = "", colData(hb.data[[sample]])$Sample)
}

dir.create("data/raw_data/", recursive = TRUE)

for (sample in 1:length(hb.data)) {
  data <- hb.data[[sample]]
  ## Calculate empty drops for each sample
  intest_calls <- emptyDrops(data, assay.type = "counts", niters = 20000, 
                             ignore = 4999, lower = 500, retain = Inf)
  ## Identify cells using a false discovery rate (FDR) of 1%
  sig_cells <- intest_calls$FDR <= 0.01 & !is.na(intest_calls$FDR)
  
  ## Subset the called cells from the SCE
  sce.drops <- data[, sig_cells]
  message("Keeping ", sum(sig_cells), " non-empty droplet barcodes")
  
  ## Save SCE as rds file for demux and cmo summary
  saveRDS(sce.drops, file = paste0("data/raw_data/", names(hb.data)[sample], ".rds"))
  gc()
}

## Attach metadata
## Load in data files
file.list <- list.files("data/raw_data", full.names = TRUE)
huh6.files <- file.list[c(3:6)]

huh6.list <- list()
for (file in 1:length(huh6.files)) {
  huh6.list[file] <- readRDS(huh6.files[file]) 
}

## Assign CellID column to metadata
for (sce in 1:length(huh6.list)) {
  colData(huh6.list[[sce]])$CellID <- paste0(colData(huh6.list[[sce]])$Sample, "_", colData(huh6.list[[sce]])$Barcode)
}

## Compute sizeFactor of each cell within samples
huh6.list <- lapply(huh6.list, computeLibraryFactors)
saveRDS(huh6.list, "data/huh6_sce_list.rds")

## Extract metadata from each sample
huh6.meta <- list()
for (sce in 1:length(huh6.list)) {
  huh6.meta[[sce]] <- as.data.frame(colData(huh6.list[[sce]]))
}
saveRDS(huh6.meta, "data/huh6_metadata_list.rds")

## [ Combine samples ] ----

huh6.list <- readRDS("data/huh6_sce_list.rds")

## Create merged SCE without batch correction
huh6.batches <- list()

for (sce in 1:length(huh6.list)) {
  huh6.batches[sce] <- unique(huh6.list[[sce]]$Sample)
}
names(huh6.list) <- huh6.batches

## Combine SCE objects without removing data
huh6.sce <- sce_cbind(huh6.list, exprs = c("counts"), colData_names = TRUE, batch_names = huh6.batches)
rm(huh6.list)

saveRDS(huh6.sce, "data/combined_huh6_sce.rds")

## [ QC filtering ] ----

## Load in data
huh6.sce <- readRDS("data/combined_huh6_sce.rds")

colnames(huh6.sce) <- huh6.sce$CellID
length(unique(colnames(huh6.sce))) ## 16454 cells

## Compute logcounts
huh6.sce <- logNormCounts(huh6.sce)

## ENSEMBL database
require(EnsDb.Hsapiens.v86)
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype == "lincRNA"]

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
#                 values = rownames(huh6.sce))
# 
# ## Save ensembl.csv for reuse
# write.csv(ens.bm, "data/huh6_ensembl_biomart.csv")
# gc()

ens.bm <- read.csv("data/huh6_ensembl_biomart.csv", row.names = NULL)
ens.bm$X <- NULL

## Identify duplicated gene names
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)] 
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])
## Duplicated gene names: EIF1B-AS1, NCMAP-DT

## Created filtered ens.bm and keep only useful information
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)

## This filters out around 50 genes
nrow(ens.filt.bm) ## 13596
nrow(huh6.sce) ## 13645

ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]

## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]

## Around 1000 genes filtered out
nrow(ens.filt.bm) ## 12476
nrow(huh6.sce) ## 13645

## Keep track of everything that didn't get annotated
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 1120

length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 1080 are novel transcripts or novel transcripts antisense to a gene

grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest are mostly either novel, mitochondrially encoded, or just empty description

## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## "EIF1B-AS1" "MATR3"     "NCMAP-DT" 

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 12473

## Identify mitochondrial, ribosomal and linc genes
mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Add QC metrics to SCE object
huh6.sce <- addPerCellQCMetrics(huh6.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))
huh6.sce <- huh6.sce[ens.filt.bm$ensembl_gene_id,]
rownames(huh6.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

## Assign rowData to SCE object with filtered ens.bm
rowData(huh6.sce) <- ens.filt.bm

saveRDS(huh6.sce, "data/huh6_sce.rds")

## Extract meta data from SCE object
huh6.meta <- as.data.frame(colData(huh6.sce))
huh6.meta$batch <- NULL

## Generate Seurat object from SCE object
huh6.seurat <- CreateSeuratObject(counts = assay(huh6.sce, "counts"),
                                  assay = "RNA",
                                  meta.data = huh6.meta)

saveRDS(huh6.seurat, "data/huh6_seurat.rds")

huh6.sce <- readRDS("data/huh6_sce.rds")

## Plot QC metrics
det.sce.gg <- plotColData(huh6.sce, y = "detected", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

umi.sce.gg <- plotColData(huh6.sce, y = "total", x = "Sample", colour_by = "Sample") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")

mt.sce.gg <- plotColData(huh6.sce, y = "subsets_mt_percent", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

rb.sce.gg <- plotColData(huh6.sce, y = "subsets_ribo_percent", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.5)

ggsave("plots/qc/huh6_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/huh6_qc_plots.png", width = 8.3, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
counts.assay <- counts(huh6.sce)
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

ggsave("plots/qc/huh6_qc_top_genes.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/huh6_qc_top_genes.png", width = 8.3, height = 5.8)

rm(counts.assay)
gc()

## Visualise MALAT1 expression and lincRNA per sample
malat1.sce.gg <- plotExpression(huh6.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Sample", 
                                colour_by = "Sample") +
  labs(x = element_blank(), y = "Log2 Expression") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

lincrna.sce.gg <- plotColData(huh6.sce, y = "subsets_linc_percent", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.7)

ggsave("plots/qc/huh6_qc_malat1_lincrna.pdf", width = 5.8, height = 8.3)
ggsave("plots/qc/huh6_qc_malat1_lincrna.png", width = 5.8, height = 8.3)

## Cell cycle scoring
huh6.seurat <- readRDS("data/huh6_seurat.rds")

huh6.seurat <- NormalizeData(huh6.seurat)
huh6.seurat <- CellCycleScoring(huh6.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

huh6.sce$Seurat.Phase <- huh6.seurat$Phase
huh6.sce$Seurat.S <- huh6.seurat$S.Score
huh6.sce$Seurat.G2M <- huh6.seurat$G2M.Score

saveRDS(huh6.sce, "data/huh6_sce_cell_cycle.rds")
saveRDS(huh6.seurat, "data/huh6_seurat_cell_cycle.rds")

huh6.sce <- readRDS("data/huh6_sce_cell_cycle.rds")
huh6.seurat <- readRDS("data/huh6_seurat_cell_cycle.rds")

## Check gene and cell meta data
colnames(colData(huh6.sce))
colnames(rowData(huh6.sce))

table(huh6.sce$Sample)
## Recovered_A Recovered_B Untreated_A Untreated_B 
## 573        4388        5211        6282 

det.sce.gg <- plotColData(huh6.sce, y = "detected", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_hline(yintercept = 800, linetype = "dotted")

mt.sce.gg <- plotColData(huh6.sce, y = "subsets_mt_percent", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dotted")

rb.sce.gg <- plotColData(huh6.sce, y = "subsets_ribo_percent", x = "Sample", colour_by = "Sample") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = huh6.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed")

det.sce.gg + mt.sce.gg + rb.sce.gg & theme(aspect.ratio = 0.5)
ggsave("plots/qc/huh6_qc_plots_filtering.pdf", width = 11.7, height = 5.8)
ggsave("plots/qc/huh6_qc_plots_filtering.png", width = 11.7, height = 5.8)

## Set the expression threshold over 500
detected.filt <- colnames(huh6.sce)[huh6.sce$detected > 800] ## Removes ~3000 cells
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(huh6.sce)[Matrix::rowSums(counts(huh6.sce)) > 5] ## Doesn't remove any genes
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(huh6.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

huh6.filt.sce <- huh6.sce[genes.filt, detected.filt]

dim(huh6.sce) - dim(huh6.filt.sce)
## We lose 2740 cells

mito.filt <- colnames(huh6.filt.sce)[huh6.filt.sce$subsets_mt_percent < 30]
#ribo.filt <- colnames(huh6.filt.sce)[huh6.filt.sce$subsets_ribo_percent > 5]
#qc.filt <- intersect(mito.filt, ribo.filt)

huh6.filt.sce <- huh6.filt.sce[, mito.filt]

dim(huh6.sce) - dim(huh6.filt.sce)
## We lose ~3000 cells in total

## Percentage of the library that we filter out
round(100 - (100 * (table(huh6.filt.sce$Sample) / table(huh6.sce$Sample))), 1)
## Recovered_A Recovered_B Untreated_A Untreated_B 
## 13.6        21.2        16.6        20.2 

table(huh6.filt.sce$Sample)
## Recovered_A Recovered_B Untreated_A Untreated_B 
## 495        3456        4345        5012  

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(huh6.seurat)[Matrix::rowSums(GetAssayData(huh6.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(huh6.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

huh6.filt.seurat <- subset(huh6.seurat, 
                           cells = WhichCells(huh6.seurat, 
                                              expression = subsets_mt_percent < 30 &
                                                detected > 800),
                           features = genes.filt)

## The dimensions are equal
dim(huh6.filt.seurat) == dim(huh6.filt.sce)

saveRDS(huh6.filt.sce, "data/huh6_sce_filtered.rds")
saveRDS(huh6.filt.seurat, "data/huh6_seurat_filtered.rds")

## [ Filter HVGs ] ----

huh6.filt.seurat <- readRDS("data/huh6_seurat_filtered.rds")

## Filter HVGs
huh6.filt.seurat <- NormalizeData(huh6.filt.seurat)
huh6.filt.seurat <- FindVariableFeatures(huh6.filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(huh6.filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(huh6.filt.seurat), 20), repel = TRUE)
ggsave("plots/filter_hvgs/huh6_filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/huh6_filter_hvgs.png", width = 8.3, height = 5.8)

## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(huh6.filt.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(huh6.filt.seurat), value = TRUE),
              grep("[\\.]", rownames(huh6.filt.seurat), value = TRUE),
              grep("^LINC", rownames(huh6.filt.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
VariableFeatures(huh6.filt.seurat) <- VariableFeatures(huh6.filt.seurat)[which(!VariableFeatures(huh6.filt.seurat) %in% rm.genes)]
## Scale data and score cell cycle
huh6.filt.seurat <- ScaleData(huh6.filt.seurat, features = VariableFeatures(huh6.filt.seurat))

tmp.seurat <- CellCycleScoring(object = huh6.filt.seurat, 
                               g2m.features = cc.genes$g2m.genes, 
                               s.features = cc.genes$s.genes)
tmp.seurat$Cycle.Score <- tmp.seurat$S.Score - tmp.seurat$G2M.Score

seurat.cycle.melt <- make_long(
  data.frame("Sample" = huh6.filt.seurat$Sample, 
             "SCTransform" = tmp.seurat$Phase,
             "Normal" = huh6.filt.seurat$Phase),
  SCTransform, Normal)

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

DefaultAssay(huh6.filt.seurat) <- "RNA"
huh6.filt.seurat$Seurat.Cycle.Score <- huh6.filt.seurat$S.Score - huh6.filt.seurat$G2M.Score

saveRDS(huh6.filt.seurat, "data/huh6_seurat_hvgs.rds")

huh6.seurat <- readRDS("data/huh6_seurat_hvgs.rds")

huh6.seurat <- RunPCA(huh6.seurat, verbose = FALSE, npcs = 100, 
                      features = VariableFeatures(huh6.seurat))
ElbowPlot(object = huh6.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- RunHarmony(huh6.seurat, group.by.vars = "Sample",
                             theta = 0.5, lambda = NULL, sigma = 0.05,
                             assay.use = "RNA", reduction = "pca",
                             dims.use = 1:30, reduction.save = "Harmony",
                             max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Sample", order = TRUE) +
  scale_colour_manual(values = huh6.cols) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("plots/filter_hvgs/umap_harmony_t05_lnull_s005.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/umap_harmony_t05_lnull_s005.png", width = 8.3, height = 5.8)

# huh6.seurat <- RunUMAP(huh6.seurat, reduction = "pca", dims = 1:30)
# huh6.seurat <- FindNeighbors(huh6.seurat, reduction = "pca", dims = 1:30)

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
saveRDS(harmony.seurat, "data/huh6_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
huh6.seurat <- readRDS("data/huh6_seurat_hvgs.rds")
obj.list <- SplitObject(huh6.seurat, split.by = "Sample")

## Rerun variable feature selection
obj.list <- lapply(obj.list, FindVariableFeatures)
## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(huh6.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(huh6.seurat), value = TRUE),
              grep("[\\.]", rownames(huh6.seurat), value = TRUE),
              grep("^LINC", rownames(huh6.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(huh6.seurat, obj.list, obj)

## Scale data
output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
## Run PCA
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 100, 
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
saveRDS(output, "data/huh6_seurat_split_list.rds")

huh6.list <- readRDS("data/huh6_seurat_split_list.rds")
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_harmony.rds")
## Look through UMAPs and choose best cluster resolution
Idents(huh6.list$Untreated_A) <- huh6.list$Untreated_A$seurat_clusters.0.4
Idents(huh6.list$Untreated_B) <- huh6.list$Untreated_B$seurat_clusters.0.4
Idents(huh6.list$Recovered_A) <- huh6.list$Recovered_A$seurat_clusters.0.4
Idents(huh6.list$Recovered_B) <- huh6.list$Recovered_B$seurat_clusters.0.4
Idents(huh6.seurat) <- huh6.seurat$seurat_clusters.0.4

## Get cluster markers - MAST takes longer to run
uta.markers <- FindAllMarkers(huh6.list$Untreated_A, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
utb.markers <- FindAllMarkers(huh6.list$Untreated_B, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
reca.markers <- FindAllMarkers(huh6.list$Recovered_A, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
recb.markers <- FindAllMarkers(huh6.list$Recovered_B, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
huh6.all.markers <- FindAllMarkers(huh6.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

dir.create("data/jaccard/", recursive = TRUE)
saveRDS(uta.markers, "data/jaccard/huh6_ut_a_markers.rds")
saveRDS(utb.markers, "data/jaccard/huh6_ut_b_markers.rds")
saveRDS(reca.markers, "data/jaccard/huh6_rec_a_markers.rds")
saveRDS(recb.markers, "data/jaccard/huh6_rec_b_markers.rds")
saveRDS(huh6.all.markers, "data/jaccard/huh6_all_markers.rds")

## Filter significant and higher expressed markers
uta.markers <- uta.markers[uta.markers$p_val_adj < 0.05, ]
utb.markers <- utb.markers[utb.markers$p_val_adj < 0.05, ]
reca.markers <- reca.markers[reca.markers$p_val_adj < 0.05, ]
recb.markers <- recb.markers[recb.markers$p_val_adj < 0.05, ]
huh6.all.markers <- huh6.all.markers[huh6.all.markers$p_val_adj < 0.05, ]

uta.markers <- uta.markers[uta.markers$avg_log2FC > 0.5, ]
utb.markers <- utb.markers[utb.markers$avg_log2FC > 0.5, ]
reca.markers <- reca.markers[reca.markers$avg_log2FC > 0.5, ]
recb.markers <- recb.markers[recb.markers$avg_log2FC > 0.5, ]
huh6.all.markers <- huh6.all.markers[huh6.all.markers$avg_log2FC > 0.5, ]
## Split into lists of genes per cluster
uta.markers <- split(uta.markers, uta.markers$cluster)
utb.markers <- split(utb.markers, utb.markers$cluster)
reca.markers <- split(reca.markers, reca.markers$cluster)
recb.markers <- split(recb.markers, recb.markers$cluster)
huh6.all.markers <- split(huh6.all.markers, huh6.all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
lapply(uta.markers, nrow) ## 0 = 81, 1 = 368, 2 = 232, 3 = 624
lapply(utb.markers, nrow) ## 0 = 218, 1 = 435, 2 = 472, 3 = 350, 4 = 699
lapply(reca.markers, nrow) ## 0 = 365, 1 = 1187, 2 = 209
lapply(recb.markers, nrow) ## 0 = 289, 1 = 235, 2 = 999, 3 = 1043
lapply(huh6.all.markers, nrow) ## 0 = 288, 1 = 1117, 2 = 282, 3 = 361, 4 = 619, 5 = 979, 6 = 3073

uta.markers <- lapply(uta.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 80)
})
utb.markers <- lapply(utb.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
reca.markers <- lapply(reca.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
recb.markers <- lapply(recb.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
huh6.all.markers <- lapply(huh6.all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 250)
})

## Give your clusters unique names
names(uta.markers) <- paste0("UT_A_", names(uta.markers))
names(utb.markers) <- paste0("UT_B_", names(utb.markers))
names(reca.markers) <- paste0("Rec_A_", names(reca.markers))
names(recb.markers) <- paste0("Rec_B_", names(recb.markers))
names(huh6.all.markers) <- paste0("HuH6_all_", names(huh6.all.markers))
## Make a big list of all the cluster markers
huh6.markers.list <- c(uta.markers, utb.markers,
                       reca.markers, recb.markers,
                       huh6.all.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
huh6.jaccard.mat <- matrix(data = NA, nrow = length(huh6.markers.list),
                           ncol = length(huh6.markers.list),
                           dimnames = list(names(huh6.markers.list),
                                           names(huh6.markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(huh6.jaccard.mat)) {
  for (j in colnames(huh6.jaccard.mat)) {
    huh6.jaccard.mat[i,j] <- jaccard(huh6.markers.list[[i]]$gene, huh6.markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Timepoint" = gsub("^(.*?_[^_]*)_.*$", "\\1", colnames(huh6.jaccard.mat)),
                      row.names = colnames(huh6.jaccard.mat))
anno.col <- list("Timepoint" = c("HuH6_all" = "black",
                                 "UT_A" = as.character(huh6.cols[3]),
                                 "UT_B" = as.character(huh6.cols[4]),
                                 "Rec_A" = as.character(huh6.cols[1]),
                                 "Rec_B" = as.character(huh6.cols[2])))
dev.off()
gt <- pheatmap(huh6.jaccard.mat,
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
ggsave("plots/jaccard/jaccard_heatmap.png", plot = gt)
ggsave("plots/jaccard/jaccard_heatmap.pdf", plot = gt)

## Run the parameter search for PAM
huh6.jaccard.dist <- as.dist(1 - huh6.jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
for (k in 2:15) {  # Assuming you want to check from 2 to 15 clusters
  pam.fit <- pam(huh6.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:15
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 10
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:16))
## Draw the geom_vline at 11 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 9, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
huh6.k10.pam <- pam(huh6.jaccard.dist, k = 10, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- huh6.k10.pam
## Assign the clusters back to your Seurat objects
huh6.clusters.df <- anno.df$PAM.Cluster
names(huh6.clusters.df) <- rownames(anno.df)
uta.set <- paste0("UT_A_", as.character(huh6.list$Untreated_A$seurat_clusters.0.4))
utb.set <- paste0("UT_B_", as.character(huh6.list$Untreated_B$seurat_clusters.0.4))
reca.set <- paste0("Rec_A_", as.character(huh6.list$Recovered_A$seurat_clusters.0.4))
recb.set <- paste0("Rec_B_", as.character(huh6.list$Recovered_B$seurat_clusters.0.4))
all.set <- paste0("HuH6_all_", as.character(huh6.seurat$seurat_clusters.0.4))

huh6.list$Untreated_A$PAM.Cluster <- factor(as.character(huh6.clusters.df[uta.set]))
huh6.list$Untreated_B$PAM.Cluster <- factor(as.character(huh6.clusters.df[utb.set]))
huh6.list$Recovered_A$PAM.Cluster <- factor(as.character(huh6.clusters.df[reca.set]))
huh6.list$Recovered_B$PAM.Cluster <- factor(as.character(huh6.clusters.df[recb.set]))
huh6.seurat$PAM.Cluster <- factor(as.character(huh6.clusters.df[all.set]))
saveRDS(huh6.list, "data/huh6_seurat_split_list_pam.rds")
saveRDS(huh6.seurat, "data/huh6_seurat_hvgs_pam.rds")

## Prepare data for cNMF
dir.create("python/cNMF/counts/UT_A/", recursive = TRUE)
dir.create("python/cNMF/counts/UT_B/", recursive = TRUE)
dir.create("python/cNMF/counts/Rec_A/", recursive = TRUE)
dir.create("python/cNMF/counts/Rec_B/", recursive = TRUE)

Matrix::writeMM(LayerData(huh6.list$Untreated_A, layer = "counts"), "python/cNMF/counts/UT_A/matrix.mtx")
write.table(cbind(Features(huh6.list$Untreated_A), Features(huh6.list$Untreated_A)),
            file = "python/cNMF/counts/UT_A/genes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(Cells(huh6.list$Untreated_A),
            file = "python/cNMF/counts/UT_A/barcodes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

Matrix::writeMM(LayerData(huh6.list$Untreated_B, layer = "counts"), "python/cNMF/counts/UT_B/matrix.mtx")
write.table(cbind(Features(huh6.list$Untreated_B), Features(huh6.list$Untreated_B)),
            file = "python/cNMF/counts/UT_B/genes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(Cells(huh6.list$Untreated_B),
            file = "python/cNMF/counts/UT_B/barcodes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

Matrix::writeMM(LayerData(huh6.list$Recovered_A, layer = "counts"), "python/cNMF/counts/Rec_A/matrix.mtx")
write.table(cbind(Features(huh6.list$Recovered_A), Features(huh6.list$Recovered_A)),
            file = "python/cNMF/counts/Rec_A/genes.tsv",
            sep = "/t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(Cells(huh6.list$Recovered_A),
            file = "python/cNMF/counts/Rec_A/barcodes.tsv",
            sep = "/t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

Matrix::writeMM(LayerData(huh6.list$Recovered_B, layer = "counts"), "python/cNMF/counts/Rec_B/matrix.mtx")
write.table(cbind(Features(huh6.list$Recovered_B), Features(huh6.list$Recovered_B)),
            file = "python/cNMF/counts/Rec_B/genes.tsv",
            sep = "/t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
write.table(Cells(huh6.list$Recovered_B),
            file = "python/cNMF/counts/Rec_B/barcodes.tsv",
            sep = "/t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

## [ Signatures ] ----

## Other signatures from literature
huh6.harmony.seurat <- readRDS("data/huh6_seurat_hvgs_pam.rds")

## Removed Hooks_NT, Hooks_C2B and CR_14q32 as not enough features in object
sig.df <- as.data.frame(read_xlsx(path = "data/hb_sigs.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")

sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})

names(sigs) <- unique(sig.df$Signature)
sigs

sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)

huh6.harmony.seurat <- AddModuleScore(huh6.harmony.seurat, features = sigs, assay = "RNA", seed = 12345, 
                                      name = sig.names)
colnames(huh6.harmony.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(huh6.harmony.seurat@meta.data))

huh6.harmony.seurat[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = huh6.harmony.seurat, vars = sig.names)))

DoHeatmap(huh6.harmony.seurat, features = sig.names, assay = "HB_sigs", slot = "data", group.by = "Sample",
          group.colors = huh6.cols, label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Condition"))

ggsave("plots/sig/huh6_hb_sigs_heatmap.pdf", width = 5.8, height = 8.7)
ggsave("plots/sig/huh6_hb_sigs_heatmap.png", width = 5.8, height = 8.7)

DoHeatmap(huh6.harmony.seurat, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "PAM.Cluster", label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM cluster"))

ggsave("plots/sig/huh6_hb_sigs_heatmap_pam_clusters.pdf", width = 5.8, height = 8.7)
ggsave("plots/sig/huh6_hb_sigs_heatmap_pam_clusters.png", width = 5.8, height = 8.7)

## HNF4A and LEF1 (Kluiver et al., 2023)
huh6.harmony.seurat$Sample <- factor(huh6.harmony.seurat$Sample,
                                     levels = c("Untreated_A", "Untreated_B", "Recovered_A", "Recovered_B"))
Idents(huh6.harmony.seurat) <- huh6.harmony.seurat$Sample

features = c("HNF4A", "LEF1")

VlnPlot(huh6.harmony.seurat, features = features, pt.size = 0, cols = huh6.cols) &
  xlab("")
ggsave("plots/sig/huh6_hnf4a_lef1_violin.pdf", width = 8.7, height = 5.8)
ggsave("plots/sig/huh6_hnf4a_lef1_violin.png", width = 8.7, height = 5.8)

FeaturePlot(huh6.harmony.seurat, features = features) &
  umap.theme()
ggsave("plots/sig/huh6_hnf4a_lef1_umap.pdf", width = 8.7, height = 5.8)
ggsave("plots/sig/huh6_hnf4a_lef1_umap.png", width = 8.7, height = 5.8)

Idents(huh6.harmony.seurat) <- huh6.harmony.seurat$PAM.Cluster

VlnPlot(huh6.harmony.seurat, features = features, pt.size = 0) &
  xlab("")
ggsave("plots/sig/huh6_hnf4a_lef1_violin_pam_clusters.pdf", width = 8.7, height = 5.8)
ggsave("plots/sig/huh6_hnf4a_lef1_violin_pam_clusters.png", width = 8.7, height = 5.8)

## GRN modules in Roegrig et al. (2024)
## Could not find enough features in object

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

huh6.harmony.seurat <- AddModuleScore(huh6.harmony.seurat, features = H.list, assay = "RNA", seed = 12345, 
                                      name = H.names)
colnames(huh6.harmony.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(huh6.harmony.seurat@meta.data))

huh6.harmony.seurat[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = huh6.harmony.seurat, vars = H.names)))

DoHeatmap(huh6.harmony.seurat, features = H.names, assay = "MSigDB_H", slot = "data", group.by = "Sample",
          group.colors = huh6.cols, size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))

ggsave("plots/sig/huh6_msigdb_H_heatmap.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/huh6_msigdb_H_heatmap.png", width = 16.5, height = 8.7)

DoHeatmap(huh6.harmony.seurat, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))

ggsave("plots/sig/huh6_msigdb_H_heatmap_pam_clusters.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/huh6_msigdb_H_heatmap_pam_clusters.png", width = 16.5, height = 8.7)

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

huh6.harmony.seurat <- AddModuleScore(huh6.harmony.seurat, features = liver.list, assay = "RNA", seed = 12345, 
                                      name = liver.names)
colnames(huh6.harmony.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(huh6.harmony.seurat@meta.data))

huh6.harmony.seurat[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = huh6.harmony.seurat, vars = liver.names)))

DoHeatmap(huh6.harmony.seurat, features = liver/liver.names, assay = "MSigDB_C2", slot = "data", group.by = "Sample",
          group.colors = huh6.cols, size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))

ggsave("plots/sig/huh6_msigdb_C2_liver_heatmap.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/huh6_msigdb_C2_liver_heatmap.png", width = 16.5, height = 8.7)

DoHeatmap(huh6.harmony.seurat, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM cluster"))

ggsave("plots/sig/huh6_msigdb_C2_liver_heatmap_pam_clusters.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/huh6_msigdb_C2_liver_heatmap_pam_clusters.png", width = 16.5, height = 8.7)

saveRDS(huh6.harmony.seurat, "data/huh6_seurat_hvgs_pam_sigs.rds")

## [ Cluster annotation ] ----

huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_sigs.rds")

umap.gg <- DimPlot(huh6.seurat, group.by = "PAM.Cluster", order = TRUE) +
  umap.theme() + labs(title = "PAM clusters")
umap.gg
ggsave("plots/filter_hvgs/umap_harmony_t05_lnull_s005_pam_clusters.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/umap_harmony_t05_lnull_s005_pam_clusters.png", width = 8.3, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk across clusters and samples //////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Pseudobulk across clusters and samples and look at expression of same pathways/markers as above
pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Sample", "PAM.Cluster"))

Idents(pseudo.huh6) <- "PAM.Cluster"

bulk.markers <- FindAllMarkers(pseudo.huh6, test.use = "DESeq2")
saveRDS(bulk.markers, "data/huh6_pseudobulk_pam_markers.rds")

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

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = H.list, assay = "RNA", seed = 12345, 
                              name = H.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = H.names)))

DoHeatmap(pseudo.huh6, features = H.names, assay = "MSigDB_H", slot = "data", group.by = "Sample",
          group.colors = huh6.cols, size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_huh6_msigdb_H_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_H_heatmap.png", width = 11.7, height = 8.3)

DoHeatmap(pseudo.huh6, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_huh6_msigdb_H_heatmap_pam_clusters.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_H_heatmap_pam_clusters.png", width = 11.7, height = 8.3)

huh6.h.mat <- as.matrix(pseudo.huh6@assays$MSigDB_H@data)

anno.df <- pseudo.huh6@meta.data["Sample"]
gt <- pheatmap(huh6.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(huh6.cols))$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_H_pheatmap.pdf", plot = gt)
ggsave("plots/sig/pseudo_huh6_msigdb_H_pheatmap.png", plot = gt)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_H_pheatmap_pam_clusters.pdf", plot = gt)
ggsave("plots/sig/pseudo_huh6_msigdb_H_pheatmap_pam_clusters.png", plot = gt)

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

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = liver.list, assay = "RNA", seed = 12345, 
                              name = liver.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = liver.names)))

DoHeatmap(pseudo.huh6, features = liver.names, assay = "MSigDB_C2", slot = "data", group.by = "Sample",
          group.colors = huh6.cols, size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 3)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_huh6_msigdb_C2_liver_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_liver_heatmap.png", width = 11.7, height = 8.3)

DoHeatmap(pseudo.huh6, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_huh6_msigdb_C2_liver_heatmap_pam_clusters.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_liver_heatmap_pam_clusters.png", width = 11.7, height = 8.3)

huh6.c2.mat <- as.matrix(pseudo.huh6@assays$MSigDB_C2@data)
huh6.reactome.mat <- huh6.c2.mat[which(grepl("REACTOME", rownames(huh6.c2.mat))), ]

anno.df <- pseudo.huh6@meta.data["Sample"]
gt <- pheatmap(huh6.reactome.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(huh6.cols))$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.reactome.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_C2_reactome_pheatmap_pam_clusters.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_reactome_pheatmap_pam_clusters.png", plot = gt, width = 11.7, height = 8.3)

huh6.cairo.mat <- huh6.c2.mat[which(grepl("CAIRO", rownames(huh6.c2.mat))), ]

anno.df <- pseudo.huh6@meta.data["Sample"]
gt <- pheatmap(huh6.cairo.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(huh6.cols))$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_C2_cairo_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_cairo_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.cairo.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_huh6_msigdb_C2_cairo_pheatmap_pam_clusters.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_huh6_msigdb_C2_cairo_pheatmap_pam_clusters.png", plot = gt, width = 11.7, height = 8.3)

saveRDS(pseudo.huh6, "data/huh6_seurat_pseudobulk_pam_sigs.rds")

## Liver differentiation markers (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/)
custom <- c("HNF4A", "EPCAM", "KRT19", "AFP", "PROM1", "SOX9", "NCAM1", "ICAM1")
custom %in% rownames(pseudo.huh6)

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])

anno.df <- pseudo.huh6@meta.data["Sample"]
gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(huh6.cols))$gtable
ggsave("plots/sig/pseudo_huh6_custom_pheatmap.pdf", plot = gt)
ggsave("plots/sig/pseudo_huh6_custom_pheatmap.png", plot = gt)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_huh6_custom_pheatmap_pam_clusters.pdf", plot = gt)
ggsave("plots/sig/pseudo_huh6_custom_pheatmap_pam_clusters.png", plot = gt)

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk across clusters //////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Pseudobulk across clusters only and look at expression of same pathways/markers as above
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_sigs.rds")

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("PAM.Cluster"))

Idents(pseudo.huh6) <- "PAM.Cluster"

bulk.markers <- FindAllMarkers(pseudo.huh6, test.use = "DESeq2")
saveRDS(bulk.markers, "data/huh6_pseudobulk_clusters_pam_markers.rds")

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

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = H.list, assay = "RNA", seed = 12345, 
                              name = H.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = H.names)))

DoHeatmap(pseudo.huh6, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_H_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_H_heatmap.png", width = 11.7, height = 8.3)

huh6.h.mat <- as.matrix(pseudo.huh6@assays$MSigDB_H@data)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.h.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_H_pheatmap.pdf", plot = gt)
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_H_pheatmap.png", plot = gt)

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

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = liver.list, assay = "RNA", seed = 12345, 
                              name = liver.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = liver.names)))

DoHeatmap(pseudo.huh6, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_liver_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_liver_heatmap.png", width = 11.7, height = 8.3)

huh6.c2.mat <- as.matrix(pseudo.huh6@assays$MSigDB_C2@data)
huh6.reactome.mat <- huh6.c2.mat[which(grepl("REACTOME", rownames(huh6.c2.mat))), ]

gt <- pheatmap(huh6.reactome.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

huh6.cairo.mat <- huh6.c2.mat[which(grepl("CAIRO", rownames(huh6.c2.mat))), ]

gt <- pheatmap(huh6.cairo.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_cairo_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_clusters_huh6_msigdb_C2_cairo_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

## Liver differentiation markers
## (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/; https://pmc.ncbi.nlm.nih.gov/articles/PMC11114060/)
custom <- c("MCAM", ## Generally not expressed in normal liver, high in HCC; PTPRC = CD45 expressed in HPCs
            "FOXA1", "FOXA2", "GATA4", ## Early hepatic specification
            "EPCAM", "NCAM1", "PROM1", ## HSC markers; EPCAM also HB marker
            "KRT8", "KRT18", "KRT19", ## HSC markers
            "CXCR4", ## Endoderm marker expressed in HSCs and HCC
            "AFP", "ICAM1", ## HB markers; AFP also fetal hepatocyte; KRT7 (CK7) also cholangiocyte marker
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
custom %in% rownames(pseudo.huh6)

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])

gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_pam_clusters_huh6_custom_pheatmap.pdf", plot = gt)
ggsave("plots/sig/pseudo_pam_clusters_huh6_custom_pheatmap.png", plot = gt)

## /////////////////////////////////////////////////////////////////////////////
## PAM cluster markers /////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Differential analysis between PAM clusters
Idents(huh6.seurat) <- huh6.seurat$PAM.Cluster
## This does a 1 vs All Differential analysis
mast.markers <- FindAllMarkers(huh6.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

table(mast.markers$cluster)
## 1   10    2    3    4    8    9 
## 4050 7415 3145 4583 2639 5409 4455 

mast.list <- lapply(unique(mast.markers$cluster), function(x){
  tmp.df <- mast.markers[mast.markers$cluster==x & mast.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})
names(mast.list) <- unique(mast.markers$cluster)

saveRDS(mast.list, "data/final_markers/huh6_pam_markers_list.rds")

mast.GO.BP.list <- lapply(mast.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

mast.C2.list <- lapply(mast.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

mast.H.list <- lapply(mast.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(mast.GO.BP.list, "data/final_markers/huh6_pam_markers_BP_ontology.rds")
saveRDS(mast.C2.list, "data/final_markers/huh6_pam_markers_C2_geneset.rds")
saveRDS(mast.H.list, "data/final_markers/huh6_pam_markers_Hallmarks.rds")

mast.lists <- list("mast.GO.BP.list" = mast.GO.BP.list,
                   "mast.C2.list" = mast.C2.list,
                   "mast.H.list" = mast.H.list,
                   "mast.list" = mast.list)

mast.names <- sub("list", "sheets", names(mast.lists))

for(i in 1:length(mast.lists)){
  mast.list <- mast.lists[[i]]
  list.sheets <- lapply(mast.list, as.data.frame)
  assign(mast.names[[i]], list.sheets)
}

write_xlsx(mast.GO.BP.sheets, "data/final_markers/huh6_pam_markers_BP_ontology.xlsx")
write_xlsx(mast.C2.sheets, "data/final_markers/huh6_pam_markers_C2_geneset.xlsx")
write_xlsx(mast.H.sheets, "data/final_markers/huh6_pam_markers_Hallmarks.xlsx")
write_xlsx(mast.sheets, "data/final_markers/huh6_pam_markers.xlsx")

## Plot expression of PAM cluster markers from scRNA-seq data
pam.list <- readRDS("data/huh6_pam_markers_list_scrna.rds")

pam.filt <- lapply(pam.list, function(x) filter(x, avg_log2FC >= 1))
names(pam.filt) <- recode(as.character(names(pam.filt)),
                          "1" = "Hepatocytic HNF1b+",
                          "2" = "HSC-like",
                          "3" = "Hepato+cholangiocytic",
                          "4" = "Fetal hepatocytic",
                          "5" = "HB-like",
                          "7" = "HB+HSC-like")
pam.genes <- lapply(pam.filt, function(x) x[, "gene"])
pam.names <- paste(names(pam.genes), ".Sig", sep = "")

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = pam.genes, assay = "RNA", seed = 12345, 
                              name = pam.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["PAM_scRNAseq"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = pam.names)))

DoHeatmap(pseudo.huh6, features = pam.names, assay = "PAM_scRNAseq", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_pam_clusters_huh6_scrna_markers_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_clusters_huh6_scrna_markers_heatmap.png", width = 11.7, height = 8.3)

huh6.scrna.mat <- as.matrix(pseudo.huh6@assays$PAM_scRNAseq@data)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.scrna.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/sig/pseudo_pam_clusters_huh6_scrna_markers_pheatmap.pdf", plot = gt)
ggsave("plots/sig/pseudo_pam_clusters_huh6_scrna_markers_pheatmap.png", plot = gt)

saveRDS(pseudo.huh6, "data/huh6_seurat_pseudobulk_pam_clusters_sigs.rds")

## Jaccard similarity between scRNA-seq and snRNA-seq clusters
scrna.markers <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/jaccard/huh6_pam_markers.rds")
snrna.markers <- readRDS("data/final_markers/huh6_pam_markers_list.rds")
snrna.markers <- bind_rows(snrna.markers)

## Filter for p-value and log2FC
scrna.markers <- scrna.markers[scrna.markers$p_val_adj < 0.05, ]
scrna.markers <- scrna.markers[scrna.markers$avg_log2FC > 0.5, ]
snrna.markers <- snrna.markers[snrna.markers$p_val_adj < 0.05, ]
snrna.markers <- snrna.markers[snrna.markers$avg_log2FC > 0.5, ]

scrna.markers$name <- recode(as.character(scrna.markers$cluster),
                            "1" = "Hepatocytic 2",
                            "2" = "Stem-like",
                            "3" = "Hepatocytic 3",
                            "4" = "Hepatocytic 1",
                            "5" = "Late progenitor",
                            "7" = "Early progenitor")

scrna.markers <- split(scrna.markers, scrna.markers$name)
snrna.markers <- split(snrna.markers, snrna.markers$cluster)

lapply(scrna.markers, nrow)
## Early progenitor = 869, Hepatocytic 1 = 1121, Hepatocytic 2 = 441,
## Hepatocytic 3 = 548, Late progenitor = 1348, Stem-like = 199
lapply(snrna.markers, nrow)
## 1 = 288, 10 = 3073, 2 = 282, 3 = 361,
## 4 = 619, 8 = 979, 9 = 1117

## Take top 200 per cluster
scrna.markers <- lapply(scrna.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
snrna.markers <- lapply(snrna.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})

hb.markers.list <- c(scrna.markers, snrna.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}

hb.jaccard.mat <- matrix(data = NA, nrow = length(hb.markers.list),
                         ncol = length(hb.markers.list),
                         dimnames = list(names(hb.markers.list),
                                         names(hb.markers.list)))

## Run your pairwise Jaccard similarity
for (i in rownames(hb.jaccard.mat)) {
  for (j in colnames(hb.jaccard.mat)) {
    hb.jaccard.mat[i,j] <- jaccard(hb.markers.list[[i]]$gene, hb.markers.list[[j]]$gene)
  }
}

dev.off()
gt <- pheatmap(hb.jaccard.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("plots/jaccard/jaccard_heatmap_clusters.png", plot = gt)
ggsave("plots/jaccard/jaccard_heatmap_clusters.pdf", plot = gt)

## Annotate clusters
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam.rds")
huh6.seurat$PAM_Name <- recode(as.character(huh6.seurat$PAM.Cluster),
                               "1" = "Hepatocytic 3",
                               "9" = "Late progenitor",
                               "2" = "Hepatocytic 2",
                               "3" = "Hepatocytic 1",
                               "4" = "Late progenitor",
                               "8" = "Early progenitor",
                               "10" = "Late progenitor")
saveRDS(huh6.seurat, "data/huh6_seurat_hvgs_pam_anno.rds")

## Cluster representation in each sample
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")

umap.gg <- DimPlot(huh6.seurat, group.by = "PAM_Name", order = TRUE) +
  umap.theme() +
  labs(title = "PAM clusters") +
  scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15))
umap.gg
ggsave("plots/sig/umap_harmony_pam_clusters_anno.pdf", width = 8.3, height = 5.8)
ggsave("plots/sig/umap_harmony_pam_clusters_anno.png", width = 8.3, height = 5.8)

anno.table <- as.data.frame(table(huh6.seurat$Sample, huh6.seurat$PAM_Name))

anno.table$Var1 <- factor(anno.table$Var1,
                          levels = c("Untreated_A", "Untreated_B",
                                     "Recovered_A", "Recovered_B"))

ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "PAM cluster name")) +
  scale_fill_manual(values = pam.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/sig/pam_cluster_proportion_condition.pdf", width = 8.3, height = 5.8)
ggsave("plots/sig/pam_cluster_proportion_condition.png", width = 8.3, height = 5.8)
