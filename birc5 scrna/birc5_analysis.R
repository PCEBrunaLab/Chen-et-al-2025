## BIRC5 scRNA-seq analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

renv::load()

library(scater)
library(scran)
library(Seurat)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(msigdbr)
library(cluster)
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggsankey)
library(pheatmap)

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

## Colours for experimental condition
group.cols <- c("NTC" = "#C4C4C4",
                "BIRC5_siRNA" = "#9FC3E6",
                "E2F4_siRNA" = "#F0AFA0")

## Colours for cell line identity
cell.cols <- c("HB303" = "#C7B6E6",
               "HepG2" = "#9ED6C1",
               "HuH6" = "#F3D7A3")

## [ Data preparation ] ----

hb303.sce <- readRDS("data/HB_BIRC5_SCE-norm_relaxed_HB303_norm.RDS")
hepg2.sce <-readRDS("data/HB_BIRC5_SCE-norm_relaxed_HepG2_norm.RDS")
huh6.sce <- readRDS("data/HB_BIRC5_SCE-norm_relaxed_HuH6_norm.RDS")

sce.all <- cbind(hb303.sce, hepg2.sce, huh6.sce)
saveRDS(sce.all, "data/birc5_sce_combined.rds")

## [ Initial filtering ] ----

sce.all <- readRDS("data/birc5_sce_combined.rds")

ensembl <- useMart("ensembl", "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
ens.bm <- getBM(attributes = c("ensembl_gene_id",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "external_gene_name",
                               "hgnc_symbol",
                               "description"),
                mart = ensembl,
                filters = "ensembl_gene_id",
                values = rownames(sce.all))
write.csv(ens.bm, "data/ensembl_biomart.csv")

## Keep track of any duplicates
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)]
## Make a vector and give elements names at the same time
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])
## Duplicates: LINC00595

## Remove unnecessary annotation from the description
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)
## This filters out around 1000 genes
nrow(ens.filt.bm) ## 35456
nrow(sce.all) ## 36601

## Keep only the simple chromosome information
table(ens.filt.bm$chromosome_name)
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]
## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]
## Now we've filtered out around 10,000 genes
nrow(ens.filt.bm) ## 25577
nrow(sce.all) ## 36601

## Keep genes that didn't get annotated so we know what gets dropped
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 9879
length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 9733 / 11082 are novel transcripts or novel transcripts antisense to a gene
grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest is a mixed bag: either novel, mitochondrially encoded, non-coding, or just empty description
## It seems unlikely to be particularly informative in our case

## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## ELFN2 = extracellular leucine rich repeat and fibronectin type III domain containing 2
## GOLGA8M = golgin A8 family member M
## RAET1E-AS1 = RAET1E antisense RNA 1
## SPATA13 = spermatogenesis associated 13

## LINC00595
## LINC01115
## LINC03021
## LINC03023
## LINC03025

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 25568

mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype=="lincRNA"]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Before we start chopping and changing, let's get some basic QC into the object
sce.all <- addPerCellQCMetrics(sce.all, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))
sce.all <- sce.all[ens.filt.bm$ensembl_gene_id, ]
rownames(sce.all) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

metadata.df <-  as.data.frame(colData(sce.all))
birc5.seurat <- CreateSeuratObject(counts = assay(sce.all, "counts"),
                                   assay = "RNA",
                                   meta.data = metadata.df)

saveRDS(sce.all, "data/birc5_sce_start.rds")
saveRDS(birc5.seurat, "data/birc5_seurat_start.rds")

## [ QC ] ----

birc5.sce <- readRDS("data/birc5_sce_start.rds")

dir.create("plots/qc/", recursive = TRUE)

table(birc5.sce$Condition)
## BIRC5_siRNA  E2F4_siRNA         NTC 
##        8167        6383       22635 
birc5.sce$Condition <- factor(birc5.sce$Condition,
                              levels = c("NTC", "BIRC5_siRNA", "E2F4_siRNA"))
birc5.sce$Sample_ID <- factor(birc5.sce$Sample_ID,
                              levels = c("HB303_NTC_B2", "HB303_NTC_B3",
                                         "HB303_BIRC5_siRNA_A2", "HB303_BIRC5_siRNA_A3",
                                         "HB303_E2F4_siRNA_B2", "HB303_E2F4_siRNA_B3",
                                         "HepG2_NTC_B2", "HepG2_NTC_B3",
                                         "HepG2_BIRC5_siRNA_A2", "HepG2_BIRC5_siRNA_A3",
                                         "HuH6_NTC_B2", "HuH6_NTC_B3",
                                         "HuH6_BIRC5_siRNA_A2", "HuH6_BIRC5_siRNA_A3",
                                         "HuH6_E2F4_siRNA_B2", "HuH6_E2F4_siRNA_B3"))

det.sce.gg <- plotColData(birc5.sce, y = "detected", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
umi.sce.gg <- plotColData(birc5.sce, y = "total", x = "Sample_ID", colour_by = "Condition") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")
mt.sce.gg <- plotColData(birc5.sce, y = "subsets_mt_percent", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
rb.sce.gg <- plotColData(birc5.sce, y = "subsets_ribo_percent", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.5)
ggsave("plots/qc/birc5_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/birc5_qc_plots.png", width = 8.3, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
## Use sparse matrix operations, as if your dataset is large, doing matrix operations the regular way will take a long time
counts.assay <- counts(birc5.sce)
counts.assay@x <- counts.assay@x/rep.int(colSums(counts.assay), diff(counts.assay@p))
top.expr <- order(Matrix::rowSums(counts.assay), decreasing = TRUE)[20:1]
top.expr <- as.matrix(t(counts.assay[top.expr, ]))
top.expr <- as.data.frame(top.expr)
ord.idx <- colnames(top.expr)
top.expr <- melt(top.expr, variable.name = "Gene", value.name = "Prop")
top.expr$Gene <- factor(top.expr$Gene, levels = ord.idx)

ggplot(top.expr,
       mapping = aes(x = Gene, y = Prop * 100)) +
  geom_boxplot() +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dotted"),
                     text = element_text(size = 15)) +
  coord_flip() +
  labs(y = "% total expression", x = "Gene", title = "Percentage expression of a single gene per total cell expression")
ggsave("plots/qc/birc5_topgenes.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/birc5_topgenes.png", width = 8.3, height = 5.8)

## Remove the counts as it's not needed.
rm(counts.assay)
gc()

## MALAT1 is a long non-coding RNA and typically comes up enriched in skim-seq like 10x
## It is routinely removed as part of preprocessing
## However, since we're doing drug treatments, it's worth checking to see if it might indicate stress response
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7099402/

malat1.sce.gg <- plotExpression(birc5.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Sample_ID", 
                                colour_by = "Condition") +
  labs(x = element_blank(), y = expression("Log"[2]*" Expression")) +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

lincrna.sce.gg <- plotColData(birc5.sce, y = "subsets_linc_percent", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.5, text = element_text(size = 15))
ggsave("plots/qc/birc5_lincrna_malat1.pdf", width = 8.3, height = 8.7)
ggsave("plots/qc/birc5_lincrna_malat1.png", width = 8.3, height = 8.7)

gc()

## Cell cycle prediction
## Tricycle assigns cell cycle along a radial trajectory between 0 and 2pi
library(tricycle)
tricycle.sce <- project_cycle_space(birc5.sce, gname.type = "SYMBOL", species = "human")
tricycle.sce <- estimate_cycle_position(tricycle.sce)

## TOP2a is a sanity check as it should peak at a specific part of the cell cycle that we can confirm is consistent with estimation
top2a.idx <- which(rowData(tricycle.sce)$Symbol == "TOP2A")
fit.l <- fit_periodic_loess(tricycle.sce$tricyclePosition,
                            assay(tricycle.sce, 'logcounts')[top2a.idx, ],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", 
                            y_lab = "log2(TOP2A)",
                            fig.title = paste0("Expression of TOP2A along \u03b8 (n=", 
                                               ncol(tricycle.sce), 
                                               ")"))
fit.l$fig + theme_bw(base_size = 14)
## TOP2A should peak at pi - it does?

## From the tricycle guide: users can approximately relate:
## 0.5pi to be the start of S stage, 
## pi to be the start of G2M stage, 
## 1.5pi to be the middle of M stage
## and 1.75pi-0.25pi to be G1/G0 stage

## Make that into some kind of maths
score.tricycle <- function(x){
  ifelse(x >= 0.5*pi & x < pi, "S", 
         ifelse(x >= pi & x < 1.75*pi, "G2M", 
                ifelse(x >= 1.75*pi | x < 0.25*pi, "G1", "NA")))
}

birc5.seurat <- readRDS("data/birc5_seurat_start.rds")
birc5.cycle.seurat <- NormalizeData(birc5.seurat)
birc5.cycle.seurat <- CellCycleScoring(birc5.cycle.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

cell.cycle.df <- data.frame("Cell.ID" = colnames(birc5.sce),
                            "Seurat" = birc5.cycle.seurat$Phase,
                            "Tricycle" = score.tricycle(tricycle.sce$tricyclePosition))
cell.cycle.melt <- make_long(cell.cycle.df, Seurat, Tricycle)

ggplot(cell.cycle.melt, 
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
  labs(x = element_blank(), title = "Cell Cycle Prediction") +
  scale_fill_manual(values = hues::iwanthue(4))
ggsave("plots/qc/cell_cycle_comparison.pdf", width = 5.8, height = 5.8)
ggsave("plots/qc/cell_cycle_comparison.png", width = 5.8, height = 5.8)

## Both results roughly agree
## Tricycle is newest so maybe we go with that?

## Add cell cycle results back into the object
birc5.sce$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
birc5.sce$Tricycle.Position <- tricycle.sce$tricyclePosition
birc5.sce$Seurat.Phase <- birc5.cycle.seurat$Phase
birc5.sce$Seurat.S <- birc5.cycle.seurat$S.Score
birc5.sce$Seurat.G2M <- birc5.cycle.seurat$G2M.Score

birc5.seurat$Seurat.Phase <- birc5.cycle.seurat$Phase
birc5.seurat$Seurat.S <- birc5.cycle.seurat$S.Score
birc5.seurat$Seurat.G2M <- birc5.cycle.seurat$G2M.Score
birc5.seurat$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
birc5.seurat$Tricycle.Position <- tricycle.sce$tricyclePosition

rm(tricycle.sce, birc5.cycle.seurat)
gc()

saveRDS(birc5.sce, "data/birc5_sce_qc.rds")
saveRDS(birc5.seurat, "data/birc5_seurat_qc.rds")

## [ QC filtering ] ----

birc5.sce <- readRDS("data/birc5_sce_qc.rds")
birc5.seurat <- readRDS("data/birc5_seurat_qc.rds")

colnames(colData(birc5.sce))
colnames(rowData(birc5.sce))

table(birc5.sce$Sample_ID, birc5.sce$Condition)

det.sce.gg <- plotColData(birc5.sce, y = "detected", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dotted")
mt.sce.gg <- plotColData(birc5.sce, y = "subsets_mt_percent", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 20, linetype = "dotted")
rb.sce.gg <- plotColData(birc5.sce, y = "subsets_ribo_percent", x = "Sample_ID", colour_by = "Condition") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dotted")

det.sce.gg + mt.sce.gg + rb.sce.gg & theme(aspect.ratio = 0.5)
ggsave("plots/qc/birc5_qc_plots_filt.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/birc5_qc_plots_filt.png", width = 8.3, height = 5.8)

## Set the expression threshold over 1000
detected.filt <- colnames(birc5.sce)[birc5.sce$detected > 1000] ## Doesn't remove any cells
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(birc5.sce)[Matrix::rowSums(counts(birc5.sce)) > 5]
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(birc5.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

filt.sce <- birc5.sce[genes.filt, detected.filt]
dim(birc5.sce) - dim(filt.sce)
## We lose ~ 5000 genes

mito.filt <- colnames(filt.sce)[filt.sce$subsets_mt_percent < 20]
ribo.filt <- colnames(filt.sce)[filt.sce$subsets_ribo_percent > 5]
qc.filt <- intersect(mito.filt, ribo.filt)

filt.sce <- filt.sce[, qc.filt]
dim(birc5.sce) - dim(filt.sce)
## We lose ~ 400 cells

## Percentage of the library that we filter out
round(100 - (100 * (table(filt.sce$Sample_ID) / table(birc5.sce$Sample_ID))), 1)
## HB303_NTC_B2         HB303_NTC_B3 HB303_BIRC5_siRNA_A2 HB303_BIRC5_siRNA_A3  HB303_E2F4_siRNA_B2 
##          0.1                  0.0                  0.0                  0.0                  0.0 
## HB303_E2F4_siRNA_B3         HepG2_NTC_B2         HepG2_NTC_B3 HepG2_BIRC5_siRNA_A2 HepG2_BIRC5_siRNA_A3 
##                 0.0                  1.3                  0.0                  0.0                  0.0 
## HuH6_NTC_B2          HuH6_NTC_B3  HuH6_BIRC5_siRNA_A2  HuH6_BIRC5_siRNA_A3   HuH6_E2F4_siRNA_B2 
##         0.2                  4.2                  0.3                  1.5                  1.4 
## HuH6_E2F4_siRNA_B3 
##                1.1 

table(filt.sce$Sample_ID)
## HB303_NTC_B2         HB303_NTC_B3 HB303_BIRC5_siRNA_A2 HB303_BIRC5_siRNA_A3  HB303_E2F4_siRNA_B2 
##         8659                  176                   70                 2683                  605 
## HB303_E2F4_siRNA_B3         HepG2_NTC_B2         HepG2_NTC_B3 HepG2_BIRC5_siRNA_A2 HepG2_BIRC5_siRNA_A3 
##                  29                 6522                 1522                  148                    5 
## HuH6_NTC_B2          HuH6_NTC_B3  HuH6_BIRC5_siRNA_A2  HuH6_BIRC5_siRNA_A3   HuH6_E2F4_siRNA_B2 
##        2702                 2829                 1153                 4042                 3190 
## HuH6_E2F4_siRNA_B3 
##               2487 

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(birc5.seurat)[Matrix::rowSums(GetAssayData(birc5.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(birc5.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

## Seurat and Scater have their own methods for calling # of detected genes so these will disagree slightly
## To be consistent we'll filter on $detected rather than $nFeature_RNA
birc5.seurat$nFeature_RNA == birc5.seurat$detected
filt.seurat <- subset(birc5.seurat, 
                      cells = WhichCells(birc5.seurat, 
                                         expression = subsets_mt_percent < 20 &
                                           subsets_ribo_percent > 5 &
                                           detected > 1000),
                      features = genes.filt)
## The dimensions are equal
dim(filt.seurat) == dim(filt.sce)

saveRDS(filt.sce, "data/birc5_sce_filt.rds")
saveRDS(filt.seurat, "data/birc5_seurat_filt.rds")

## [ Filter HVGs ] ----

filt.seurat <- readRDS("data/birc5_seurat_filt.rds")

## Yura's script to filter HVGs
filt.seurat <- NormalizeData(filt.seurat)
filt.seurat <- FindVariableFeatures(filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(filt.seurat), 20), repel = TRUE)
ggsave("plots/filter_hvgs/filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/filter_hvgs.png", width = 8.3, height = 5.8)

## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(filt.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(filt.seurat), value = TRUE),
              grep("[\\.]", rownames(filt.seurat), value = TRUE),
              grep("^LINC", rownames(filt.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
VariableFeatures(filt.seurat) <- VariableFeatures(filt.seurat)[which(!VariableFeatures(filt.seurat) %in% rm.genes)]
## Scale data and score cell cycle
filt.seurat <- ScaleData(filt.seurat, features = VariableFeatures(filt.seurat))

tmp.seurat <- CellCycleScoring(object = filt.seurat, 
                               g2m.features = cc.genes$g2m.genes, 
                               s.features = cc.genes$s.genes)
tmp.seurat$Cycle.Score <- tmp.seurat$S.Score - tmp.seurat$G2M.Score

seurat.cycle.melt <- make_long(data.frame("Sample" = filt.seurat$Sample, 
                                          "SCTransform" = tmp.seurat$Phase,
                                          "Normal" = filt.seurat$Seurat.Phase),
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

DefaultAssay(filt.seurat) <- "RNA"
filt.seurat$Seurat.Cycle.Score <- filt.seurat$Seurat.S - filt.seurat$Seurat.G2M

saveRDS(filt.seurat, "data/birc5_seurat_hvgs.rds")

birc5.seurat <- readRDS("data/birc5_seurat_hvgs.rds")

birc5.seurat <- RunPCA(birc5.seurat, verbose = FALSE, npcs = 100, 
                       features = VariableFeatures(birc5.seurat))
ElbowPlot(object = birc5.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate UMAPs
library(harmony)
harmony.seurat <- RunHarmony(birc5.seurat, group.by.vars = "Sample",
                             theta = 1, lambda = 1, sigma = 0.07,
                             assay.use = "RNA", reduction = "pca",
                             dims.use = 1:30, reduction.save = "Harmony",
                             max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = group.cols) +
  umap.theme() + labs(title = "Condition")
umap.gg
ggsave("plots/filter_hvgs/umap_harmony_t1_l1_s007.png", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/umap_harmony_t1_l1_s007.pdf", width = 8.3, height = 5.8)

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
saveRDS(harmony.seurat, "data/birc5_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
birc5.seurat <- readRDS("data/birc5_seurat_hvgs.rds")
obj.list <- SplitObject(birc5.seurat, split.by = "Cell_Line")

## Rerun variable feature selection
obj.list <- lapply(obj.list, FindVariableFeatures)
## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(birc5.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(birc5.seurat), value = TRUE),
              grep("[\\.]", rownames(birc5.seurat), value = TRUE),
              grep("^LINC", rownames(birc5.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(birc5.seurat, obj.list, obj)

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
saveRDS(output, "data/birc5_seurat_split_list.rds")

birc5.list <- readRDS("data/birc5_seurat_split_list.rds")
birc5.seurat <- readRDS("data/birc5_seurat_hvgs_harmony.rds")
## Look through UMAPs and choose best cluster resolution
# DimPlot(birc5.seurat, group.by = "seurat_clusters.0.6", order = TRUE)
Idents(birc5.list$HB303) <- birc5.list$HB303$seurat_clusters.0.4
Idents(birc5.list$HepG2) <- birc5.list$HepG2$seurat_clusters.0.4
Idents(birc5.list$HuH6) <- birc5.list$HuH6$seurat_clusters.0.4
Idents(birc5.seurat) <- birc5.seurat$seurat_clusters.0.6

## Get cluster markers - MAST takes longer to run
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
hb303.markers <- FindAllMarkers(birc5.list$HB303, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
hepg2.markers <- FindAllMarkers(birc5.list$HepG2, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
huh6.markers <- FindAllMarkers(birc5.list$HuH6, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
all.markers <- FindAllMarkers(birc5.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

dir.create("data/jaccard/", recursive = TRUE)
saveRDS(hb303.markers, "data/jaccard/hb303_markers.rds")
saveRDS(hepg2.markers, "data/jaccard/hepg2_markers.rds")
saveRDS(huh6.markers, "data/jaccard/huh6_markers.rds")
saveRDS(all.markers, "data/jaccard/all_markers.rds")

## Filter significant and higher expressed markers
hb303.markers <- hb303.markers[hb303.markers$p_val_adj < 0.05, ]
hepg2.markers <- hepg2.markers[hepg2.markers$p_val_adj < 0.05, ]
huh6.markers <- huh6.markers[huh6.markers$p_val_adj < 0.05, ]
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]

hb303.markers <- hb303.markers[hb303.markers$avg_log2FC > 0.5, ]
hepg2.markers <- hepg2.markers[hepg2.markers$avg_log2FC > 0.5, ]
huh6.markers <- huh6.markers[huh6.markers$avg_log2FC > 0.5, ]
all.markers <- all.markers[all.markers$avg_log2FC > 0.5, ]
## Split into lists of genes per cluster
hb303.markers <- split(hb303.markers, hb303.markers$cluster)
hepg2.markers <- split(hepg2.markers, hepg2.markers$cluster)
huh6.markers <- split(huh6.markers, huh6.markers$cluster)
all.markers <- split(all.markers, all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(hb303.markers, nrow) ## 0 = 35, 1 = 743, 2 = 114, 3 = 245, 4 = 40, 5 = 331
lapply(hepg2.markers, nrow) ## 0 = 73, 1 = 267, 2 = 313, 3 = 132, 4 = 362
lapply(huh6.markers, nrow) ## 0 = 338, 1 = 187, 2 = 422, 3 = 149, 4 = 285, 5 = 276, 6 = 1335, 7 = 1284, 8 = 499
lapply(all.markers, nrow) ## 0 = 1170, 1 = 920, 2 = 1022, 3 = 586, 4 = 966, 5 = 1290, 6 = 1109, 7 = 1121, 8 = 994, 9 = 892, 10 = 688, 11 = 1093

hb303.markers <- lapply(hb303.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
hepg2.markers <- lapply(hepg2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
huh6.markers <- lapply(huh6.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
all.markers <- lapply(all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})

## Give your clusters unique names
names(hb303.markers) <- paste0("HB303_", names(hb303.markers))
names(hepg2.markers) <- paste0("HepG2_", names(hepg2.markers))
names(huh6.markers) <- paste0("HuH6_", names(huh6.markers))
names(all.markers) <- paste0("All_", names(all.markers))
## Make a big list of all the cluster markers
markers.list <- c(hb303.markers, hepg2.markers, huh6.markers, all.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
jaccard.mat <- matrix(data = NA, nrow = length(markers.list),
                      ncol = length(markers.list),
                      dimnames = list(names(markers.list),
                                      names(markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(jaccard.mat)) {
  for (j in colnames(jaccard.mat)) {
    jaccard.mat[i,j] <- jaccard(markers.list[[i]]$gene, markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Cell_Line" = gsub("_.+", "", colnames(jaccard.mat)),
                      row.names = colnames(jaccard.mat))
anno.col <- list("Cell_Line" = c("All" = "black",
                                 "HB303" = as.character(cell.cols[1]),
                                 "HepG2" = as.character(cell.cols[2]),
                                 "HuH6" = as.character(cell.cols[3])))
dev.off()
gt <- pheatmap(jaccard.mat,
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
ggsave("plots/jaccard/jaccard_heatmap.png", plot = gt, height = 8.3, width = 8.3)
ggsave("plots/jaccard/jaccard_heatmap.pdf", plot = gt, height = 8.3, width = 8.3)

## /////////////////////////////////////////////////////////////////////////////
## PAM clustering //////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Run the parameter search for PAM
jaccard.dist <- as.dist(1 - jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
for (k in 2:15) {  # Assuming you want to check from 2 to 15 clusters
  pam.fit <- pam(jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:15
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 9
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:16))
## Draw the geom_vline at 9 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 8, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
k9.pam <- pam(jaccard.dist, k = 9, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- k9.pam
## Assign the clusters back to your Seurat objects
clusters.df <- anno.df$PAM.Cluster
names(clusters.df) <- rownames(anno.df)
hb303.set <- paste0("HB303_", as.character(birc5.list$HB303$seurat_clusters.0.4))
hepg2.set <- paste0("HepG2_", as.character(birc5.list$HepG2$seurat_clusters.0.4))
huh6.set <- paste0("HuH6_", as.character(birc5.list$HuH6$seurat_clusters.0.4))
all.set <- paste0("All_", as.character(birc5.seurat$seurat_clusters.0.6))

birc5.list$HB303$PAM.Cluster <- factor(as.character(clusters.df[hb303.set]))
birc5.list$HepG2$PAM.Cluster <- factor(as.character(clusters.df[hepg2.set]))
birc5.list$HuH6$PAM.Cluster <- factor(as.character(clusters.df[huh6.set]))
birc5.seurat$PAM.Cluster <- factor(as.character(clusters.df[all.set]))
saveRDS(birc5.list, "data/birc5_seurat_split_list_pam.rds")
saveRDS(birc5.seurat, "data/birc5_seurat_hvgs_pam.rds")

## [ Meta-clustering ] ----

## Meta-clustering using Jaccard similarity and PAM clustering had limited success
## in finding similarity between clusters in different cell lines
## Try manual meta-clustering by comparing the expression similarity between clusters in different samples

birc5.list <- readRDS("data/birc5_seurat_split_list_pam.rds")

## Assign inital clusters from separated objects into combined
birc5.list$HB303$Init.Cluster <- paste0("HB303_", birc5.list$HB303$seurat_clusters.0.4)
birc5.list$HepG2$Init.Cluster <- paste0("HepG2_", birc5.list$HepG2$seurat_clusters.0.4)
birc5.list$HuH6$Init.Cluster <- paste0("HuH6_", birc5.list$HuH6$seurat_clusters.0.4)

birc5.seurat <- Reduce(function(x, y) merge(x, y), birc5.list)

clusters <- unique(birc5.seurat$Init.Cluster)
order <- clusters[order(match(sub("_.*", "", clusters), c("HB303", "HepG2", "HuH6")), as.integer(sub(".*_", "", clusters)))]
birc5.seurat$Init.Cluster <- factor(birc5.seurat$Init.Cluster, levels = order)
saveRDS(birc5.seurat, "data/birc5_seurat_hvgs_init.rds")

birc5.list$HB303$Init.Cluster <- factor(birc5.list$HB303$Init.Cluster, levels = order[order %in% unique(birc5.list$HB303$Init.Cluster)])
birc5.list$HepG2$Init.Cluster <- factor(birc5.list$HepG2$Init.Cluster, levels = order[order %in% unique(birc5.list$HepG2$Init.Cluster)])
birc5.list$HuH6$Init.Cluster <- factor(birc5.list$HuH6$Init.Cluster, levels = order[order %in% unique(birc5.list$HuH6$Init.Cluster)])
saveRDS(birc5.list, "data/birc5_seurat_split_list_init.rds")

init.cols <- hues::iwanthue(length(unique(birc5.seurat$Init.Cluster)))
names(init.cols) <- order
saveRDS(init.cols, "data/initial_cluster_cols.rds")

for (cell in names(birc5.list)) {
  DimPlot(birc5.list[[cell]], group.by = "Init.Cluster", order = TRUE) +
    scale_colour_manual(values = init.cols) +
    umap.theme() + labs(title = paste(cell, "Initial Clusters"))
  ggsave(paste0("plots/umaps/umap_init_", tolower(cell),".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_init_", tolower(cell),".png"), width = 8.3, height = 5.8)
}

init.df <- birc5.seurat@meta.data %>%
  count(Condition, Init.Cluster) %>%
  group_by(Condition) %>%                 
  mutate(prop = n / sum(n))  
ggplot(init.df, aes(x = Condition, y = prop, fill = Init.Cluster)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = init.cols) +
  guides(fill = guide_legend(title = "Initial clusters", ncol = 2)) +
  theme_bw() +
  theme(text = element_text(size = 15), aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/anno/proportion_init.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/proportion_init.png", width = 8.3, height = 5.8)

init.df <- init.df %>%
  mutate(Cell_Line = sub("_.*$", "", Init.Cluster)) %>%
  group_by(Cell_Line, Condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
for (cell in unique(init.df$Cell_Line)) {
  ggplot(filter(init.df, Cell_Line == cell),
         aes(x = Condition, y = prop, fill = Init.Cluster)) +
    geom_bar(stat = "identity") +
    ylab("Proportion") +
    xlab("") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = init.cols) +
    guides(fill = guide_legend(title = "Initial clusters")) +
    ggtitle(cell) +
    theme_bw() +
    theme(text = element_text(size = 15), aspect.ratio = 0.75,
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0("plots/anno/proportion_init_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/anno/proportion_init_", tolower(cell), ".png"), width = 8.3, height = 5.8)
}

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk //////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

pseudo.init <- AggregateExpression(birc5.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Init.Cluster"))
saveRDS(pseudo.init, "data/birc5_pseudobulk_init.rds")

## MSigDB Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", collection = "H")
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
init.levels <- gsub("_", "-", levels(birc5.seurat$Init.Cluster))
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
ggsave("plots/anno/pseudo_init_msigdb_H_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_msigdb_H_pheatmap.png", plot = gt, width = 8.3, height = 8.3)

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
ggsave("plots/anno/pseudo_init_msigdb_C2_liver_heatmap.pdf", width = 8.3, height = 11.7)
ggsave("plots/anno/pseudo_init_msigdb_C2_liver_heatmap.png", width = 8.3, height = 11.7)
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
ggsave("plots/anno/pseudo_init_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 8.3, height = 8.3)
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
ggsave("plots/anno/pseudo_init_msigdb_C2_cairo_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_msigdb_C2_cairo_pheatmap.png", plot = gt, width = 8.3, height = 8.3)

## HB signatures from literature
sig.df <- as.data.frame(readxl::read_xlsx(path = "../huh6/data/hb_sigs.xlsx", col_names = FALSE))
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
ggsave("plots/anno/pseudo_init_hb_sigs_heatmap.pdf", width = 8.3, height = 11.7)
ggsave("plots/anno/pseudo_init_hb_sigs_heatmap.png", width = 8.3, height = 11.7)

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
ggsave("plots/anno/pseudo_init_hb_sigs_pheatmap.pdf", plot = gt, width = 5.8, height = 8.3)
ggsave("plots/anno/pseudo_init_hb_sigs_pheatmap.png", plot = gt, width = 5.8, height = 8.3)

saveRDS(pseudo.init, "data/birc5_pseudobulk_init_sigs.rds")

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
ggsave("plots/anno/pseudo_init_custom_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_init_custom_pheatmap.png", plot = gt, width = 8.3, height = 8.3)

## Heatmaps by cell line
sample.id <- sub("-.*$", "", colnames(init.custom.mat))
for (cell in unique(sample.id)) {
  cols <- which(sample.id == cell)
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
                 main = cell)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_init_custom_pheatmap_", tolower(cell), ".pdf"), plot = gt, width = 5.8, height = 8.3)
  ggsave(filename = paste0("plots/anno/pseudo_init_custom_pheatmap_", tolower(cell), ".png"), plot = gt, width = 5.8, height = 8.3)
}

## /////////////////////////////////////////////////////////////////////////////
## Jaccard /////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

hb303.markers <- readRDS("data/jaccard/hb303_markers.rds")
hepg2.markers <- readRDS("data/jaccard/hepg2_markers.rds")
huh6.markers <- readRDS("data/jaccard/huh6_markers.rds")
ref.markers <- readRDS("../huh6/data/jaccard/huh6_pam_markers.rds")
ref.markers$name <- recode(as.character(ref.markers$cluster),
                           "1" = "Hepatocytic 2",
                           "2" = "Stem-like",
                           "3" = "Hepatocytic 3",
                           "4" = "Hepatocytic 1",
                           "5" = "Late progenitor",
                           "7" = "Early progenitor")

## Filter significant and higher expressed markers
hb303.markers <- hb303.markers[hb303.markers$p_val_adj < 0.05, ]
hepg2.markers <- hepg2.markers[hepg2.markers$p_val_adj < 0.05, ]
huh6.markers <- huh6.markers[huh6.markers$p_val_adj < 0.05, ]
ref.markers <- ref.markers[ref.markers$p_val_adj < 0.05, ]

hb303.markers <- hb303.markers[hb303.markers$avg_log2FC > 0.5, ]
hepg2.markers <- hepg2.markers[hepg2.markers$avg_log2FC > 0.5, ]
huh6.markers <- huh6.markers[huh6.markers$avg_log2FC > 0.5, ]
ref.markers <- ref.markers[ref.markers$avg_log2FC > 0.5, ]
## Split into lists of genes per cluster
hb303.markers <- split(hb303.markers, hb303.markers$cluster)
hepg2.markers <- split(hepg2.markers, hepg2.markers$cluster)
huh6.markers <- split(huh6.markers, huh6.markers$cluster)
ref.markers <- split(ref.markers, ref.markers$name)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(hb303.markers, nrow) ## 0 = 35, 1 = 743, 2 = 114, 3 = 245, 4 = 40, 5 = 331
lapply(hepg2.markers, nrow) ## 0 = 73, 1 = 267, 2 = 313, 3 = 132, 4 = 362
lapply(huh6.markers, nrow) ## 0 = 338, 1 = 187, 2 = 422, 3 = 149, 4 = 285, 5 = 276, 6 = 1335, 7 = 1284, 8 = 499
lapply(ref.markers, nrow) ## Early progenitor = 869, Hepatocytic 1 = 1121, Hepatocytic 2 = 441, Hepatocytic 3 = 548, = Late progenitor = 1348, = Stem-like = 199

hb303.markers <- lapply(hb303.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
hepg2.markers <- lapply(hepg2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
huh6.markers <- lapply(huh6.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
ref.markers <- lapply(ref.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})

## Give your clusters unique names
names(hb303.markers) <- paste0("HB303_", names(hb303.markers))
names(hepg2.markers) <- paste0("HepG2_", names(hepg2.markers))
names(huh6.markers) <- paste0("HuH6_", names(huh6.markers))
names(ref.markers) <- paste0("Ref_", names(ref.markers))
## Make a big list of all the cluster markers
markers.list <- c(hb303.markers, hepg2.markers, huh6.markers, ref.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
jaccard.mat <- matrix(data = NA, nrow = length(markers.list),
                      ncol = length(markers.list),
                      dimnames = list(names(markers.list),
                                      names(markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(jaccard.mat)) {
  for (j in colnames(jaccard.mat)) {
    jaccard.mat[i,j] <- jaccard(markers.list[[i]]$gene, markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Sample" = gsub("_.+", "", colnames(jaccard.mat)),
                      row.names = colnames(jaccard.mat))
anno.col <- list("Sample" = c("Ref" = "black",
                              "HB303" = as.character(cell.cols[1]),
                              "HepG2" = as.character(cell.cols[2]),
                              "HuH6" = as.character(cell.cols[3])))
dev.off()
gt <- pheatmap(jaccard.mat,
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
ggsave("plots/anno/ref_init_jaccard_heatmap.png", plot = gt, height = 8.3, width = 8.3)
ggsave("plots/anno/ref_init_jaccard_heatmap.pdf", plot = gt, height = 8.3, width = 8.3)

## /////////////////////////////////////////////////////////////////////////////
## Final annotation ////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

birc5.list$HB303$Final.Cluster <- recode(as.character(birc5.list$HB303$Init.Cluster),
                                         "HB303_0" = "Progenitor",
                                         "HB303_1" = "Hepatocytic",
                                         "HB303_2" = "Hepatocytic",
                                         "HB303_3" = "Progenitor",
                                         "HB303_4" = "Hepatocytic",
                                         "HB303_5" = "Progenitor")
birc5.list$HepG2$Final.Cluster <- recode(as.character(birc5.list$HepG2$Init.Cluster),
                                         "HepG2_0" = "Hepatocytic",
                                         "HepG2_1" = "Hepatocytic",
                                         "HepG2_2" = "Progenitor 2",
                                         "HepG2_3" = "Progenitor 1",
                                         "HepG2_4" = "Progenitor 1")
birc5.list$HuH6$Final.Cluster <- recode(as.character(birc5.list$HuH6$Init.Cluster),
                                        "HuH6_0" = "Progenitor 1",
                                        "HuH6_1" = "Progenitor 1",
                                        "HuH6_2" = "Progenitor 1",
                                        "HuH6_3" = "Progenitor 2",
                                        "HuH6_4" = "Hepatocytic 1",
                                        "HuH6_5" = "Hepatocytic 1",
                                        "HuH6_6" = "Hepatocytic 2",
                                        "HuH6_7" = "Hepatocytic 2",
                                        "HuH6_8" = "Progenitor 2")
saveRDS(birc5.list, "data/birc5_seurat_split_list_final.rds")

## Similarity with HuH6 cisplatin scRNA-seq clusters
ref.seurat <- readRDS("../huh6/data/huh6_seurat_hvgs_pam_anno.rds")
ref.seurat$Final.Cluster <- ref.seurat$PAM_Name
seurat.list <- c(birc5.list, ref.seurat)
names(seurat.list)[4] <- "Ref"
pseudo.list <- list()
for (sample in names(seurat.list)) {
  obj <- seurat.list[[sample]]
  pseudo <- AggregateExpression(obj, assays = "RNA", return.seurat = T,
                                group.by = c("Final.Cluster"))
  pseudo.list[[sample]] <- pseudo
}
## Combine pseudobulk objects
mat.list <- lapply(pseudo.list, function(obj) {
  GetAssayData(obj, assay = "RNA", slot = "data")
})
common.genes <- Reduce(intersect, lapply(mat.list, rownames))
mat.list <- lapply(mat.list, function(m) m[common.genes, , drop = FALSE])
mat.list <- Map(function(m, nm) {
  colnames(m) <- paste0(nm, "_", colnames(m))
  m
}, mat.list, names(mat.list))
mat.all <- do.call(cbind, mat.list)
## Combine annotations
anno.list <- Map(function(obj, nm) {
  md <- obj[[]]
  md$Cell_Line <- nm
  md$Cluster <- paste0(nm, "_", rownames(md))
  md
}, pseudo.list, names(pseudo.list))
anno.all <- bind_rows(anno.list)
rownames(anno.all) <- anno.all$Cluster
anno.all <- anno.all[colnames(mat.all), , drop = FALSE]
anno.df <- anno.all[, "Cell_Line", drop = FALSE]
anno.col <- list("Cell_Line" = c("Ref" = "black",
                                 "HB303" = as.character(cell.cols[1]),
                                 "HepG2" = as.character(cell.cols[2]),
                                 "HuH6" = as.character(cell.cols[3])))
## Add z-score per gene
mat.z <- t(scale(t(as.matrix(mat.all))))
mat.z[is.na(mat.z)] <- 0

custom.mat <- mat.z[custom, , drop = FALSE]
gt <- pheatmap(custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = anno.col)$gtable
ggsave("plots/anno/pseudo_ref_custom_pheatmap.pdf", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/anno/pseudo_ref_custom_pheatmap.png", plot = gt, width = 8.3, height = 8.3)

## Heatmaps by cell line
sample.id <- sub("_.*$", "", colnames(custom.mat))
for (cell in unique(sample.id)) {
  cols <- which(sample.id == cell)
  mat.sub  <- custom.mat[, cols, drop = FALSE]
  gt <- pheatmap(mat.sub,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 main = cell)$gtable
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", tolower(cell), ".pdf"), plot = gt, width = 5.8, height = 8.3)
  ggsave(filename = paste0("plots/anno/pseudo_final_custom_pheatmap_", tolower(cell), ".png"), plot = gt, width = 5.8, height = 8.3)
}

## Find final cluster markers
final.markers <- list()
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
for (cell in names(birc5.list)) {
  obj <- birc5.list[[cell]]
  Idents(obj) <- obj$Final.Cluster
  markers <- FindAllMarkers(obj, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
  final.markers[[cell]] <- markers
}
saveRDS(final.markers, "data/final_markers/final_markers_list.rds")

## Jaccard similarity
hb303.markers <- final.markers$HB303
hepg2.markers <- final.markers$HepG2
huh6.markers <- final.markers$HuH6
ref.markers <- readRDS("../huh6/data/jaccard/huh6_pam_markers.rds")
ref.markers$name <- recode(as.character(ref.markers$cluster),
                           "1" = "Hepatocytic 2",
                           "2" = "Stem-like",
                           "3" = "Hepatocytic 3",
                           "4" = "Hepatocytic 1",
                           "5" = "Late progenitor",
                           "7" = "Early progenitor")

hb303.markers <- hb303.markers[hb303.markers$p_val_adj < 0.05, ]
hepg2.markers <- hepg2.markers[hepg2.markers$p_val_adj < 0.05, ]
huh6.markers <- huh6.markers[huh6.markers$p_val_adj < 0.05, ]
ref.markers <- ref.markers[ref.markers$p_val_adj < 0.05, ]

hb303.markers <- hb303.markers[hb303.markers$avg_log2FC > 0.5, ]
hepg2.markers <- hepg2.markers[hepg2.markers$avg_log2FC > 0.5, ]
huh6.markers <- huh6.markers[huh6.markers$avg_log2FC > 0.5, ]
ref.markers <- ref.markers[ref.markers$avg_log2FC > 0.5, ]

hb303.markers <- split(hb303.markers, hb303.markers$cluster)
hepg2.markers <- split(hepg2.markers, hepg2.markers$cluster)
huh6.markers <- split(huh6.markers, huh6.markers$cluster)
ref.markers <- split(ref.markers, ref.markers$name)

lapply(hb303.markers, nrow) ## Hepatocytic = 553, Progenitor = 115
lapply(hepg2.markers, nrow) ## Hepatocytic = 174, Progenitor 2 = 129, Progenitor 1 = 313
lapply(huh6.markers, nrow) ## Hepatocytic 1 = 272, Hepatocytic 2 = 1356, Progenitor 1 = 785, Progenitor 2 = 154
lapply(ref.markers, nrow) ## Early progenitor = 869, Hepatocytic 1 = 1121, Hepatocytic 2 = 441, Hepatocytic 3 = 548, = Late progenitor = 1348, = Stem-like = 199

hb303.markers <- lapply(hb303.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
hepg2.markers <- lapply(hepg2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
huh6.markers <- lapply(huh6.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
ref.markers <- lapply(ref.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})

names(hb303.markers) <- paste0("HB303_", names(hb303.markers))
names(hepg2.markers) <- paste0("HepG2_", names(hepg2.markers))
names(huh6.markers) <- paste0("HuH6_", names(huh6.markers))
names(ref.markers) <- paste0("Ref_", names(ref.markers))
markers.list <- c(hb303.markers, hepg2.markers, huh6.markers, ref.markers)

jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
jaccard.mat <- matrix(data = NA, nrow = length(markers.list),
                      ncol = length(markers.list),
                      dimnames = list(names(markers.list),
                                      names(markers.list)))
for (i in rownames(jaccard.mat)) {
  for (j in colnames(jaccard.mat)) {
    jaccard.mat[i,j] <- jaccard(markers.list[[i]]$gene, markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Sample" = gsub("_.+", "", colnames(jaccard.mat)),
                      row.names = colnames(jaccard.mat))
anno.col <- list("Sample" = c("Ref" = "black",
                              "HB303" = as.character(cell.cols[1]),
                              "HepG2" = as.character(cell.cols[2]),
                              "HuH6" = as.character(cell.cols[3])))
dev.off()
gt <- pheatmap(jaccard.mat,
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
ggsave("plots/anno/ref_final_jaccard_heatmap.png", plot = gt, height = 5.8, width = 8.3)
ggsave("plots/anno/ref_final_jaccard_heatmap.pdf", plot = gt, height = 5.8, width = 8.3)

## [ Cell line UMAPs ] ----

## Split data
birc5.list <- readRDS("data/birc5_seurat_split_list_final.rds")
birc5.list <- lapply(birc5.list, FindVariableFeatures)
birc5.seurat <- readRDS("data/birc5_seurat_hvgs_init.rds")
rm.genes <- c(grep("^MT-", rownames(birc5.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(birc5.seurat), value = TRUE),
              grep("[\\.]", rownames(birc5.seurat), value = TRUE),
              grep("^LINC", rownames(birc5.seurat), value = TRUE),
              c("MALAT1"))
output <- list()
for (i in 1:length(birc5.list)) {
  obj <- birc5.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(birc5.list)
rm(birc5.list, birc5.seurat, obj)

output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 100,
         features = VariableFeatures(obj))
})
output <- lapply(output, RunUMAP, reduction = "pca", dims = 1:30)
output <- lapply(output, FindNeighbors, reduction = "pca", dims = 1:30)
saveRDS(output, "data/birc5_seurat_split_list_hvgs.rds")

birc5.list <- readRDS("data/birc5_seurat_split_list_hvgs.rds")
for (cell in names(birc5.list)) {
  DimPlot(birc5.list[[cell]], group.by = "Condition", order = TRUE) +
    scale_colour_manual(values = group.cols) +
    umap.theme() + labs(title = paste(cell, " - Conditions"))
  ggsave(paste0("plots/umaps/umap_condition_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_condition_", tolower(cell), ".png"), width = 8.3, height = 5.8)
}

hb303.cols <- c("Progenitor" = "#E37A2B",
                "Hepatocytic" = "#ED7474")
hepg2.cols <- c("Progenitor 1" = "#E37A2B",
                "Progenitor 2" = "#EDB074",
                "Hepatocytic" = "#ED7474")
huh6.cols <- c("Hepatocytic 1" = "#ED7474",
               "Hepatocytic 2" = "#E4A1A1",
               "Progenitor 1" = "#E37A2B",
               "Progenitor 2" = "#EDB074")
final.cols <- c(hb303.cols, hepg2.cols, huh6.cols)
final.cols <- final.cols[!duplicated(names(final.cols))]
saveRDS(final.cols, "data/final_cluster_cols.rds")

all.clusters <-c("Progenitor", "Progenitor 1", "Progenitor 2", "Progenitor 3",
                 "Hepatocytic", "Hepatocytic 1", "Hepatocytic 2")
for (cell in names(birc5.list)) {
  seurat <- birc5.list[[cell]]
  order <- all.clusters[all.clusters %in% unique(seurat$Final.Cluster)]
  seurat$Final.Cluster <- factor(seurat$Final.Cluster, levels = order)
  DimPlot(seurat, group.by = "Final.Cluster", order = TRUE) +
    scale_colour_manual(values = final.cols) +
    umap.theme() + labs(title = paste(cell, " - Final Cluster"))
  ggsave(paste0("plots/umaps/umap_final_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/umaps/umap_final_", tolower(cell), ".png"), width = 8.3, height = 5.8)
}

meta.combined <- do.call(rbind, lapply(seq_along(birc5.list), function(i) {
  df <- birc5.list[[i]]@meta.data[, c("Condition", "Final.Cluster", "Cell_Line")]
  return(df)}))
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
ggsave("plots/anno/final_cluster_bar.pdf", width = 8.3, height = 5.8)
ggsave("plots/anno/final_cluster_bar.png", width = 8.3, height = 5.8)

cluster.filt <- meta.combined %>%
  group_by(Cell_Line, Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
cluster.filt$Cell_Line <- sub("_[^_]+$", "", cluster.filt$Cell_Line)
cell.lines <- unique(cluster.filt$Cell_Line)
for (cell in cell.lines) {
  df <- cluster.filt[cluster.filt$Cell_Line == cell, ]
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
  ggsave(paste0("plots/anno/", tolower(cell),"_final_cluster_bar.pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/anno/", tolower(cell),"_final_cluster_bar.png"), width = 5.8, height = 5.8)
}

## BIRC5 and E2F4 plots
for (cell in names(birc5.list)) {
  seurat <- birc5.list[[cell]]
  Idents(seurat) <- "Condition"
  VlnPlot(seurat, features = "BIRC5", pt.size = 0, cols = group.cols) +
    xlab("") &
    labs(title = paste(cell, " - BIRC5 Expression")) &
    theme(legend.position = "none",
          text = element_text(size = 12))
  ggsave(paste0("plots/birc5/birc5_condition_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/birc5/birc5_condition_", tolower(cell), ".png"), width = 8.3, height = 5.8)
  VlnPlot(seurat, features = "E2F4", pt.size = 0, cols = group.cols) +
    xlab("") &
    labs(title = paste(cell, " - E2F4 Expression")) &
    theme(legend.position = "none",
          text = element_text(size = 12))
  ggsave(paste0("plots/birc5/e2f4_condition_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/birc5/e2f4_condition_", tolower(cell), ".png"), width = 8.3, height = 5.8)
}

for (cell in names(birc5.list)) {
  seurat <- birc5.list[[cell]]
  order <- all.clusters[all.clusters %in% unique(seurat$Final.Cluster)]
  seurat$Final.Cluster <- factor(seurat$Final.Cluster, levels = order)
  Idents(seurat) <- "Final.Cluster"
  VlnPlot(seurat, features = "BIRC5", pt.size = 0, cols = final.cols) +
    xlab("") &
    labs(title = paste(cell, " - BIRC5 Expression")) &
    theme(legend.position = "none",
          text = element_text(size = 12))
  ggsave(paste0("plots/birc5/birc5_final_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/birc5/birc5_final_", tolower(cell), ".png"), width = 8.3, height = 5.8)
  VlnPlot(seurat, features = "E2F4", pt.size = 0, cols = final.cols) +
    xlab("") &
    labs(title = paste(cell, " - E2F4 Expression")) &
    theme(legend.position = "none",
          text = element_text(size = 12))
  ggsave(paste0("plots/birc5/e2f4_final_", tolower(cell), ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/birc5/e2f4_final_", tolower(cell), ".png"), width = 8.3, height = 5.8)
}

## [ Figure plots ] ----

## Remove E2F4 siRNA samples from UMAPs and bar charts
birc5.list <- readRDS("data/birc5_seurat_split_list_final.rds")
birc5.list <- lapply(birc5.list, function(obj) {
  subset(obj, subset = Condition %in% c("NTC", "BIRC5_siRNA"))
})
birc5.list <- lapply(birc5.list, FindVariableFeatures)
birc5.seurat <- readRDS("data/birc5_seurat_hvgs_init.rds")
rm.genes <- c(grep("^MT-", rownames(birc5.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(birc5.seurat), value = TRUE),
              grep("[\\.]", rownames(birc5.seurat), value = TRUE),
              grep("^LINC", rownames(birc5.seurat), value = TRUE),
              c("MALAT1"))
output <- list()
for (i in 1:length(birc5.list)) {
  obj <- birc5.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(birc5.list)
rm(birc5.list, birc5.seurat, obj)

output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 100,
         features = VariableFeatures(obj))
})
output <- lapply(output, RunUMAP, reduction = "pca", dims = 1:30)
output <- lapply(output, FindNeighbors, reduction = "pca", dims = 1:30)
saveRDS(output, "data/birc5_seurat_split_list_filt.rds")

birc5.list <- readRDS("data/birc5_seurat_split_list_filt.rds")
for (cell in names(birc5.list)) {
  DimPlot(birc5.list[[cell]], group.by = "Condition", order = TRUE) +
    scale_colour_manual(values = group.cols) +
    umap.theme() + labs(title = paste(cell, " - Conditions")) +
    theme(text = element_text(size = 20))
  ggsave(paste0("plots/figures/umap_condition_", tolower(cell), ".pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/figures/umap_condition_", tolower(cell), ".png"), width = 5.8, height = 5.8)
}

final.cols <- readRDS("data/final_cluster_cols.rds")
all.clusters <-c("Progenitor", "Progenitor 1", "Progenitor 2", "Progenitor 3",
                 "Hepatocytic", "Hepatocytic 1", "Hepatocytic 2")
for (cell in names(birc5.list)) {
  seurat <- birc5.list[[cell]]
  order <- all.clusters[all.clusters %in% unique(seurat$Final.Cluster)]
  seurat$Final.Cluster <- factor(seurat$Final.Cluster, levels = order)
  DimPlot(seurat, group.by = "Final.Cluster", order = TRUE) +
    scale_colour_manual(values = final.cols) +
    umap.theme() + labs(title = paste(cell, " - Final Cluster")) +
    theme(text = element_text(size = 20))
  ggsave(paste0("plots/figures/umap_final_", tolower(cell), ".pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/figures/umap_final_", tolower(cell), ".png"), width = 5.8, height = 5.8)
}

meta.combined <- do.call(rbind, lapply(seq_along(birc5.list), function(i) {
  df <- birc5.list[[i]]@meta.data[, c("Condition", "Final.Cluster", "Cell_Line")]
  return(df)}))
cluster.df <- meta.combined %>%
  group_by(Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
cluster.df$Condition <- factor(cluster.df$Condition, levels = c("NTC", "BIRC5_siRNA"))
ggplot(cluster.df, aes(x = Condition, y = prop, fill = Final.Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = final.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 1.5)
ggsave("plots/figures/final_cluster_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/figures/final_cluster_bar.png", width = 5.8, height = 5.8)

cluster.filt <- meta.combined %>%
  group_by(Cell_Line, Condition, Final.Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(prop = n / sum(n))
cluster.filt$Condition <- factor(cluster.filt$Condition, levels = c("NTC", "BIRC5_siRNA"))
cluster.filt$Cell_Line <- sub("_[^_]+$", "", cluster.filt$Cell_Line)
cell.lines <- unique(cluster.filt$Cell_Line)
for (cell in cell.lines) {
  df <- cluster.filt[cluster.filt$Cell_Line == cell, ]
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
  ggsave(paste0("plots/figures/", tolower(cell),"_final_cluster_bar.pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/figures/", tolower(cell),"_final_cluster_bar.png"), width = 5.8, height = 5.8)
}

## Remove some genes from annotation heatmaps
birc5.list <- readRDS("data/birc5_seurat_split_list_final.rds")
custom <- c("EPCAM", "NCAM1", "CLDN3", "PROM1", "SOX17", ## HSC markers
            "KRT8", "KRT18", "KRT19", "CXCR4", ## HSC markers
            "AFP", "KRT7", "ICAM1", "ONECUT2", "TBX3",  ## HB markers
            "HNF4A", "CEBPA", "ALB", "RBPJ", "NR5A2", ## Hepatocytic cell fate
            "MAT1A", "NR1I2", "MAT2A", ## Adult liver; MAT1A = fetal liver + replaces MAT1A in HCC
            "SOX9", "HNF1B", "SALL4", ## Cholangiocyte fate regulators
            "TGFB1", "TGFB2", "TGFBR2") ## High TGFb = cholangiocyte differentiation

## Similarity with HuH6 cisplatin scRNA-seq clusters
ref.seurat <- readRDS("../huh6/data/huh6_seurat_hvgs_pam_anno.rds")
ref.seurat$Final.Cluster <- ref.seurat$PAM_Name
seurat.list <- c(birc5.list, ref.seurat)
names(seurat.list)[4] <- "Ref"
pseudo.list <- list()
for (sample in names(seurat.list)) {
  obj <- seurat.list[[sample]]
  pseudo <- AggregateExpression(obj, assays = "RNA", return.seurat = T,
                                group.by = c("Final.Cluster"))
  pseudo.list[[sample]] <- pseudo
}
## Combine pseudobulk objects
mat.list <- lapply(pseudo.list, function(obj) {
  GetAssayData(obj, assay = "RNA", slot = "data")
})
common.genes <- Reduce(intersect, lapply(mat.list, rownames))
mat.list <- lapply(mat.list, function(m) m[common.genes, , drop = FALSE])
mat.list <- Map(function(m, nm) {
  colnames(m) <- paste0(nm, "_", colnames(m))
  m
}, mat.list, names(mat.list))
mat.all <- do.call(cbind, mat.list)
## Combine annotations
anno.list <- Map(function(obj, nm) {
  md <- obj[[]]
  md$Cell_Line <- nm
  md$Cluster <- paste0(nm, "_", rownames(md))
  md
}, pseudo.list, names(pseudo.list))
anno.all <- bind_rows(anno.list)
rownames(anno.all) <- anno.all$Cluster
anno.all <- anno.all[colnames(mat.all), , drop = FALSE]
anno.df <- anno.all[, "Cell_Line", drop = FALSE]
anno.col <- list("Cell_Line" = c("Ref" = "black",
                                 "HB303" = as.character(cell.cols[1]),
                                 "HepG2" = as.character(cell.cols[2]),
                                 "HuH6" = as.character(cell.cols[3])))
## Add z-score per gene
mat.z <- t(scale(t(as.matrix(mat.all))))
mat.z[is.na(mat.z)] <- 0

custom.mat <- mat.z[custom, , drop = FALSE]
gt <- pheatmap(custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = anno.col)$gtable
ggsave("plots/figures/pseudo_ref_custom_pheatmap.pdf", plot = gt)
ggsave("plots/figures/pseudo_ref_custom_pheatmap.png", plot = gt)

## Heatmaps by cell line
final.cols <- readRDS("data/final_cluster_cols.rds")
sample.id <- sub("_.*$", "", colnames(custom.mat))
lines <- c("HB303", "HepG2", "HuH6")
for (cell in lines) {
  cols <- which(sample.id == cell)
  mat.sub  <- custom.mat[, cols, drop = FALSE]
  anno.df <- data.frame(row.names = colnames(mat.sub), "Cluster" = sub("^.*?_", "", colnames(mat.sub)))
  anno.col <- list("Cluster" = final.cols[names(final.cols) %in% anno.df$Cluster])
  gt <- pheatmap(mat.sub,
                 border_color = NA,
                 cellwidth = 8, cellheight = 8,
                 fontsize_row = 8, fontsize_col = 8,
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "ward.D2",
                 color = Seurat:::SpatialColors(100),
                 main = cell,
                 annotation_col = anno.df,
                 annotation_colors = anno.col)$gtable
  ggsave(filename = paste0("plots/figures/pseudo_final_custom_pheatmap_", tolower(cell), ".pdf"), plot = gt)
  ggsave(filename = paste0("plots/figures/pseudo_final_custom_pheatmap_", tolower(cell), ".png"), plot = gt)
}
