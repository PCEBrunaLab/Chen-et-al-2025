## HuH6 scRNA-seq analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(scater)
library(scran)
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(TSCAN)
library(slingshot)
library(biomaRt)
library(harmony)
library(tricycle)
library(cluster)
library(clusterProfiler)
library(msigdbr)
library(RcisTarget)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ggbeeswarm)
library(ggsankey)
library(pheatmap)
library(reshape2)
library(patchwork)
library(viridis)
library(hues)
library(writexl)
library(readxl)

## Labels for experimental conditions
conditions.list <- c("POT", "Untreated", "Cisplatin", bquote("Recovery t"[1]), bquote("Recovery t"[2]), bquote("Recovery t"[3]))

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
group.cols <- c("POT" = "#228833", 
                "Untreated" = "#4477AA", 
                "Cisplatin_off" = "#CCBB44", 
                "Cisplatin_recovery1" = "#EE6677", 
                "Cisplatin_recovery2" = "#AA3377", 
                "Cisplatin_recovery3" = "#66CCEE")

cool.warm.pal <- colorRampPalette(c("#1d4877", "#1b8a5a", "#fbb021", "#f68838", "#ee3e32"))

## Colours for PAM clusters
require(paletteer)
pal <- paletteer_d("rcartocolor::Pastel")
pam.cols <- c("Hepatocytic 1" = pal[2], 
              "Hepatocytic 2" = pal[4], 
              "Hepatocytic 3" = pal[3], 
              "Late progenitor" = pal[5],
              "Early progenitor" = pal[1],
              "Stem-like" = pal[6])

## [ Data preparation ] ----

huh6.sce <- readRDS("data/HuH6_SCE-norm.RDS")

## Demux samples
demux.list <- list.files("data/data prep/", pattern = "demuxed", full.names = TRUE)
d.list <- list()
for(x in seq_along(demux.list)){
  x.d <- read.table(demux.list[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x.name <- gsub(demux.list[x], pattern = "(.*/)(Multiplex_[0-9]*)(_demuxed\\.tsv)", replacement = "\\2")
  x.d$Sample <- x.name
  d.list[[x.name]] <- x.d
}

## Make demux data frame
demux.df <- do.call(rbind.data.frame, d.list)
demux.df$CellID <- paste(demux.df$Sample, demux.df$ID, sep = "_")

## Remove path format from metadata values
col.df <- as.data.frame(colData(huh6.sce))
col.df$Sample <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", col.df$Sample)
col.df$Sample <- gsub(pattern = "(/)", "", col.df$Sample)
col.df$CellID <- paste(col.df$Sample, col.df$Barcode, sep="_")

unique(col.df$Sample)

## Demux sample conditions
config.df <- read.table("data/data prep/HuH6_config.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
## Recovery column
config.df$Rec <- FALSE
rec.patterns <- c("off", "recovery")
config.df$Rec[grepl(config.df$description, pattern = paste(rec.patterns, collapse = "|"))] <- TRUE
## Condition column
config.df$Condition <- gsub(config.df$description, pattern = "(\\S+)(_[A-C])", replacement="\\1")
## Fix sample ID mistake
config.df$sample_id <- gsub(config.df$sample_id, pattern="HuH_sample9", replacement="HuH6_sample9")

## Get size factors for QC
sf.df <- data.frame("CellID" = colnames(huh6.sce), "Sample" = colData(huh6.sce)$Sample, "sizeFactor" = sizeFactors(huh6.sce))
sf.df$Sample <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", sf.df$Sample)
sf.df$Sample <- gsub(pattern = "(/)", "", sf.df$Sample)

meta.df <- Reduce(x = list(col.df, demux.df, sf.df), f = function(x, y) merge(x, y, by = c("CellID", "Sample")))
meta.df <- merge(meta.df, config.df, by.x = c("Sample", "CMO"), by.y = c("Sample", "cmo_ids"))

## CMO and cell counts metadata
cmo.sum.files <- list.files("data/data prep/", pattern = "CMOsummary", full.names = TRUE)
cell.sum.files <- list.files("data/data prep/", pattern = "Cellsummary", full.names = TRUE)
cmo.list <- list()
cell.list <- list()
for(x in seq_along(cmo.sum.files)){
  x.cmo <- read.table(cmo.sum.files[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x.cell <- read.table(cell.sum.files[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cmo.list[[x]] <- x.cmo
  cell.list[[x]] <- x.cell
}

cmo.sum.df <- do.call(rbind.data.frame, cmo.list)
cell.sum.df <- do.call(rbind.data.frame, cell.list)

cell.sum.df$CellID <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", cell.sum.df$CellID)
cell.sum.df$CellID <- gsub(pattern = "(/)", "", cell.sum.df$CellID)

cell.cmo.df <- merge(demux.df, cell.sum.df, by = c("CellID", "Sample"))

intersect(colnames(cell.cmo.df), colnames(meta.df))
test.df <- merge(cell.cmo.df, meta.df, by = c("CellID", "Sample", "ID", "Class", "CMO"))
test.df <- test.df %>%
  dplyr::select(-c(sizeFactor.y)) %>%
  mutate(sizeFactor = sizeFactor.x) %>%
  dplyr::select(-c(sizeFactor.x))
sum(is.na(test.df$CellID)) ## 0

saveRDS(test.df, "data/HuH6_metadata.rds")

## Add metadata to SCE object
final.meta <- readRDS("data/HuH6_metadata.rds")
colData(huh6.sce)$CellID <- rownames(colData(huh6.sce))
huh6.sce$Sample <- NULL
## Remove cells not in metadata
huh6.sce$CellID[!(huh6.sce$CellID %in% final.meta$CellID)]
## "Multiplex_5_AATCACGAGATGCTGG-1"
## "Multiplex_5_ATAGGCTCATTCTGTT-1"
## "Multiplex_5_CCAAGCGGTTCGAACT-1"
## "Multiplex_5_CGGCAGTCATGGGCAA-1"
## "Multiplex_5_GCCAGGTCAAGGTCGA-1"
## These 5 cells have CMO11, which was not used for the HuH6

cell.nms <- colnames(huh6.sce)[colnames(huh6.sce) %in% final.meta$CellID]
huh6.sce <- huh6.sce[ , final.meta$CellID]
dim(huh6.sce) ## 36613 14001
colData(huh6.sce) <- merge(colData(huh6.sce), final.meta, by = c("CellID", "Barcode", "sizeFactor"))
colnames(huh6.sce) <- cell.nms 

saveRDS(huh6.sce, "data/HuH6_sce_meta.rds")

## [ Initial filtering ] ----

huh6.sce <- readRDS("data/HuH6_sce_meta.rds")

#ensembl <- useMart("ensembl", "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
#ens.bm <- getBM(attributes = c("ensembl_gene_id", 
#                               "chromosome_name",
#                               "start_position",
#                               "end_position",
#                               "external_gene_name", 
#                               "hgnc_symbol", 
#                               "description"), 
#                mart = ensembl, 
#                filters = "ensembl_gene_id", 
#                values = rownames(huh6.sce))
#write.csv(ens.bm, "data/ensembl_biomart.csv")
ens.bm <- read.csv("data/ensembl_biomart.csv", row.names = 1)

## Keep track of any duplicates
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)]
## Make a vector and give elements names at the same time
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])

## Two duplicates
## ENSG00000230417 = LINC00595
## ENSG00000276085 = CCL3L3

## Remove unnecessary annotation from the description
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)

## This filters out around 200 genes
nrow(ens.filt.bm) ## 36392
nrow(huh6.sce) ## 36613

## Keep only the simple chromosome information
table(ens.filt.bm$chromosome_name)
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]

## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]

## Now we've filtered out around 10,000 genes
nrow(ens.filt.bm) ## 25310
nrow(huh6.sce) ## 36613

## Keep genes that didn't get annotated so we know what gets dropped
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 11082

length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 10924 / 11082 are novel transcripts or novel transcripts antisense to a gene

grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest is a mixed bag: ~160 genes which are either novel, mitochondrially encoded, non-coding, or just empty description
## It seems unlikely to be particularly informative in our case

## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## ARMCX5-GPRASP2
## CCL3L1 - part of chemokine cluster on chr17
## DNAJC9-AS1 - heatshock protein homolog -AS1 is readthrough
## ELFN2 - postsynaptic cell adhesion molecule in group III mGluRs
## GOLGA8M - something to do with the Golgi membrane
## GPR84-AS1 - a g-protein coupled receptor complex member readthrough
## MATR3 - nuclear matrix protein          
## MKKS - centrosomal shuttling protein
## NPIPA9 - nuclear pore complex related
## PRICKLE2-AS1 - homolog of prickle readthrough
## RAET1E-AS1 - MHC class I related genes member readthrough
## SPATA13 - spermatogenesis associated
## TNFRSF10A-DT - TNF-receptor superfamily member; DT = divergent transcript (?)

## LINC00484
## LINC00595
## LINC01115
## LINC01238
## LINC01605
## LINC03021
## LINC03023
## LINC03025

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 25288

mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]

ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype=="lincRNA"]

linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Before we start chopping and changing, let's get some basic QC into the object
huh6.sce <- addPerCellQCMetrics(huh6.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))

huh6.sce <- huh6.sce[ens.filt.bm$ensembl_gene_id, ]
rownames(huh6.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

metadata.df <-  as.data.frame(colData(huh6.sce))

huh6.seurat <- CreateSeuratObject(counts = assay(huh6.sce, "counts"),
                                  assay = "RNA",
                                  meta.data = metadata.df)

saveRDS(huh6.sce, "data/huh6_sce_start.rds")
saveRDS(huh6.seurat, "data/huh6_seurat_start.rds")

## [ QC ] ----

huh6.sce <- readRDS("data/huh6_sce_start.rds")

dir.create("plots/qc/", recursive = TRUE)

table(huh6.sce$Condition)

huh6.sce$Condition <- factor(huh6.sce$Condition,
                             levels = c("POT", "Untreated", "Cisplatin_off", "Cisplatin_recovery1", "Cisplatin_recovery2", "Cisplatin_recovery3"))

det.sce.gg <- plotColData(huh6.sce, y = "detected", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

umi.sce.gg <- plotColData(huh6.sce, y = "total", x = "Condition", colour_by = "Condition") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")

mt.sce.gg <- plotColData(huh6.sce, y = "subsets_mt_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

rb.sce.gg <- plotColData(huh6.sce, y = "subsets_ribo_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.5)

ggsave("plots/qc/huh6_QC_plots.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/huh6_QC_plots.png", width = 8.7, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
## Use sparse matrix operations, as if your dataset is large, doing matrix operations the regular way will take a long time
counts.assay <- counts(huh6.sce)
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

ggsave("plots/qc/huh6_QC_topgenes.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/huh6_QC_topgenes.png", width = 8.7, height = 5.8)

## Remove the counts as it's not needed.
rm(counts.assay)
gc()

## MALAT1 is a long non-coding RNA and typically comes up enriched in skim-seq like 10x
## It is routinely removed as part of preprocessing
## However, since we're doing drug treatments, it's worth checking to see if it might indicate stress response
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7099402/

malat1.sce.gg <- plotExpression(huh6.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Condition", 
                                colour_by = "Condition") +
  labs(x = element_blank(), y = expression("Log"[2]*" Expression")) +
  scale_x_discrete(labels = conditions.list) +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

## During preprocessing we added a QC annotation for lincRNA

lincrna.sce.gg <- plotColData(huh6.sce, y = "subsets_linc_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.5, text = element_text(size = 15))

## Cisplatin is having some kind of increase in lincRNA, and this is increased in recovery 1
## However, decreases in recovery 2

ggsave("plots/qc/huh6_QC_lincrna_malat1.pdf", width = 8.7, height = 8.7)
ggsave("plots/qc/huh6_QC_lincrna_malat1.png", width = 8.7, height = 8.7)

lincrna.sce.gg <- lincrna.sce.gg +
  labs(title = "Proportion of lincRNA in total cell expression")

lincrna.sce.gg +
  theme(text = element_text(size = 14))

ggsave("plots/qc/huh6_lincrna.pdf", width = 5.8, height = 5.8)
ggsave("plots/qc/huh6_lincrna.png", width = 5.8, height = 5.8)

lincrna.sce.gg <- lincrna.sce.gg +
  labs(title = "Proportion of lincRNA\nin total cell expression") +
  theme(axis.text.x = element_text(size = 6))

mt.sce.gg + rb.sce.gg + lincrna.sce.gg + theme(aspect.ratio = 0.5)

ggsave("plots/qc/huh6_QC_mt_rb_lincrna.pdf", width = 11.6, height = 5.8)
ggsave("plots/qc/huh6_QC_mt_rb_lincrna.png", width = 11.6, height = 5.8)

gc()

## Cell cycle prediction
## Tricycle assigns cell cycle along a radial trajectory between 0 and 2pi
tricycle.sce <- project_cycle_space(huh6.sce, gname.type = "SYMBOL", species = "human")
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

huh6.seurat <- readRDS("data/huh6_seurat_start.rds")

huh6.cycle.seurat <- NormalizeData(huh6.seurat)
huh6.cycle.seurat <- CellCycleScoring(huh6.cycle.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

cell.cycle.df <- data.frame("Cell.ID" = colnames(huh6.sce),
                            "Seurat" = huh6.cycle.seurat$Phase,
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
huh6.sce$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
huh6.sce$Tricycle.Position <- tricycle.sce$tricyclePosition

huh6.sce$Seurat.Phase <- huh6.cycle.seurat$Phase
huh6.sce$Seurat.S <- huh6.cycle.seurat$S.Score
huh6.sce$Seurat.G2M <- huh6.cycle.seurat$G2M.Score

huh6.seurat$Seurat.Phase <- huh6.cycle.seurat$Phase
huh6.seurat$Seurat.S <- huh6.cycle.seurat$S.Score
huh6.seurat$Seurat.G2M <- huh6.cycle.seurat$G2M.Score

huh6.seurat$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
huh6.seurat$Tricycle.Position <- tricycle.sce$tricyclePosition

rm(tricycle.sce, huh6.cycle.seurat)
gc()

saveRDS(huh6.sce, "data/huh6_sce_working.rds")
saveRDS(huh6.seurat, "data/huh6_seurat_working.rds")

## [ QC filtering ] ----

huh6.sce <- readRDS("data/huh6_sce_working.rds")
huh6.seurat <- readRDS("data/huh6_seurat_working.rds")

colnames(colData(huh6.sce))
colnames(rowData(huh6.sce))

table(huh6.sce$Class)
## All singlet

table(huh6.sce$sample_id, huh6.sce$Condition)

sample.cols <- iwanthue(length(names(table(huh6.sce$sample_id))))
names(sample.cols) <- names(table(huh6.sce$sample_id))

det.sce.gg <- plotColData(huh6.sce, y = "detected", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dotted")

mt.sce.gg <- plotColData(huh6.sce, y = "subsets_mt_percent", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dotted")

rb.sce.gg <- plotColData(huh6.sce, y = "subsets_ribo_percent", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dotted")

det.sce.gg + mt.sce.gg + rb.sce.gg & theme(aspect.ratio = 0.5)

ggsave("plots/qc/huh6_QC_plots_sample.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/huh6_QC_plots_sample.png", width = 8.7, height = 5.8)

## Set the expression threshold over 1000
detected.filt <- colnames(huh6.sce)[huh6.sce$detected > 1000] ## Doesn't remove any cells
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(huh6.sce)[Matrix::rowSums(counts(huh6.sce)) > 5]
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(huh6.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

ut.sce <- huh6.sce[genes.filt, detected.filt]

dim(huh6.sce) - dim(ut.sce)
## We lose ~ 4000 genes

mito.filt <- colnames(ut.sce)[ut.sce$subsets_mt_percent < 10]
ribo.filt <- colnames(ut.sce)[ut.sce$subsets_ribo_percent > 5]
qc.filt <- intersect(mito.filt, ribo.filt)

ut.sce <- ut.sce[, qc.filt]

dim(huh6.sce) - dim(ut.sce)
## We lose ~ 3000 cells

## Percentage of the library that we filter out
round(100 - (100 * (table(ut.sce$sample_id) / table(huh6.sce$sample_id))), 1)
## HuH6_sample1 HuH6_sample10 HuH6_sample14 HuH6_sample16 HuH6_sample17 HuH6_sample18 HuH6_sample19  HuH6_sample2  
## 24.3          17.6           0.0          32.1          21.5          45.0          28.9          21.2     
## HuH6_sample3  HuH6_sample4  HuH6_sample5 HuH6_sample6  HuH6_sample7  HuH6_sample8  HuH6_sample9 
## 23.5          25.2           7.4          27.6          21.1          25.3          26.3 

table(huh6.filt.sce$sample_id)
## HuH6_sample1 HuH6_sample10 HuH6_sample14 HuH6_sample16 HuH6_sample17 HuH6_sample18 HuH6_sample19  HuH6_sample2
## 395          1836             5           146           278           94            734           368   
## HuH6_sample3  HuH6_sample4  HuH6_sample5 HuH6_sample6  HuH6_sample7  HuH6_sample8  HuH6_sample9 
## 364           243           831          1172          1728          1255          1342 

## The cisplatin recovery 1 and 2 samples got hit the hardest
## But relatively evenly hit looking at the percentages
round(100 - (100 * (table(huh6.filt.sce$Condition) / table(huh6.sce$Condition))), 1)
## POT           Untreated       Cisplatin_off Cisplatin_recovery1 Cisplatin_recovery2 Cisplatin_recovery3 
## 24.1                26.6                23.1                32.1                29.1                18.9 

table(huh6.filt.sce$Condition)
## POT           Untreated       Cisplatin_off Cisplatin_recovery1 Cisplatin_recovery2 Cisplatin_recovery3 
## 400                2638                3226                 146                 372                4009

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(huh6.seurat)[Matrix::rowSums(GetAssayData(huh6.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(huh6.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

## Seurat and Scater have their own methods for calling # of detected genes so these will disagree slightly
## To be consistent we'll filter on $detected rather than $nFeature_RNA
huh6.seurat$nFeature_RNA == huh6.seurat$detected

huh6.filt.seurat <- subset(huh6.seurat, 
                           cells = WhichCells(huh6.seurat, 
                                              expression = subsets_mt_percent < 10 &
                                                subsets_ribo_percent > 5 &
                                                detected > 1000),
                           features = genes.filt)

## The dimensions are equal
dim(huh6.filt.seurat) == dim(huh6.filt.sce)

saveRDS(huh6.filt.sce, "data/huh6_sce_filtered.rds")
saveRDS(huh6.filt.seurat, "datahuh6_seurat_filtered.rds")

## [ Filter HVGs ] ----

huh6.filt.seurat <- readRDS("data/huh6_seurat_filtered.rds")

## Filter HVGs
huh6.filt.seurat <- NormalizeData(huh6.filt.seurat)
huh6.filt.seurat <- FindVariableFeatures(huh6.filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(huh6.filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(huh6.filt.seurat), 20), repel = TRUE)
ggsave("plots/norm/filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/norm/filter_hvgs.png", width = 8.3, height = 5.8)

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
             "Normal" = huh6.filt.seurat$Seurat.Phase),
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
huh6.filt.seurat$Seurat.Cycle.Score <- huh6.filt.seurat$Seurat.S - huh6.filt.seurat$Seurat.G2M

saveRDS(huh6.filt.seurat, "data/huh6_seurat_hvgs.rds")

huh6.seurat <- readRDS("data/huh6_seurat_hvgs.rds")

huh6.seurat <- RunPCA(huh6.seurat, verbose = FALSE, npcs = 100, 
                      features = VariableFeatures(huh6.seurat))
ElbowPlot(object = huh6.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- RunHarmony(huh6.seurat, group.by.vars = "Sample",
                             theta = 1, lambda = 1, sigma = 0.07,
                             assay.use = "RNA", reduction = "pca",
                             dims.use = 1:30, reduction.save = "Harmony",
                             max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

harmony.seurat$Condition <- factor(harmony.seurat$Condition,
                                   levels = c("POT",
                                              "Untreated",
                                              "Cisplatin_off",
                                              "Cisplatin_recovery1",
                                              "Cisplatin_recovery2",
                                              "Cisplatin_recovery3"))

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = group.cols) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("plots/norm/umap_harmony_t1_l1_s007.png", width = 8.3, height = 5.8)
ggsave("plots/norm/umap_harmony_t1_l1_s007.pdf", width = 8.3, height = 5.8)

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
obj.list <- SplitObject(huh6.seurat, split.by = "Condition")

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
Idents(huh6.list$POT) <- huh6.list$POT$seurat_clusters.0.4
Idents(huh6.list$Untreated) <- huh6.list$Untreated$seurat_clusters.0.4
Idents(huh6.list$Cisplatin_off) <- huh6.list$Cisplatin_off$seurat_clusters.0.4
Idents(huh6.list$Cisplatin_recovery1) <- huh6.list$Cisplatin_recovery1$seurat_clusters.0.8
Idents(huh6.list$Cisplatin_recovery2) <- huh6.list$Cisplatin_recovery2$seurat_clusters.0.4
Idents(huh6.list$Cisplatin_recovery3) <- huh6.list$Cisplatin_recovery3$seurat_clusters.0.6
Idents(huh6.seurat) <- huh6.seurat$seurat_clusters.0.4

## Get cluster markers - MAST takes longer to run
pot.markers <- FindAllMarkers(huh6.list$POT, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
ut.markers <- FindAllMarkers(huh6.list$Untreated, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
cis.markers <- FindAllMarkers(huh6.list$Cisplatin_off, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
rec1.markers <- FindAllMarkers(huh6.list$Cisplatin_recovery1, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
rec2.markers <- FindAllMarkers(huh6.list$Cisplatin_recovery2, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
rec3.markers <- FindAllMarkers(huh6.list$Cisplatin_recovery3, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
huh6.all.markers <- FindAllMarkers(huh6.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

dir.create("data/jaccard/", recursive = TRUE)
saveRDS(pot.markers, "data/jaccard/pot_markers.rds")
saveRDS(ut.markers, "data/jaccard/ut_markers.rds")
saveRDS(cis.markers, "data/jaccard/cis_markers.rds")
saveRDS(rec1.markers, "data/jaccard/rec1_markers.rds")
saveRDS(rec2.markers, "data/jaccard/rec2_markers.rds")
saveRDS(rec3.markers, "data/jaccard/rec3_markers.rds")
saveRDS(huh6.all.markers, "data/jaccard/huh6_all_markers.rds")

## Filter significant and higher expressed markers
pot.markers <- pot.markers[pot.markers$p_val_adj < 0.05, ]
ut.markers <- ut.markers[ut.markers$p_val_adj < 0.05, ]
cis.markers <- cis.markers[cis.markers$p_val_adj < 0.05, ]
rec1.markers <- rec1.markers[rec1.markers$p_val_adj < 0.05, ]
rec2.markers <- rec2.markers[rec2.markers$p_val_adj < 0.05, ]
rec3.markers <- rec3.markers[rec3.markers$p_val_adj < 0.05, ]
huh6.all.markers <- huh6.all.markers[huh6.all.markers$p_val_adj < 0.05, ]

pot.markers <- pot.markers[pot.markers$avg_log2FC > 0.5, ]
ut.markers <- ut.markers[ut.markers$avg_log2FC > 0.5, ]
cis.markers <- cis.markers[cis.markers$avg_log2FC > 0.5, ]
rec1.markers <- rec1.markers[rec1.markers$avg_log2FC > 0.5, ]
rec2.markers <- rec2.markers[rec2.markers$avg_log2FC > 0.5, ]
rec3.markers <- rec3.markers[rec3.markers$avg_log2FC > 0.5, ]
huh6.all.markers <- huh6.all.markers[huh6.all.markers$avg_log2FC > 0.5, ]
## Split into lists of genes per cluster
pot.markers <- split(pot.markers, pot.markers$cluster)
ut.markers <- split(ut.markers, ut.markers$cluster)
cis.markers <- split(cis.markers, cis.markers$cluster)
rec1.markers <- split(rec1.markers, rec1.markers$cluster)
rec2.markers <- split(rec2.markers, rec2.markers$cluster)
rec3.markers <- split(rec3.markers, rec3.markers$cluster)
huh6.all.markers <- split(huh6.all.markers, huh6.all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
lapply(pot.markers, nrow) ## 0 = 565, 1 = 1591
lapply(ut.markers, nrow) ## 0 = 516, 1 = 636, 2 = 643, 3 = 1018, 4 = 844, 5 = 1210
lapply(cis.markers, nrow) ## 0 = 402, 1 = 239, 2 = 358, 3 = 945, 4 = 238, 5 = 660
lapply(rec1.markers, nrow) ## 0 = 34, 1 = 236, 2 = 32
lapply(rec2.markers, nrow) ## 0 = 362, 1 = 2377, 2 = 1174, 3 = 1467
lapply(rec3.markers, nrow) ## 0 = 228, 1 = 231, 2 = 1073, 3 = 357, 4 = 705, 5 = 435, 6 = 1780, 7 = 1002
lapply(huh6.all.markers, nrow) ## 0 = 548, 1 = 1176, 2 = 441, 3 = 869, 4 = 955, 5 = 1121, 6 = 199

pot.markers <- lapply(pot.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})
ut.markers <- lapply(ut.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})
cis.markers <- lapply(cis.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
rec1.markers <- lapply(rec1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 30)
})
rec2.markers <- lapply(rec2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
rec3.markers <- lapply(rec3.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
huh6.all.markers <- lapply(huh6.all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})

## Give your clusters unique names
names(pot.markers) <- paste0("POT_", names(pot.markers))
names(ut.markers) <- paste0("UT_", names(ut.markers))
names(cis.markers) <- paste0("Cis_", names(cis.markers))
names(rec1.markers) <- paste0("Rec1_", names(rec1.markers))
names(rec2.markers) <- paste0("Rec2_", names(rec2.markers))
names(rec3.markers) <- paste0("Rec3_", names(rec3.markers))
names(huh6.all.markers) <- paste0("HuH6_all_", names(huh6.all.markers))
## Make a big list of all the cluster markers
huh6.markers.list <- c(pot.markers, ut.markers, cis.markers,
                       rec1.markers, rec2.markers, rec3.markers,
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

rownames(huh6.jaccard.mat) <- gsub("HuH6_", "", rownames(huh6.jaccard.mat))
colnames(huh6.jaccard.mat) <- gsub("HuH6_", "", colnames(huh6.jaccard.mat))
rownames(huh6.jaccard.mat) <- gsub("all", "All", rownames(huh6.jaccard.mat))
colnames(huh6.jaccard.mat) <- gsub("all", "All", colnames(huh6.jaccard.mat))
anno.df <- data.frame("Timepoint" = gsub("_.+", "", colnames(huh6.jaccard.mat)),
                      row.names = colnames(huh6.jaccard.mat))
anno.col <- list("Timepoint" = c("All" = "black",
                                 "POT" = as.character(group.cols[1]),
                                 "UT" = as.character(group.cols[2]),
                                 "Cis" = as.character(group.cols[3]),
                                 "Rec1" = as.character(group.cols[4]),
                                 "Rec2" = as.character(group.cols[5]),
                                 "Rec3" = as.character(group.cols[6])))
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

## /////////////////////////////////////////////////////////////////////////////
## PAM clustering //////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

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
## For this data it was 11
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:16))
## Draw the geom_vline at 11 for this data
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
huh6.k11.pam <- pam(huh6.jaccard.dist, k = 11, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- huh6.k11.pam
## Assign the clusters back to your Seurat objects
huh6.clusters.df <- anno.df$PAM.Cluster
names(huh6.clusters.df) <- rownames(anno.df)
pot.set <- paste0("POT_", as.character(huh6.list$POT$seurat_clusters.0.4))
ut.set <- paste0("UT_", as.character(huh6.list$Untreated$seurat_clusters.0.4))
cis.set <- paste0("Cis_", as.character(huh6.list$Cisplatin_off$seurat_clusters.0.4))
rec1.set <- paste0("Rec1_", as.character(huh6.list$Cisplatin_recovery1$seurat_clusters.0.8))
rec2.set <- paste0("Rec2_", as.character(huh6.list$Cisplatin_recovery2$seurat_clusters.0.4))
rec3.set <- paste0("Rec3_", as.character(huh6.list$Cisplatin_recovery3$seurat_clusters.0.6))
all.set <- paste0("All_", as.character(huh6.seurat$seurat_clusters.0.4))

huh6.list$POT$PAM.Cluster <- factor(as.character(huh6.clusters.df[pot.set]))
huh6.list$Untreated$PAM.Cluster <- factor(as.character(huh6.clusters.df[ut.set]))
huh6.list$Cisplatin_off$PAM.Cluster <- factor(as.character(huh6.clusters.df[cis.set]))
huh6.list$Cisplatin_recovery1$PAM.Cluster <- factor(as.character(huh6.clusters.df[rec1.set]))
huh6.list$Cisplatin_recovery2$PAM.Cluster <- factor(as.character(huh6.clusters.df[rec2.set]))
huh6.list$Cisplatin_recovery3$PAM.Cluster <- factor(as.character(huh6.clusters.df[rec3.set]))
huh6.seurat$PAM.Cluster <- factor(as.character(huh6.clusters.df[all.set]))
saveRDS(huh6.list, "data/huh6_seurat_split_list_pam.rds")
saveRDS(huh6.seurat, "data/huh6_seurat_hvgs_pam.rds")

## [ Markers ] ----

## PAM cluster markers
Idents(huh6.seurat) <- huh6.seurat$PAM.Cluster
huh6.pam.markers <- FindAllMarkers(huh6.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
saveRDS(huh6.pam.markers, "data/jaccard/huh6_pam_markers.rds")

table(huh6.pam.markers$cluster)
## 1    2    3    4    5    7 
## 6588 2614 6007 7056 6807 6914

huh6.pam.list <- lapply(unique(huh6.pam.markers$cluster), function(x){
  tmp.df <- huh6.pam.markers[huh6.pam.markers$cluster==x & huh6.pam.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})
names(huh6.pam.list) <- unique(huh6.pam.markers$cluster)

saveRDS(huh6.pam.list, "data/jaccard/huh6_pam_markers_list.rds")

huh6.pam.GO.BP.list <- lapply(huh6.pam.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

huh6.pam.C2.list <- lapply(huh6.pam.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

huh6.pam.H.list <- lapply(huh6.pam.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(huh6.pam.GO.BP.list, "data/jaccard/huh6_pam_markers_BP_ontology.rds")
saveRDS(huh6.pam.C2.list, "data/jaccard/huh6_pam_markers_C2_geneset.rds")
saveRDS(huh6.pam.H.list, "data/jaccard/huh6_pam_markers_Hallmarks.rds")

huh6.pam.lists <- list("huh6.pam.GO.BP.list" = huh6.pam.GO.BP.list,
                       "huh6.pam.C2.list" = huh6.pam.C2.list,
                       "huh6.pam.H.list" = huh6.pam.H.list,
                       "huh6.pam.list" = huh6.pam.list)

huh6.pam.names <- sub("list", "sheets", names(huh6.pam.lists))

for(i in 1:length(huh6.pam.lists)){
  huh6.pam.list <- huh6.pam.lists[[i]]
  list.sheets <- lapply(huh6.pam.list, as.data.frame)
  assign(huh6.pam.names[[i]], list.sheets)
}

write_xlsx(huh6.pam.GO.BP.sheets, "data/jaccard/huh6_pam_markers_BP_ontology.xlsx")
write_xlsx(huh6.pam.C2.sheets, "data/jaccard/huh6_pam_markers_C2_geneset.xlsx")
write_xlsx(huh6.pam.H.sheets, "data/jaccard/huh6_pam_markers_Hallmarks.xlsx")
write_xlsx(huh6.pam.sheets, "data/jaccard/huh6_pam_markers.xlsx")

## [ Cluster annotation ] ----

huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam.rds")

umap.gg <- DimPlot(huh6.seurat, group.by = "Condition", order = TRUE) +
  umap.theme() + labs(title = "Condition") +
  scale_colour_manual(values = group.cols, labels = conditions.list)
umap.gg
ggsave("plots/jaccard/umap_harmony_conditions.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/umap_harmony_conditions.png", width = 8.3, height = 5.8)

umap.gg <- DimPlot(huh6.seurat, group.by = "PAM.Cluster", order = TRUE) +
  umap.theme() + labs(title = "PAM clusters")
umap.gg
ggsave("plots/jaccard/umap_harmony_pam_clusters.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/umap_harmony_pam_clusters.png", width = 8.3, height = 5.8)

umap.gg <- DimPlot(huh6.seurat, group.by = "seurat_clusters.0.4", order = TRUE) +
  umap.theme() + labs(title = "Seurat clusters")
umap.gg
ggsave("plots/jaccard/umap_harmony_seurat_clusters.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/umap_harmony_seurat_clusters.png", width = 8.3, height = 5.8)

## Use signature scores to annotate clusters
## HB signatures from literature
sig.df <- as.data.frame(read_xlsx(path = "data/hb_sigs.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")

sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})

names(sigs) <- unique(sig.df$Signature)
sigs

sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)

huh6.seurat <- AddModuleScore(huh6.seurat, features = sigs, assay = "RNA", seed = 12345, 
                              name = sig.names)
colnames(huh6.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(huh6.seurat@meta.data))

huh6.seurat[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = huh6.seurat, vars = sig.names)))

DoHeatmap(huh6.seurat, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "PAM.Cluster", label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/jaccard/hb_sigs_heatmap_pam.pdf", width = 5.8, height = 8.7)
ggsave("plots/jaccard/hb_sigs_heatmap_pam.png", width = 5.8, height = 8.7)

DoHeatmap(huh6.seurat, features = sig.names, assay = "HB_sigs", slot = "data",
          group.by = "seurat_clusters.0.4", label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Seurat"))
ggsave("plots/jaccard/hb_sigs_heatmap_clusters.pdf", width = 5.8, height = 8.7)
ggsave("plots/jaccard/hb_sigs_heatmap_clusters.png", width = 5.8, height = 8.7)

## HNF4A and LEF1 (Kluiver et al., 2023)
features = c("HNF4A", "LEF1")
FeaturePlot(huh6.seurat, features = features) &
  umap.theme()
ggsave("plots/jaccard/hnf4a_lef1_umap.pdf", width = 8.7, height = 5.8)
ggsave("plots/jaccard/hnf4a_lef1_umap.png", width = 8.7, height = 5.8)

Idents(huh6.seurat) <- huh6.seurat$PAM.Cluster
VlnPlot(huh6.seurat, features = features, pt.size = 0) &
  xlab("")
ggsave("plots/jaccard/hnf4a_lef1_violin_pam.pdf", width = 8.7, height = 5.8)
ggsave("plots/jaccard/hnf4a_lef1_violin_pam.png", width = 8.7, height = 5.8)

Idents(huh6.seurat) <- huh6.seurat$seurat_clusters.0.4
VlnPlot(huh6.seurat, features = features, pt.size = 0) &
  xlab("")
ggsave("plots/jaccard/hnf4a_lef1_violin_clusters.pdf", width = 8.7, height = 5.8)
ggsave("plots/jaccard/hnf4a_lef1_violin_clusters.png", width = 8.7, height = 5.8)

## GRN modules in Roegrig et al. (2024)
roehrig.sig <- read_xlsx("data/roehrig2024_sig.xlsx")
unique(roehrig.sig$`GRN module`)
states <- c("scH", "sc-epi", "scLP", "scM")

roehrig.filt <- roehrig.sig[(roehrig.sig$`GRN module` %in% states), ]
unique(roehrig.filt$`GRN module`)

grn.sigs <- lapply(unique(roehrig.filt$`GRN module`), function(x){
  roehrig.filt[roehrig.filt$`GRN module`==x, "Gene"]
})

names(grn.sigs) <- unique(roehrig.filt$`GRN module`)
grn.sigs

grn.names <- paste(names(grn.sigs), ".Sig", sep = "")

huh6.seurat <- AddModuleScore(huh6.seurat, features = grn.sigs, assay = "RNA", seed = 12345,
                              name = grn.names)
colnames(huh6.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(huh6.seurat@meta.data))

huh6.seurat[["Roehrig_sigs"]] <- CreateAssayObject(data = t(FetchData(object = huh6.seurat, vars = grn.names)))

DoHeatmap(huh6.seurat, features = grn.names, assay = "Roehrig_sigs", slot = "data",
          group.by = "PAM.Cluster", label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/jaccard/roehrig_sigs_heatmap_pam.pdf", width = 5.8, height = 8.7)
ggsave("plots/jaccard/roehrig_sigs_heatmap_pam.png", width = 5.8, height = 8.7)

DoHeatmap(huh6.seurat, features = grn.names, assay = "Roehrig_sigs", slot = "data",
          group.by = "seurat_clusters.0.4", label = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Seurat"))
ggsave("plots/jaccard/roehrig_sigs_heatmap_clusters.pdf", width = 5.8, height = 8.7)
ggsave("plots/jaccard/roehrig_sigs_heatmap_clusters.png", width = 5.8, height = 8.7)

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

huh6.seurat <- AddModuleScore(huh6.seurat, features = H.list, assay = "RNA", seed = 12345, 
                              name = H.names)
colnames(huh6.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(huh6.seurat@meta.data))

huh6.seurat[["MSigDB_H"]] <- CreateAssayObject(data = t(FetchData(object = huh6.seurat, vars = H.names)))

DoHeatmap(huh6.seurat, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/sig/msigdb_H_heatmap_pam.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/msigdb_H_heatmap_pam.png", width = 16.5, height = 8.7)

DoHeatmap(huh6.seurat, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "seurat_clusters.0.4", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Seurat"))
ggsave("plots/sig/msigdb_H_heatmap_clusters.pdf", width = 16.5, height = 8.7)
ggsave("plots/sig/msigdb_H_heatmap_clusters.png", width = 16.5, height = 8.7)

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

huh6.seurat <- AddModuleScore(huh6.seurat, features = liver.list, assay = "RNA", seed = 12345, 
                              name = liver.names)
colnames(huh6.seurat@meta.data) <- gsub("\\.Sig[0-9]+$", "\\.Sig", colnames(huh6.seurat@meta.data))

huh6.seurat[["MSigDB_C2"]] <- CreateAssayObject(data = t(FetchData(object = huh6.seurat, vars = liver.names)))

DoHeatmap(huh6.seurat, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "PAM.Cluster", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "PAM"))
ggsave("plots/jaccard/msigdb_C2_liver_heatmap_pam.pdf", width = 16.5, height = 8.7)
ggsave("plots/jaccard/msigdb_C2_liver_heatmap_pam.png", width = 16.5, height = 8.7)

DoHeatmap(huh6.seurat, features = liver.names, assay = "MSigDB_C2", slot = "data",
          group.by = "seurat_clusters.0.4", size = 4) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = "Seurat"))
ggsave("plots/jaccard/msigdb_C2_liver_heatmap_clusters.pdf", width = 16.5, height = 8.7)
ggsave("plots/jaccard/msigdb_C2_liver_heatmap_clusters.png", width = 16.5, height = 8.7)

saveRDS(huh6.seurat, "data/huh6_seurat_hvgs_pam_sigs.rds")

## /////////////////////////////////////////////////////////////////////////////
## Pseudobulk //////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Pseudobulk across PAM clusters and samples and look at expression of same pathways/markers as above
pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Condition", "PAM.Cluster"))
Idents(pseudo.huh6) <- "PAM.Cluster"

bulk.markers <- FindAllMarkers(pseudo.huh6, test.use = "DESeq2")
saveRDS(bulk.markers, "data/pseudobulk_pam_markers.rds")

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
ggsave("plots/jaccard/pseudo_msigdb_H_heatmap_pam.pdf", width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_msigdb_H_heatmap_pam.png", width = 11.7, height = 8.3)

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
ggsave("plots/jaccard/pseudo_msigdb_H_pheatmap_pam.pdf", plot = gt)
ggsave("plots/jaccard/pseudo_msigdb_H_pheatmap_pam.png", plot = gt)

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
ggsave("plots/jaccard/pseudo_msigdb_C2_liver_heatmap_pam.pdf", width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_msigdb_C2_liver_heatmap_pam.png", width = 11.7, height = 8.3)

huh6.c2.mat <- as.matrix(pseudo.huh6@assays$MSigDB_C2@data)
huh6.reactome.mat <- huh6.c2.mat[which(grepl("REACTOME", rownames(huh6.c2.mat))), ]

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
ggsave("plots/jaccard/pseudo_msigdb_C2_reactome_pheatmap_pam.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_msigdb_C2_reactome_pheatmap_pam.png", plot = gt, width = 11.7, height = 8.3)

huh6.cairo.mat <- huh6.c2.mat[which(grepl("CAIRO", rownames(huh6.c2.mat))), ]

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
ggsave("plots/jaccard/pseudo_msigdb_C2_cairo_pheatmap_pam.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_msigdb_C2_cairo_pheatmap_pam.png", plot = gt, width = 11.7, height = 8.3)

saveRDS(pseudo.huh6, "data/huh6_seurat_pseudobulk_sigs.rds")

## Liver differentiation markers
## (https://pmc.ncbi.nlm.nih.gov/articles/PMC4999623/; https://pmc.ncbi.nlm.nih.gov/articles/PMC11114060/)
custom <- c("AFP", "ALB", "CEBPA", "EPCAM", "FOXA1", "FOXA2", "GATA4",
            "HNF1A", "HNF1B", "HNF4A", "ICAM1", "KRT19", "MAT1A", "MAT2A",
            "NCAM1", "NR1I2", "NR5A2", "ONECUT2", "PROM1", "PROX1",
            "SOX9", "TAF10", "TAF4", "TBP", "TBX3", "TGFB1", "TGFB2")
custom %in% rownames(pseudo.huh6)

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/jaccard/pseudo_custom_pheatmap.pdf", plot = gt)
ggsave("plots/jaccard/pseudo_custom_pheatmap.png", plot = gt)

## Pseudobulk across clusters only and look at expression of same pathways/markers as above
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_sigs.rds")

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("PAM.Cluster"))
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

DoHeatmap(pseudo.huh6, features = H.names, assay = "MSigDB_H", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))
ggsave("plots/sig/pseudo_pam_msigdb_H_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_msigdb_H_heatmap.png", width = 11.7, height = 8.3)

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
ggsave("plots/jaccard/pseudo_pam_msigdb_H_pheatmap.pdf", plot = gt)
ggsave("plots/jaccard/pseudo_pam_msigdb_H_pheatmap.png", plot = gt)

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
ggsave("plots/sig/pseudo_pam_msigdb_C2_liver_heatmap.pdf", width = 11.7, height = 8.3)
ggsave("plots/sig/pseudo_pam_msigdb_C2_liver_heatmap.png", width = 11.7, height = 8.3)

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
ggsave("plots/jaccard/pseudo_pam_msigdb_C2_reactome_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_pam_msigdb_C2_reactome_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

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
ggsave("plots/jaccard/pseudo_pam_msigdb_C2_cairo_pheatmap.pdf", plot = gt, width = 11.7, height = 8.3)
ggsave("plots/jaccard/pseudo_pam_msigdb_C2_cairo_pheatmap.png", plot = gt, width = 11.7, height = 8.3)

saveRDS(pseudo.huh6, "data/huh6_seurat_pseudobulk_pam_sigs.rds")

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
custom %in% rownames(pseudo.huh6)

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])

gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/jaccard/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggsave("plots/jaccard/pseudo_pam_custom_pheatmap.png", plot = gt)

## Liver differentiation signature (Wesley et al., 2022)
wesley.df <- as.data.frame(read_xlsx(path = "data/welsey2022_sig.xlsx", col_names = FALSE))
colnames(wesley.df) <- c("Gene", "Signature")

wesley.sig <- lapply(unique(wesley.df$Signature), function(x){
  wesley.df[wesley.df$Signature==x, "Gene"]
})
names(wesley.sig) <- unique(wesley.df$Signature)
wesley.sig

wesley.names <- paste(names(wesley.sig), ".Sig", sep = "")

pseudo.huh6 <- AddModuleScore(pseudo.huh6, features = wesley.sig, assay = "RNA", seed = 12345, 
                              name = wesley.names)
colnames(pseudo.huh6@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.huh6@meta.data))

pseudo.huh6[["Wesley_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.huh6, vars = wesley.names)))

DoHeatmap(pseudo.huh6, features = wesley.names, assay = "Wesley_sigs", slot = "data",
          group.by = "PAM.Cluster", size = 4, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 12), legend.key.size = unit(1, "cm"), axis.text.y = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = ""))

huh6.wes.mat <- as.matrix(pseudo.huh6@assays$Wesley_sigs@data)

anno.df <- pseudo.huh6@meta.data["PAM.Cluster"]
gt <- pheatmap(huh6.wes.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df)$gtable
ggsave("plots/jaccard/pseudo_pam_wesley_pheatmap.pdf", plot = gt)
ggsave("plots/jaccard/pseudo_pam_wesley_pheatmap.png", plot = gt)

## Annotate PAM clusters
huh6.seurat$PAM_Name <- recode(as.character(huh6.seurat$PAM.Cluster),
                               "1" = "Hepatocytic 2",
                               "2" = "Stem-like",
                               "3" = "Hepatocytic 3",
                               "4" = "Hepatocytic 1",
                               "5" = "Late progenitor",
                               "7" = "Early progenitor")
saveRDS(huh6.seurat, "data/huh6_seurat_hvgs_pam_anno.rds")

umap.gg <- DimPlot(huh6.seurat, group.by = "PAM_Name", order = TRUE) +
  umap.theme() +
  labs(title = "PAM clusters") +
  scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15))
umap.gg
ggsave("plots/jaccard/umap_harmony_pam_clusters_anno.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/umap_harmony_pam_clusters_anno.png", width = 8.3, height = 5.8)

anno.data <- huh6.seurat@meta.data[c("description", "Condition", "PAM_Name")]
anno.table <- as.data.frame(table(anno.data$Condition, anno.data$PAM_Name))

ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "PAM cluster name")) +
  scale_x_discrete(labels = conditions.list) +
  scale_fill_manual(values = pam.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/jaccard/pam_cluster_proportion_condition.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_cluster_proportion_condition.png", width = 8.3, height = 5.8)

anno.table <- as.data.frame(table(anno.data$description, anno.data$PAM_Name))
anno.table$Var1 <- factor(anno.table$Var1,
                          levels = c("POT", "Untreated_A", "Untreated_B", "Untreated_C",
                                     "Cisplatin_off_A", "Cisplatin_off_B", "Cisplatin_off_C",
                                     "Cisplatin_recovery1_A", "Cisplatin_recovery1_B", "Cisplatin_recovery1_C",
                                     "Cisplatin_recovery2_A", "Cisplatin_recovery2_B", "Cisplatin_recovery2_C",
                                     "Cisplatin_recovery3_A", "Cisplatin_recovery3_B", "Cisplatin_recovery3_C"))

ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "PAM cluster name")) +
  scale_fill_manual(values = pam.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/jaccard/pam_cluster_proportion_replicate.pdf", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_cluster_proportion_replicate.png", width = 8.3, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Immunofluorescence targets and markers //////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("Condition"))
Idents(pseudo.huh6) <- "Condition"

custom <- c("PROX1", ## HB proliferation/migration
            "NCAM1", "CLDN3", ## HSC marker
            "EPCAM", ## HSC and hepatoblast marker
            "ICAM1", "KRT7", ## Hepatoblast
            "AFP", ## Hepatoblast and hepatocyte
            "HNF4A", ## Hepatocyte marker
            "SOX9", ## Hepatoblast and cholangiocyte marker
            "HNF1B", ## Cholangiocyte marker
            "AURKA", "IGF2", "CORIN", ## Hepatocytic cluster 
            "SHH", "KDM5B", ## Stem-like cluster 
            "PMEPA1", "NNMT", ## Early progenitor cluster
            "TGFB2") ## Late progenitor cluster
custom %in% rownames(pseudo.huh6)

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])
anno.df <- pseudo.huh6@meta.data["Condition"]
anno.cols <- group.cols
names(anno.cols) <- gsub("_", "-", names(anno.cols))
gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(Condition = anno.cols))$gtable
ggsave("plots/targets/condition_pheatmap.pdf", plot = gt)
ggsave("plots/targets/condition_pheatmap.png", plot = gt)

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = c("PAM_Name"))
Idents(pseudo.huh6) <- "PAM_Name"

huh6.custom.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data[custom, ])
anno.df <- pseudo.huh6@meta.data["PAM_Name"]
anno.cols <- pam.cols
gt <- pheatmap(huh6.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = list(PAM_Name = anno.cols))$gtable
ggsave("plots/targets/cluster_pheatmap.pdf", plot = gt)
ggsave("plots/targets/cluster_pheatmap.png", plot = gt)

## [ Entropy ] ----

huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")
huh6.sce <- as.SingleCellExperiment(huh6.seurat)
saveRDS(huh6.sce, "data/huh6_sce_hvgs_pam_anno.rds")
huh6.sce <- readRDS("data/huh6_sce_hvgs_pam_anno.rds")

dir.create("plots/entropy/", recursive = TRUE)

## Entropy is a measure of gene expression variability
## Higher entropy = greater diversity
## It can be used to root trajectories often to the most undifferentiated state
## but in this case we want to root it to the most differentiated state

## Computing the entropy of each cell's expression profile
entropy <- perCellEntropy(huh6.sce)
entropy.df <- data.frame(cluster = huh6.sce$PAM_Name, entropy = entropy)

ggplot(entropy.df, aes(x = cluster, y = entropy, fill = cluster)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size = 12)) +
  scale_fill_manual(values = pam.cols) +
  labs(x = "Cluster", y = "Entropy")
ggsave("plots/entropy/entropy_cluster_violin.pdf", width = 8.3, height = 5.8)
ggsave("plots/entropy/entropy_cluster_violin.png", width = 8.3, height = 5.8)

## [ Trajectory inference ] ----

## Trajectory inference was performed using Slingshot
huh6.slingshot.sce <- slingshot(huh6.sce, reducedDim = "PCA",
                                clusterLabels = "PAM_Name", start.clus = "HB-like",
                                stretch = 0)
saveRDS(huh6.slingshot.sce, "data/huh6_sce_slingshot.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(huh6.slingshot.sce)
head(slingshot.lineages) ## 3 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

umap.gg <- plotUMAP(huh6.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(huh6.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = umap_1, y=umap_2), linewidth = 1.2)
}
umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))
ggsave("plots/slingshot/umap_slingshot.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/umap_slingshot.png", width = 8.3, height = 5.8)

saveRDS(slingshot.embedded.all, "data/slingshot_embedded_all.rds")

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(huh6.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/slingshot_sds_embedded_all.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(huh6.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot/slingshot_pseudotime_lineages.pdf", width = 8.3, height = 8.3)
png("plots/slingshot/slingshot_pseudotime_lineages.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(huh6.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(huh6.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(huh6.slingshot.sce, colour_by = "PAM_Name") +
  umap.theme() + scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

dir.create("plots/slingshot/by_lineage/", recursive = TRUE)

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/slingshot_embedded_", lineage.names[[i]], ".rds"))
  
  plotUMAP(huh6.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 5.8)
}

## [ Gene dynamics ] ----

slingshot.sce <- readRDS("data/huh6_sce_slingshot.rds")
slingshot.cols <- grep("^slingPseudotime_", colnames(colData(slingshot.sce)), value = TRUE)
cluster.cols <- pam.cols
names(cluster.cols) <- names(pam.cols)

for (pseudotime in slingshot.cols) {
  lineage <- TSCAN::testPseudotime(slingshot.sce, pseudotime = slingshot.sce[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/huh6_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/huh6_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(slingshot.sce, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "PAM_Name") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(slingshot.sce[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(slingshot.sce[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "PAM_Name",
                                  column_annotation_colours = list(PAM_Name = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/huh6_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/huh6_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(slingshot.sce, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "PAM_Name") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(slingshot.sce[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "PAM_Name",
                                column_annotation_colours = list(PAM_Name = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/slingshot/expression/huh6_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## [ Gene annotation ] ----

file.paths <- list.files(path = here::here("data/trajectory_genes/"), pattern = "\\.rds", full.names = TRUE)
file.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(file.paths)))
dynamic.genes <- lapply(file.paths, readRDS)
names(dynamic.genes) <- file.names

## GSEA of dynamic genes for each lineage
dynamic.genes.GO.BP.list <- lapply(dynamic.genes, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

dynamic.genes.C2.list <- lapply(dynamic.genes, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

dynamic.genes.H.list <- lapply(dynamic.genes, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})
saveRDS(dynamic.genes.GO.BP.list, "data/trajectory_genes/huh6_dynamic_genes_BP_ontology.rds")
saveRDS(dynamic.genes.C2.list, "data/trajectory_genes/huh6_dynamic_genes_C2_geneset.rds")
saveRDS(dynamic.genes.H.list, "data/trajectory_genes/huh6_dynamic_genes_Hallmarks.rds")

dynamic.genes.lists <- list("dynamic.genes.GO.BP.list" = dynamic.genes.GO.BP.list,
                            "dynamic.genes.C2.list" = dynamic.genes.C2.list,
                            "dynamic.genes.H.list" = dynamic.genes.H.list)
dynamic.genes.names <- sub("list", "sheets", names(dynamic.genes.lists))

for(i in 1:length(dynamic.genes.lists)){
  gs.list <- dynamic.genes.lists[[i]]
  list.sheets <- lapply(gs.list, as.data.frame)
  assign(dynamic.genes.names[[i]], list.sheets)
}
write_xlsx(dynamic.genes.GO.BP.sheets, "data/trajectory_genes/huh6_dynamic_genes_BP_ontology.xlsx")
write_xlsx(dynamic.genes.C2.sheets, "data/trajectory_genes/huh6_dynamic_genes_C2_geneset.xlsx")
write_xlsx(dynamic.genes.H.sheets, "data/trajectory_genes/huh6_dynamic_genes_Hallmarks.xlsx")

## Plot GSEA results for hallmark and gene ontology biological processes gene sets
for (i in 1:length(dynamic.genes.GO.BP.list)) {
  dotplot(dynamic.genes.GO.BP.list[[i]], showCategory = 20) +
    ggtitle(paste0("GSEA ", gsub("_", " ", names(dynamic.genes.GO.BP.list)[i])))
  ggsave(paste0("plots/slingshot/expression/huh6_gobp_", names(dynamic.genes.GO.BP.list)[i], ".pdf"), width = 8.3, height = 8.3)
  ggsave(paste0("plots/slingshot/expression/huh6_gobp_", names(dynamic.genes.GO.BP.list)[i], ".png"), width = 8.3, height = 8.3)
}

for (i in 1:length(dynamic.genes.H.list)) {
  dotplot(dynamic.genes.H.list[[i]], showCategory = 20) +
    ggtitle(paste0("GSEA ", gsub("_", " ", names(dynamic.genes.H.list)[i]))) +
    theme(text = element_text(size = 20),
          axis.text.y = element_text(size = 15))
  ggsave(paste0("plots/slingshot/expression/huh6_hallmarks_", names(dynamic.genes.H.list)[i], ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/expression/huh6_hallmarks_", names(dynamic.genes.H.list)[i], ".png"), width = 8.3, height = 5.8)
}

## /////////////////////////////////////////////////////////////////////////////
## TF enrichment ///////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

file.paths <- list.files(path = "data/trajectory_genes/", pattern = "\\.rds", full.names = TRUE)
file.paths <- file.paths[grep("slingpseudotime", file.paths, ignore.case = TRUE)]
file.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(file.paths)))
dynamic.genes <- lapply(file.paths, readRDS)
names(dynamic.genes) <- file.names
gene.lists <- lapply(dynamic.genes, `[[`, "gene")

data(motifAnnotations_hgnc)
# motifRankings <- importRankings("~/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motifRankings <- importRankings("~/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

motif.enrichment <- cisTarget(gene.lists, motifRankings,
                              motifAnnot = motifAnnotations)
saveRDS(motif.enrichment, "data/trajectory_genes/huh6_motif_enrichment_500bp.rds")
split.list <- split(motif.enrichment, motif.enrichment$geneSet)
write_xlsx(split.list, "data/trajectory_genes/huh6_motif_enrichment_500bp.xlsx")

# results.subset <- motif.enrichment[motif.enrichment$geneSet == "pseudotime2_down", ]
# showLogo(results.subset)

down.list <- split.list[grep("_down", split.list, ignore.case = TRUE)]
up.list <- split.list[grep("_up", split.list, ignore.case = TRUE)]

clean_and_split <- function(df) {
  cleaned <- gsub("\\(.*?\\)", "", df$TF_highConf)  
  split_genes <- unlist(strsplit(cleaned, "\\.|;")) 
  split_genes <- unique(trimws(split_genes[split_genes != ""]))
  split_genes[nchar(split_genes) > 0]
}

down.tfs <- lapply(down.list, clean_and_split)
down.common <- Reduce(intersect, down.tfs)
up.tfs <- lapply(up.list, clean_and_split)
up.common <- Reduce(intersect, up.tfs)

saveRDS(down.common, "data/trajectory_genes/tfs_common_down.rds")
saveRDS(up.common, "data/trajectory_genes/tfs_common_up.rds")

## Expression of TFs of interest
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")
up.common <- readRDS("data/trajectory_genes/tfs_common_up.rds")

filt.tfs <- up.common[up.common %in% rownames(huh6.seurat)]

pdf("plots/slingshot/expression/tfs_expression.pdf", width = 8.3, height = 5.8)
for (tf in filt.tfs) {
  print(VlnPlot(huh6.seurat, features = tf, group.by = "Condition", pt.size = 0) +
    scale_x_discrete(labels = conditions.list) +
    scale_fill_manual(values = group.cols, labels = conditions.list) +
    xlab(""))
}
dev.off()

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = "Condition")
Idents(pseudo.huh6) <- "Condition"
huh6.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data)
plot.mat <- huh6.mat[rownames(huh6.mat) %in% up.common, ]
anno.df <- data.frame(Condition = c("POT", "Untreated", "Cisplatin", "Recovery t1", "Recovery t2", "Recovery t3"),
                      row.names = colnames(pseudo.huh6))
anno.cols <- list(Condition = group.cols)
names(anno.cols$Condition) <- c("POT", "Untreated", "Cisplatin", "Recovery t1", "Recovery t2", "Recovery t3")

gt <- pheatmap(plot.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               cluster_cols = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = anno.cols)$gtable
ggsave("plots/slingshot/expression/tfs_expression_heatmap.pdf", plot = gt)
ggsave("plots/slingshot/expression/tfs_expression_heatmap.png", plot = gt)

pdf("plots/slingshot/expression/tfs_expression_cluster.pdf", width = 8.3, height = 5.8)
for (tf in filt.tfs) {
  print(VlnPlot(huh6.seurat, features = tf, group.by = "PAM_Name", pt.size = 0) +
          scale_fill_manual(values = pam.cols) +
          xlab(""))
}
dev.off()

pseudo.huh6 <- AggregateExpression(huh6.seurat, assays = "RNA", return.seurat = T,
                                   group.by = "PAM_Name")
Idents(pseudo.huh6) <- "PAM_Name"
huh6.mat <- as.matrix(pseudo.huh6@assays$RNA$scale.data)
plot.mat <- huh6.mat[rownames(huh6.mat) %in% up.common, ]
anno.df <- data.frame(Cluster = colnames(pseudo.huh6),
                      row.names = colnames(pseudo.huh6))
anno.cols <- list(Cluster = pam.cols)

gt <- pheatmap(plot.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               cluster_cols = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_col = anno.df,
               annotation_colors = anno.cols)$gtable
ggsave("plots/slingshot/expression/tfs_expression_heatmap_cluster.pdf", plot = gt)
ggsave("plots/slingshot/expression/tfs_expression_heatmap_cluster.png", plot = gt)

## Filter pseudotime genes for targets of enriched TFs
## Outputs from: https://www.grnpedia.org/trrust/
dir <- "data/trajectory_genes/trrust/"
tsv.files <- list.files(dir, pattern = "\\.tsv$", full.names = TRUE)

trrust.df <- tsv.files %>%
  lapply(function(file) read_tsv(file, col_types = cols(.default = "c"))) %>%
  bind_rows()
saveRDS(trrust.df, "data/trajectory_genes/trrust/trrust_ref.rds")

huh6.paths <- list.files(path = "data/trajectory_genes/", pattern = "\\.rds", full.names = TRUE)
huh6.paths <- huh6.paths[grep("slingpseudotime[0-9]_up", huh6.paths, ignore.case = TRUE)]
huh6.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(huh6.paths)))
huh6.list <- lapply(huh6.paths, readRDS)
names(huh6.list) <- huh6.names

huh6.list <- lapply(huh6.list, function(df) {
  df$abs_logFC <- abs(df$logFC)
  return(df)
})
huh6.markers <- lapply(huh6.list, \(x) {
  x <- x[order(x$abs_logFC, decreasing = TRUE), ]
  head(x, 5000)
})
huh6.markers <- lapply(huh6.markers, `[[`, "gene")
common.markers <- Reduce(intersect, huh6.markers)
test.genes <- common.markers[common.markers %in% trrust.df$Target]
saveRDS(test.genes, "data/trajectory_genes/test_genes.rds")

## /////////////////////////////////////////////////////////////////////////////
## Trajectory gene plots ///////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Expression plots
huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")

huh6.paths <- list.files(path = "data/trajectory_genes/", pattern = "\\.rds", full.names = TRUE)
huh6.paths <- huh6.paths[grep("slingpseudotime", huh6.paths, ignore.case = TRUE)]
huh6.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(huh6.paths)))
huh6.list <- lapply(huh6.paths, readRDS)
names(huh6.list) <- huh6.names

huh6.list <- lapply(huh6.list, function(df) {
  df$abs_logFC <- abs(df$logFC)
  return(df)
})
huh6.markers <- lapply(huh6.list, \(x) {
  x <- x[order(x$abs_logFC, decreasing = TRUE), ]
  head(x, 5000)
})
huh6.markers <- lapply(huh6.markers, `[[`, "gene")

marker.names <- paste(names(huh6.markers), ".sig", sep = "")
marker.names <- gsub("_", "-", marker.names)

huh6.seurat <- AddModuleScore(huh6.seurat, features = huh6.markers, assay = "RNA", seed = 12345, 
                              name = marker.names)
colnames(huh6.seurat@meta.data) <- gsub("\\.sig[0-9]+$", "\\.sig", colnames(huh6.seurat@meta.data))

Idents(huh6.seurat) <- huh6.seurat$Condition
VlnPlot(huh6.seurat, features = marker.names,
        pt.size = 0, cols = group.cols, ncol = 2) &
  xlab("")
ggsave("plots/slingshot/expression/violin_condition.pdf", width = 8.7, height = 11.3)
ggsave("plots/slingshot/expression/violin_condition.png", width = 8.7, height = 11.3)

Idents(huh6.seurat) <- huh6.seurat$PAM_Name
VlnPlot(huh6.seurat, features = marker.names,
        pt.size = 0, cols = pam.cols, ncol = 2) &
  xlab("")
ggsave("plots/slingshot/expression/violin_cluster.pdf", width = 8.7, height = 11.3)
ggsave("plots/slingshot/expression/violin_cluster.png", width = 8.7, height = 11.3)

## Similarity between HuH6 and HepG2 trajectory gene lists
huh6.paths <- list.files(path = "data/trajectory_genes/", pattern = "\\.rds", full.names = TRUE)
huh6.paths <- huh6.paths[grep("slingpseudotime", huh6.paths, ignore.case = TRUE)]
huh6.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(huh6.paths)))
huh6.list <- lapply(huh6.paths, readRDS)
names(huh6.list) <- huh6.names

hepg2.paths <- list.files(path = "/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/hepg2/data/slingshot_late_progenitor/trajectory_genes", pattern = "\\.rds", full.names = TRUE)
hepg2.paths <- hepg2.paths[grep("slingpseudotime1", hepg2.paths, ignore.case = TRUE)]
hepg2.names <-  gsub("_genes.rds$", "",
                     gsub("hepg2_sling", "", x = basename(hepg2.paths)))
hepg2.list <- lapply(hepg2.paths, readRDS)
names(hepg2.list) <- hepg2.names

lapply(huh6.list, nrow)
## pseudotime1_down = 7818, pseudotime1_up = 8221, 
## pseudotime2_down = 7859, pseudotime2_up = 7749,
## pseudotime3_down = 8196, pseudotime3_up = 7016
lapply(hepg2.list, nrow)
## pseudotime1_down = 9858, pseudotime1_up = 6293

## Convert Log2FC values to absolute values
huh6.list <- lapply(huh6.list, function(df) {
  df$abs_logFC <- abs(df$logFC)
  return(df)
})
hepg2.list <- lapply(hepg2.list, function(df) {
  df$abs_logFC <- abs(df$logFC)
  return(df)
})

## Take top 5000 per pseudotime lineage
huh6.markers <- lapply(huh6.list, \(x) {
  x <- x[order(x$abs_logFC, decreasing = TRUE), ]
  head(x, 5000)
})
hepg2.markers <- lapply(hepg2.list, \(x) {
  x <- x[order(x$abs_logFC, decreasing = TRUE), ]
  head(x, 5000)
})

names(huh6.markers) <- paste0("HuH6_", names(huh6.markers))
names(hepg2.markers) <- paste0("HepG2_", names(hepg2.markers))
hb.markers.list <- c(huh6.markers, hepg2.markers)

## Jaccard similarity
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

rownames(hb.jaccard.mat) <- gsub("_", " ", rownames(hb.jaccard.mat))
colnames(hb.jaccard.mat) <- gsub("_", " ", colnames(hb.jaccard.mat))

anno.df <- data.frame("Type" = gsub("HuH6 pseudotime(\\d+) ", "",
                                    gsub("HepG2 pseudotime(\\d+) ", "",
                                         gsub(" \\d+$", "", colnames(hb.jaccard.mat)))),
                      row.names = colnames(hb.jaccard.mat))
pal <- paletteer_d("PNWColors::Shuksan2")
anno.col <- list("Type" = c("up" = pal[5],
                            "down" = pal[1]))

dev.off()
gt <- pheatmap(hb.jaccard.mat,
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
ggsave("plots/jaccard/hb_jaccard_heatmap_pseudotime.png", plot = gt)
ggsave("plots/jaccard/hb_jaccard_heatmap_pseudotime.pdf", plot = gt)

## Overlap with HepG2 trajectory genes
all.tfs <- readRDS("data/trajectory_genes/huh6_hepg2_up_tfs.rds")
all.genes <- unique(unlist(all.tfs))
gene.matrix <- data.frame(Gene = all.genes,
                          A = as.integer(all.genes %in% all.tfs[[1]]),
                          B = as.integer(all.genes %in% all.tfs[[2]]),
                          C = as.integer(all.genes %in% all.tfs[[3]]),
                          D = as.integer(all.genes %in% all.tfs[[4]]))
colnames(gene.matrix)[-1] <- names(all.tfs)

png("plots/figures/shared_tfs_upset.png", units = "in", width = 8.7, height = 5.8, res = 200)
pdf("plots/figures/shared_tfs_upset.pdf", width = 8.7, height = 5.8)
upset(gene.matrix,
      sets = names(all.tfs),
      order.by = "freq",   
      mainbar.y.label = "Gene Set Intersections",    
      sets.x.label = "Number of Genes",   
      text.scale = 1.5,   
      keep.order = TRUE)
dev.off()

## Roehrig data
roehrig.seurat <- readRDS("data/seurat_object.RDS")

DimPlot(roehrig.seurat, group.by = "sample_identity", order = TRUE) +
  umap.theme() + labs(title = "Sample")

DimPlot(roehrig.seurat, group.by = "cell_identity", order = TRUE) +
  umap.theme() + labs(title = "Cell type")

DimPlot(roehrig.seurat, group.by = "H_LP_M", order = TRUE) +
  umap.theme() + labs(title = "Tumour cell state")

huh6.paths <- list.files(path = "data/trajectory_genes/", pattern = "\\.rds", full.names = TRUE)
huh6.paths <- huh6.paths[grep("slingpseudotime", huh6.paths, ignore.case = TRUE)]
huh6.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(huh6.paths)))
huh6.list <- lapply(huh6.paths, readRDS)
names(huh6.list) <- huh6.names

huh6.list <- lapply(huh6.list, function(df) {
  df$abs_logFC <- abs(df$logFC)
  return(df)
})
huh6.markers <- lapply(huh6.list, \(x) {
  x <- x[order(x$abs_logFC, decreasing = TRUE), ]
  head(x, 5000)
})
huh6.markers <- lapply(huh6.markers, `[[`, "gene")

marker.names <- paste(names(huh6.markers), ".sig", sep = "")
marker.names <- gsub("_", "-", marker.names)

## Pseudobulk across sample
sample.roehrig<- AggregateExpression(roehrig.seurat, assays = "SCT", return.seurat = T,
                                     group.by = "sample_identity")

sample.roehrig <- AddModuleScore(sample.roehrig, features = huh6.markers, assay = "SCT", seed = 12345, 
                                 name = marker.names)
colnames(sample.roehrig@meta.data) <- gsub("\\.sig[0-9]+$", "\\.sig", colnames(sample.roehrig@meta.data))
sample.roehrig[["Markers"]] <- CreateAssayObject(data = t(FetchData(object = sample.roehrig, vars = marker.names)))

sample.mat <- as.matrix(sample.roehrig@assays$Markers@data)
anno.row <- data.frame("Pseudotime" = gsub("pseudotime[1-9]-", "",
                                           gsub(".sig", "", rownames(sample.mat))),
                       row.names = rownames(sample.mat))
anno.col <- data.frame(Sample = c("Tumour", "Non-tumour"),
                       row.names = colnames(sample.mat))
anno.cols <- list("Sample" = c("Tumour" = "#B1283AFF",
                               "Non-tumour" = "#A8A6A7FF"),
                  "Pseudotime" = c("up" = "#FFB900FF",
                                   "down" = "#5773CCFF"))
gt <- pheatmap(sample.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_row = anno.row,
               annotation_col = anno.col,
               annotation_colors = anno.cols)$gtable
ggsave("plots/roehrig/pheatmap_sample_5000markers.pdf", plot = gt)
ggsave("plots/roehrig/pheatmap_sample_5000markers.png", plot = gt)

## Pseudobulk across cell type
cell.roehrig<- AggregateExpression(roehrig.seurat, assays = "SCT", return.seurat = T,
                                   group.by = "cell_identity")

cell.roehrig <- AddModuleScore(cell.roehrig, features = huh6.markers, assay = "SCT", seed = 12345, 
                               name = marker.names)
colnames(cell.roehrig@meta.data) <- gsub("\\.sig[0-9]+$", "\\.sig", colnames(cell.roehrig@meta.data))
cell.roehrig[["Markers"]] <- CreateAssayObject(data = t(FetchData(object = cell.roehrig, vars = marker.names)))

cell.mat <- as.matrix(cell.roehrig@assays$Markers@data)
anno.row <- data.frame("Pseudotime" = gsub("pseudotime[1-9]-", "",
                                           gsub(".sig", "", rownames(cell.mat))),
                       row.names = rownames(cell.mat))
anno.col <- data.frame(Cell_type = c("Endothelial", "Hepatocytes", "Liver stellate", "Macrophages",
                                     "Other non-tumor", "T cells", "Tumour cells"),
                       row.names = colnames(cell.mat))
pal <- paletteer_d("Redmonder::qMSO12")
anno.cols <- list("Cell_type" = c("Endothelial" = pal[6],
                                  "Hepatocytes" = pal[2],
                                  "Liver stellate" = pal[5],
                                  "Macrophages" = pal[8],
                                  "Other non-tumor" = pal[1],
                                  "T cells" = pal[7],
                                  "Tumour cells" = pal[4]),
                  "Pseudotime" = c("up" = "#FFB900FF",
                                   "down" = "#5773CCFF"))
gt <- pheatmap(cell.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_row = anno.row,
               annotation_col = anno.col,
               annotation_colors = anno.cols)$gtable
ggsave("plots/roehrig/pheatmap_celltype_5000markers.pdf", plot = gt)
ggsave("plots/roehrig/pheatmap_celltype_5000markers.png", plot = gt)

## Pseudobulk across tumour phenotypes
pheno.roehrig<- AggregateExpression(roehrig.seurat, assays = "SCT", return.seurat = T,
                                    group.by = "H_LP_M")

pheno.roehrig <- AddModuleScore(pheno.roehrig, features = huh6.markers, assay = "SCT", seed = 12345, 
                                name = marker.names)
colnames(pheno.roehrig@meta.data) <- gsub("\\.sig[0-9]+$", "\\.sig", colnames(pheno.roehrig@meta.data))
pheno.roehrig[["Markers"]] <- CreateAssayObject(data = t(FetchData(object = pheno.roehrig, vars = marker.names)))

pheno.mat <- as.matrix(pheno.roehrig@assays$Markers@data)
anno.row <- data.frame("Pseudotime" = gsub("pseudotime[1-9]-", "",
                                           gsub(".sig", "", rownames(pheno.mat))),
                       row.names = rownames(pheno.mat))
anno.col <- data.frame(Phenotype = colnames(pheno.mat),
                       row.names = colnames(pheno.mat))
pal <- paletteer::paletteer_d("RColorBrewer::Set2")
anno.cols <- list("Phenotype" = c("H" = "#F4CAE4FF",
                                  "H+LP" = pal[4],
                                  "LP" = pal[3],
                                  "M" = pal[7]),
                  "Pseudotime" = c("up" = "#FFB900FF",
                                   "down" = "#5773CCFF"))
gt <- pheatmap(pheno.mat,
               border_color = NA,
               cellwidth = 6, cellheight = 6,
               fontsize_row = 5, fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_row = anno.row,
               annotation_col = anno.col,
               annotation_colors = anno.cols)$gtable
ggsave("plots/roehrig/pheatmap_phenotype_5000markers.pdf", plot = gt)
ggsave("plots/roehrig/pheatmap_phenotype_5000markers.png", plot = gt)

## Violin plots
filt.seurat <- subset(roehrig.seurat, idents = c("H", "H+LP", "LP", "M"))
filt.seurat <- AddModuleScore(filt.seurat, features = huh6.markers, assay = "SCT", seed = 12345, 
                              name = marker.names)
colnames(filt.seurat@meta.data) <- gsub("\\.sig[0-9]+$", "\\.sig", colnames(filt.seurat@meta.data))

pal <- paletteer::paletteer_d("RColorBrewer::Set2")
roehrig.cols <- c("H" = "#F4CAE4FF",
                  "H+LP" = pal[4],
                  "LP" = pal[3],
                  "M" = pal[7])
Idents(filt.seurat) <- filt.seurat$H_LP_M
VlnPlot(filt.seurat, features = marker.names,
        pt.size = 0, cols = roehrig.cols, ncol = 2) &
  xlab("")
ggsave("plots/roehrig/violin_phenotype_5000markers.pdf", width = 8.7, height = 11.3)
ggsave("plots/roehrig/violin_phenotype_5000markers.png", width = 8.7, height = 11.3)
