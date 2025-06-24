## HepG2 Cellecta Barcode Analysis

## /////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggforce)
library(patchwork)
library(plotrix)

## Labels for experimental conditions
conditions.list <- c("POT", "Untreated", "Cisplatin", bquote("Recovery t"[1]), bquote("Recovery t"[2]))

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

## Colours for PAM clusters
require(paletteer)
pal <- paletteer_d("rcartocolor::Pastel")
pam.cols <- c("Hepatocytic 1" = pal[2], 
              "Stem-like" = pal[6], 
              "Hepatocytic 2" = pal[3], 
              "Early progenitor 1" = pal[1],
              "Late progenitor" = pal[5],
              "Early progenitor 2" = pal[9])

## [ Add barcodes ] ----

hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs_pam_anno.rds")

multi3.df <- read.csv("data/multiplex_3_cellecta_metadata.csv", row.names = 1)
multi4.df <- read.csv("data/multiplex_4_cellecta_metadata.csv", row.names = 1)
multi7.df <- read.csv("data/multiplex_7_cellecta_metadata.csv", row.names = 1)
multi8.df <- read.csv("data/multiplex_8_cellecta_metadata.csv", row.names = 1)

hepg2.df <- do.call("rbind", list(multi3.df, multi4.df, multi7.df, multi8.df))
rownames(hepg2.df) <- hepg2.df$CellID
saveRDS(hepg2.df, "data/hepg2_cellecta_metadata.rds")

hepg2.seurat$CellID <- gsub("-1", "", hepg2.seurat$CellID)
colnames(hepg2.seurat) <- hepg2.seurat$CellID

hepg2.filt.seurat <- hepg2.seurat[ , colnames(hepg2.seurat) %in% rownames(hepg2.df)]

test.order <- function(x,y) {
  if (all(x==y)) print('Perfect match in same order')
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}

test.order(hepg2.df[colnames(hepg2.filt.seurat), ]$CellID, hepg2.filt.seurat$CellID)
## "Perfect match in same order"

hepg2.filt.seurat$Full.BCS <- hepg2.df[colnames(hepg2.filt.seurat), ]$Full.BCS

saveRDS(hepg2.filt.seurat, "data/hepg2_seurat_barcode.rds")

## [ POT clones identity ] ----

hepg2.seurat <- readRDS("data/hepg2_seurat_barcode.rds")
length(unique(hepg2.seurat$Full.BCS)) ## 101

condition.df <- as.data.frame.matrix(table(hepg2.seurat$Full.BCS, hepg2.seurat$Condition))

pot.df <- condition.df[condition.df$POT >= 1, ] ## 36/101 barcodes in the POT

Idents(hepg2.seurat) <- hepg2.seurat$Condition
pot.cells <- WhichCells(hepg2.seurat, idents = "POT") ## Cells in POT with barcode

DimPlot(hepg2.seurat, cells.highlight = pot.cells, group.by = "Full.BCS") +
  scale_colour_manual(values = c("#DDDDDD", "#228833"),
                      labels = c("All other conditions", "POT cells")) +
  labs(title = "Barcodes in the POT", colour = "Barcoded cells") +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "#16161D", linewidth = 0.8),
        axis.ticks = element_line(colour = "#16161D", linewidth = 0.8),
        text = element_text(size = 15)) +
  guides(colour = guide_legend(reverse = TRUE, override.aes = list(size = 5)))
ggsave("plots/umap_pot.pdf", width = 8.3, height = 5.8)
ggsave("plots/umap_pot.png", width = 8.3, height = 5.8)

hepg2.seurat$PAM_POT <- NA
hepg2.seurat$PAM_POT[pot.cells] <- hepg2.seurat$PAM_Name[pot.cells]

DimPlot(hepg2.seurat, group.by = "PAM_POT", order = TRUE) +
  umap.theme() + labs(title = "PAM Cluster of POT cells") +
  scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(reverse = TRUE, override.aes = list(size = 5)))
ggsave("plots/umap_pam_clusters_pot.pdf", width = 8.3, height = 5.8)
ggsave("plots/umap_pam_clusters_pot.png", width = 8.3, height = 5.8)

pot.data <- hepg2.seurat@meta.data[pot.cells, c("Full.BCS", "Condition", "PAM_Name")]
pot.table <- as.data.frame(table(pot.data$Full.BCS, pot.data$PAM_Name))
barcode.anno <- rownames(pot.df)
barcode.anno <- setNames(1:length(barcode.anno), barcode.anno)
pot.table$Name <- recode(as.character(pot.table$Var1), !!!barcode.anno)
pot.table$Name <- as.character(pot.table$Name)
pot.table$Name <- factor(pot.table$Name,
                         levels = c(1:length(barcode.anno)))

ggplot(pot.table, aes(x = Name, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Barcode") +
  ylab("Cluster proportion") +
  ggtitle("POT clones in POT") +
  guides(fill = guide_legend(title = "PAM cluster name")) +
  scale_fill_manual(values = pam.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pam_clusters_pot_proportion.pdf", width = 8.3, height = 5.8)
ggsave("plots/pam_clusters_pot_proportion.png", width = 8.3, height = 5.8)

ggplot(pot.table, aes(x = Name, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("Barcode") +
  ylab("Cluster Frequency") +
  ggtitle("POT clones in POT") +
  guides(fill = guide_legend(title = "PAM cluster name")) +
  scale_fill_manual(values = pam.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pam_clusters_pot_freq.pdf", width = 8.3, height = 5.8)
ggsave("plots/pam_clusters_pot_freq.png", width = 8.3, height = 5.8)

## [ POT clones over time ] ----

hepg2.seurat <- readRDS("data/hepg2_seurat_barcode.rds")

condition.df <- as.data.frame.matrix(table(hepg2.seurat$Full.BCS, hepg2.seurat$Condition))
pot.df <- condition.df[condition.df$POT >= 1, ] 
pot.barcodes <- rownames(pot.df)
saveRDS(pot.barcodes, "data/pot_barcodes.rds")

conditions.list <- as.character(unique(hepg2.seurat$Condition))
conditions.list <- conditions.list[conditions.list != "POT"] 
cells.list <- list()

## Identify cells with a barcode in the POT
for (condition in 1:length(conditions.list)) {
  cells <- WhichCells(hepg2.seurat, cells = rownames(hepg2.seurat@meta.data)[
    hepg2.seurat@meta.data$Full.BCS %in% pot.barcodes & 
      hepg2.seurat@meta.data$Condition == conditions.list[condition]])
  cells.list[[condition]] <- cells
  names(cells.list)[condition] <- conditions.list[condition]
}

barcode.anno <- rownames(pot.df)
barcode.anno <- setNames(1:length(barcode.anno), barcode.anno)
saveRDS(barcode.anno, "data/barcode_anno.rds")

for (i in 1:length(cells.list)) {
  cells <- cells.list[[i]]
  data <- hepg2.seurat@meta.data[cells, c("Full.BCS", "Condition", "PAM_Name")]
  table <- as.data.frame(table(data$Full.BCS, data$PAM_Name))
  table$Name <- recode(as.character(table$Var1), !!!barcode.anno)
  table$Name <- as.character(table$Name)
  table$Name <- factor(table$Name,
                       levels = c(1:length(barcode.anno)))
  
  file.name <- paste0("plots/pam_clusters_", gsub("_", "",
                                                  gsub("cisplatin", "", tolower(names(cells.list)[i]))),"_proportion")
  ggplot(table, aes(x = Name, y = Freq, fill = Var2)) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Barcode") +
    ylab("Cluster proportion") +
    ggtitle(paste0("POT clones in ", gsub("_", " ", tolower(names(cells.list)[i])))) +
    guides(fill = guide_legend(title = "PAM cluster name")) +
    scale_fill_manual(values = pam.cols) +
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(file.name, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0(file.name, ".png"), width = 8.3, height = 5.8)
  
  file.name <- paste0("plots/pam_clusters_", gsub("_", "",
                                                  gsub("cisplatin", "",tolower(names(cells.list)[i]))),"_freq")
  ggplot(table, aes(x = Name, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", width = 0.5) +
    xlab("Barcode") +
    ylab("Cluster Frequency") +
    ggtitle(paste0("POT clones in ", gsub("_", " ", tolower(names(cells.list)[i])))) +
    guides(fill = guide_legend(title = "PAM cluster name")) +
    scale_fill_manual(values = pam.cols) +
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(file.name, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0(file.name, ".png"), width = 8.3, height = 5.8)
}

## /////////////////////////////////////////////////////////////////////////////
## Bubble pie chart ////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

meta.df <- readRDS("data/pot_barcodes_meta_df.rds")
pie.df <- meta.df %>%
  group_by(Condition, Barcode, PAM_Name) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Condition, Barcode, PAM_Name, fill = list(Count = 0)) %>%
  group_by(Condition, Barcode) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(x = as.numeric(factor(Condition)),
         y = as.numeric(factor(Barcode)))
pie.df[is.na(pie.df)] <- 0

## Compute total counts per condition and barcode
pie.summary <- pie.df %>%
  group_by(Condition, Barcode) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), 
            x_coord = first(x), 
            y_coord = first(y), 
            .groups = 'drop')

## Calculate radii to normalise bubble size
max_count <- max(pie.summary$Total_Count, na.rm = TRUE)
# if (max_count == 0) max_count <- 1 ## Prevent division by zero
max_radius <- 0.3  
pie.summary$radius <- ifelse(
  pie.summary$Total_Count > 0, 
  max_radius * sqrt(pie.summary$Total_Count / max_count), 
  0.001  ## Assign small nonzero radius for empty groups
)

## Plotting function for bubble pie chart
pie_bubbles <- function(xpos, ypos, radii, sectors, pam.cols, pam.names, 
                        x_labels, y_labels, main = "", xlab = "", ylab = "") {
  
  xlim <- c(min(xpos - radii, na.rm = TRUE), max(xpos + radii, na.rm = TRUE))
  ylim <- c(min(ypos - radii, na.rm = TRUE), max(ypos + radii, na.rm = TRUE))
  
  nbubbles <- length(xpos)
  sector_col <- list()
  
  ## Assign sector colors based on cluster
  for (bubble in 1:nbubbles) {
    sector_col[[bubble]] <- pam.cols[pam.names[[bubble]]]
  }
  
  ## Create empty plot
  plot(0, xlim = xlim, ylim = ylim, type = "n", main = main, xlab = xlab, ylab = ylab, axes = FALSE)
  
  ## Draw pie charts
  for (bubble in 1:nbubbles) {
    if (length(sectors[[bubble]]) > 0 && sum(sectors[[bubble]], na.rm = TRUE) > 0) {
      floating.pie(
        xpos = xpos[bubble], ypos = ypos[bubble], 
        x = sectors[[bubble]] / sum(sectors[[bubble]], na.rm = TRUE),  # Normalize proportions
        radius = radii[bubble], 
        col = sector_col[[bubble]]
      )
    }
  }
  
  # Add y-axis labels (1 to 26 for Barcode)
  axis(2, at = seq_along(y_labels), labels = y_labels, las = 1)
  
  # Add x-axis labels (Condition names)
  axis(1, at = seq_along(x_labels), labels = FALSE)  # Hide default labels
  text(x = seq_along(x_labels), y = par("usr")[3] - 1.5, labels = x_labels, 
       srt = 45, adj = 1, xpd = TRUE)
}

## Extract variables
xpos <- pie.summary$x_coord
ypos <- pie.summary$y_coord
radii <- pie.summary$radius

sector.list <- split(pie.df$Proportion, list(pie.df$Condition, pie.df$Barcode))
order <- paste(pie.summary$Condition, pie.summary$Barcode, sep = ".") ## Expected order
sector.list <- sector.list[order]

pam.names <- split(pie.df$PAM_Name, list(pie.df$Condition, pie.df$Barcode))
pam.names <- pam.names[order]
x.labels <- c("POT", "Untreated", "Cisplatin", "Recovery t1", "Recovery t2")
y.labels <- 1:36
png("plots/barcode_cluster_condition_pie_bubble.png", width = 5.8, height = 16.5,
    units = "in", res = 200)
pdf("plots/barcode_cluster_condition_pie_bubble.pdf", width = 5.8, height = 16.5)
par(mar = c(10, 4, 4, 2))
pie_bubbles(xpos, ypos, radii, sector.list, pam.cols, pam.names, 
            x.labels, y.labels, main = "Pie Bubbles with Proper Axis Labels", ylab = "Barcode")
dev.off()

## Make pie-bubbles for defined Total_Count values to use for the legend
head(pie.summary)
pie.summary$Total_Count[1] = 1
pie.summary$Total_Count[2] = 10
pie.summary$Total_Count[3] = 100
pie.summary$Total_Count[4] = 1000
head(pie.summary)
## Recalculate radii
max_count <- max(pie.summary$Total_Count, na.rm = TRUE)
if (max_count == 0) max_count <- 1
max_radius <- 0.3  
pie.summary$radius <- ifelse(
  pie.summary$Total_Count > 0, 
  max_radius * sqrt(pie.summary$Total_Count / max_count), 
  0.001
)
## Extract variables
xpos <- pie.summary$x_coord
ypos <- pie.summary$y_coord
radii <- pie.summary$radius
sector.list <- split(pie.df$Proportion, list(pie.df$Condition, pie.df$Barcode))
order <- paste(pie.summary$Condition, pie.summary$Barcode, sep = ".") ## Expected order
sector.list <- sector.list[order]
pam.names <- split(pie.df$PAM_Name, list(pie.df$Condition, pie.df$Barcode))
pam.names <- pam.names[order]
x.labels <- c("POT", "Untreated", "Cisplatin", "Recovery t1", "Recovery t2")
y.labels <- 1:36
## Make palette white
pam.cols <- c("Hepatocytic 1" = "white", 
              "Stem-like" = "white", 
              "Hepatocytic 3" = "white", 
              "Hepatocytic 2" = "white", 
              "Late progenitor" = "white", 
              "Early progenitor" = "white")
pdf("plots/barcode_pie_bubble_legend.pdf", width = 5.8, height = 16.5)
par(mar = c(10, 4, 4, 2))
pie_bubbles(xpos, ypos, radii, sector.list, pam.cols, pam.names, 
            x.labels, y.labels, main = "Pie Bubbles with Proper Axis Labels", ylab = "Barcode")
dev.off()
