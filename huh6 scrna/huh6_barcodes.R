## HuH6 Cellecta Barcode Analysis

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

## Colours for PAM clusters
require(paletteer)
pal <- paletteer_d("rcartocolor::Pastel")
pam.cols <- c("Hepatocytic 1" = pal[2], 
              "Stem-like" = pal[6], 
              "Hepatocytic 3" = pal[3], 
              "Hepatocytic 2" = pal[4], 
              "Late progenitor" = pal[5], 
              "Early progenitor" = pal[1])

## [ Add barcodes ] ----

huh6.seurat <- readRDS("data/huh6_seurat_hvgs_pam_anno.rds")

multi1.df <- read.csv("data/multiplex_1_cellecta_metadata.csv", row.names = 1)
multi2.df <- read.csv("data/multiplex_2_cellecta_metadata.csv", row.names = 1)
multi5.df <- read.csv("data/multiplex_5_cellecta_metadata.csv", row.names = 1)
multi6.df <- read.csv("data/multiplex_6_cellecta_metadata.csv", row.names = 1)

huh6.df <- do.call("rbind", list(multi1.df, multi2.df, multi5.df, multi6.df))
rownames(huh6.df) <- huh6.df$CellID
saveRDS(huh6.df, "data/huh6_cellecta_metadata.rds")

huh6.seurat$CellID <- gsub("-1", "", huh6.seurat$CellID)
colnames(huh6.seurat) <- huh6.seurat$CellID
huh6.filt.seurat <- huh6.seurat[ , colnames(huh6.seurat) %in% rownames(huh6.df)]

test.order <- function(x,y) {
  if (all(x==y)) print('Perfect match in same order')
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}
test.order(huh6.df[colnames(huh6.filt.seurat), ]$CellID, huh6.filt.seurat$CellID)
## "Perfect match in same order"

huh6.filt.seurat$Full.BCS <- huh6.df[colnames(huh6.filt.seurat), ]$Full.BCS
saveRDS(huh6.filt.seurat, "data/huh6_seurat_barcode.rds")

## [ POT clones identity ] ----

huh6.seurat <- readRDS("data/huh6_seurat_barcode.rds")
length(unique(huh6.seurat$Full.BCS)) ## 130

condition.df <- as.data.frame.matrix(table(huh6.seurat$Full.BCS, huh6.seurat$Condition))

pot.df <- condition.df[condition.df$POT >= 1, ] ## 26/130 barcodes in the POT

Idents(huh6.seurat) <- huh6.seurat$Condition
pot.cells <- WhichCells(huh6.seurat, idents = "POT") ## Cells in POT with barcode

DimPlot(huh6.seurat, cells.highlight = pot.cells, group.by = "Full.BCS") +
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

huh6.seurat$PAM_POT <- NA
huh6.seurat$PAM_POT[pot.cells] <- huh6.seurat$PAM_Name[pot.cells]

DimPlot(huh6.seurat, group.by = "PAM_POT", order = TRUE) +
  umap.theme() + labs(title = "PAM Cluster of POT cells") +
  scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(reverse = TRUE, override.aes = list(size = 5)))
ggsave("plots/umap_pam_clusters_pot.pdf", width = 8.3, height = 5.8)
ggsave("plots/umap_pam_clusters_pot.png", width = 8.3, height = 5.8)

pot.data <- huh6.seurat@meta.data[pot.cells, c("Full.BCS", "Condition", "PAM_Name")]
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

huh6.seurat <- readRDS("data/huh6_seurat_barcode.rds")

condition.df <- as.data.frame.matrix(table(huh6.seurat$Full.BCS, huh6.seurat$Condition))
pot.df <- condition.df[condition.df$POT >= 1, ] 
pot.barcodes <- rownames(pot.df)
saveRDS(pot.barcodes, "data/pot_barcodes.rds")

conditions.list <- as.character(unique(huh6.seurat$Condition))
conditions.list <- conditions.list[conditions.list != "POT"] 
cells.list <- list()

## Identify cells with a barcode in the POT
for (condition in 1:length(conditions.list)) {
  cells <- WhichCells(huh6.seurat, cells = rownames(huh6.seurat@meta.data)[
    huh6.seurat@meta.data$Full.BCS %in% pot.barcodes & 
      huh6.seurat@meta.data$Condition == conditions.list[condition]])
  cells.list[[condition]] <- cells
  names(cells.list)[condition] <- conditions.list[condition]
}

barcode.anno <- rownames(pot.df)
barcode.anno <- setNames(1:length(barcode.anno), barcode.anno)
saveRDS(barcode.anno, "data/barcode_anno.rds")

for (i in 1:length(cells.list)) {
  cells <- cells.list[[i]]
  data <- huh6.seurat@meta.data[cells, c("Full.BCS", "Condition", "PAM_Name")]
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
## Grouped-stacked bar chart ///////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

huh6.seurat <- readRDS("data/huh6_seurat_barcode.rds")
pot.barcodes <- readRDS("data/pot_barcodes.rds")
barcode.anno <- readRDS("data/barcode_anno.rds")

filt.seurat <- subset(huh6.seurat, subset = Full.BCS %in% pot.barcodes)
meta.df <- filt.seurat@meta.data[c("Condition", "Full.BCS", "PAM_Name")]
meta.df$Barcode <- recode(as.character(meta.df$Full.BCS), !!!barcode.anno)
saveRDS(meta.df, "data/pot_barcodes_meta_df.rds")

plot.df <- meta.df %>%
  group_by(Barcode, Condition, PAM_Name) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Barcode, Condition)

plot.df$Condition <- recode(as.character(plot.df$Condition), !!!barcode.anno)

ggplot(plot.df) +
  geom_bar(aes(x = Condition, y = Count, fill = PAM_Name),
           stat = "identity",
           position = "stack") +
  facet_grid(~ Barcode, switch = "x") +
  scale_fill_manual(values = pam.cols) +
  scale_x_discrete(labels = conditions.list) +
  ggtitle("Barcode Cluster Proportions Over Time") +
  labs(x = "",
       y = "Cluster proportion",
       fill = "Cluster name") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = "white"),
        panel.spacing = unit(-0.01, "cm"),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 12))
ggsave("plots/barcode_cluster_proportion_condition.pdf", width = 16.5, height = 5.8)
ggsave("plots/barcode_cluster_proportion_condition.png", width = 16.5, height = 5.8)

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
if (max_count == 0) max_count <- 1 ## Prevent division by zero
max_radius <- 0.5  
pie.summary$radius <- ifelse(
  pie.summary$Total_Count > 0, 
  max_radius * sqrt(pie.summary$Total_Count / max_count), 
  0.1  ## Assign small nonzero radius for empty groups
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
x.labels <- c("POT", "Untreated", "Cisplatin", "Recovery t1", "Recovery t2", "Recovery t3")
y.labels <- 1:26
png("plots/barcode_cluster_condition_pie_bubble.png", width = 5.8, height = 16.5,
    units = "in", res = 200)
pdf("plots/barcode_cluster_condition_pie_bubble.pdf", width = 5.8, height = 16.5)
par(mar = c(10, 4, 4, 2))
pie_bubbles(xpos, ypos, radii, sector.list, pam.cols, pam.names, 
            x.labels, y.labels, main = "Pie Bubbles with Proper Axis Labels", ylab = "Barcode")
dev.off()
