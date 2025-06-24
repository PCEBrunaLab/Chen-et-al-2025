## HuH6 Cellecta barcodes DNA analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(tidyverse)
library(ggplot2)
library(vegan)
library(Polychrome)
library(viridis)
library(patchwork)
library(cowplot)

## Replicate colours
rep.cols <- c("Cisplatin A" = "#dfd48a", "Cisplatin B" = "#CCBB44", "Cisplatin C" = "#80752b",
              "UT A" = "#8aaaca", "UT B" = "#4477AA", "UT C" = "#2b4a6a",
              "Cisplatin R 1week A" = "#f49faa", "Cisplatin R 1week B" = "#EE6677", "Cisplatin R 1week C" = "#95404a",
              "Cisplatin R 3weeks A" = "#ca80aa", "Cisplatin R 3weeks B" = "#aa3377", "Cisplatin R 3weeks C" = "#6a204a",
              "Cisplatin R 3months A" = "#9fdff4", "Cisplatin R 3months B" = "#66CCEE", "Cisplatin R 3months C" = "#408095")

## Replicate plot labels
rep.list <- c("Untreated A", "Untreated B", "Untreated C",
              "Cisplatin A", "Cisplatin B", "Cisplatin C",
              expression("Recovery t"[1]*" A"), expression("Recovery t"[1]*" B"), expression("Recovery t"[1]*" C"),
              expression("Recovery t"[2]*" A"), expression("Recovery t"[2]*" B"), expression("Recovery t"[2]*" C"),
              expression("Recovery t"[3]*" A"), expression("Recovery t"[3]*" B"), expression("Recovery t"[3]*" C"))

## [ Data prep ] ----

all.barcodes <- read.csv("data/all_barcodes_HB.tsv")

huh6.barcodes <- all.barcodes[grep("HUH6", all.barcodes$sample_id), ]

length(unique(huh6.barcodes$real_bc44)) ## 102171 unique barcodes

meta.df <- read.csv("data/HB_scRNA_samples.csv")
huh6.barcodes <- merge(huh6.barcodes, meta.df, by = c("sample_id"))

write.csv(huh6.barcodes, "data/HuH6_barcodes.csv", row.names = FALSE)

## [ Initial processing ] ----

huh6.barcodes <- read.csv("data/HuH6_barcodes.csv")

## Change proportion to percentage
huh6.barcodes$proportion <- huh6.barcodes$perc
huh6.barcodes$perc <- NULL
huh6.barcodes <- huh6.barcodes %>%
  mutate(perc = proportion * 100)

## Define full_sample variable
huh6.barcodes <- huh6.barcodes %>%
  mutate(full_sample = ifelse(recovery == TRUE, paste0(condition, " ", "R", " ", recovery_timepoint, " ", replicate), 
                              paste0(condition, " ", replicate)))
huh6.barcodes <- huh6.barcodes %>%
  mutate(sample = ifelse(recovery == TRUE, paste0(condition, " ", "R", " ", recovery_timepoint), paste0(condition)))

saveRDS(huh6.barcodes, "data/huh6_barcodes_unfiltered.rds")

## Summarise the number of unique barcodes within each sample
barcode.summary <- huh6.barcodes %>%
  select(c(real_bc44, full_sample, sample)) %>%
  group_by(full_sample, sample) %>%
  summarise(count = n_distinct(real_bc44))
barcode.summary$sample <- factor(barcode.summary$sample, 
                                 levels = c("UT", "Cisplatin", "Cisplatin R 1week",
                                            "Cisplatin R 3weeks", "Cisplatin R 3months"))

## Plot number of unique barcodes observed across samples
options(scipen = 999)

ggplot(data = barcode.summary, aes(x = sample, y = count)) +
  geom_boxplot() +
  ggtitle("Clone Frequency Across Experimental\nSamples") + 
  xlab("Condition")+
  ylab("Number of Unique Cellecta Clones")+
  expand_limits(y=c(0, 40000)) +
  scale_y_continuous(breaks = seq(0, 40000, 5000)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15))

ggsave("plots/HuH6_clones_per_condition.pdf", width = 5.8, height = 5.8)
ggsave("plots/HuH6_clones_per_condition.png", width = 5.8, height = 5.8)

## [ Initial analysis ] ----

## Calculate Shannon diversity
shannon.summary <- huh6.barcodes %>%
  group_by(full_sample) %>%
  mutate(shannon = (-1)* (sum(proportion * log(proportion))),
         richness = n_distinct(real_bc44)) %>%
  ungroup() %>%
  select(c(shannon, sample, richness)) %>%
  distinct()

## Visualise Shannon diversity
shannon.summary$sample <- factor(shannon.summary$sample, 
                                 levels = c("UT", "Cisplatin", "Cisplatin R 1week",
                                            "Cisplatin R 3weeks", "Cisplatin R 3months"))

## Plot Shannon diveristy scores
ggplot(shannon.summary, aes(x = sample, y = shannon)) +
  geom_bar(stat = "summary", fill = "white", colour = "black") +
  geom_point()+
  xlab("Condition") +
  ylab("Shannon Diversity Index") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15))

ggsave("plots/HuH6_shannon_diversity_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/HuH6_shannon_diversity_bar.png", width = 5.8, height = 5.8)

## Cumulative frequency plots for barcode diversity
group.cols <- c("POT" = "#228833", 
                "UT" = "#4477AA", 
                "Cisplatin" = "#EE6677", 
                "Cisplatin R 1week" = "#CCBB44", 
                "Cisplatin R 3weeks" = "#66CCEE", 
                "Cisplatin R 3months" = "#AA3377")

huh6.barcodes %>%
  group_by(full_sample, real_bc44, sample) %>%
  summarise(abundance = sum(N)) %>%
  ungroup() %>%
  group_by(full_sample, sample) %>%
  arrange(desc(abundance)) %>%
  mutate(cumulative_prop = cumsum(abundance) / sum(abundance),
         id = row_number()) %>%
  ungroup() %>%
  group_by(sample, id) %>%
  summarise(min = min(cumulative_prop),
            max = max(cumulative_prop),
            avg = mean(cumulative_prop)) %>%
  ungroup() %>%
  ggplot(aes(x = id, colour = sample))+
  geom_ribbon(aes(ymin = min, ymax = max, fill = sample), alpha = 0.5, lty = 0) +
  geom_line(aes(y = avg))+
  scale_colour_manual(values = group.cols, name = "Condition") +
  scale_fill_manual(values = group.cols, name = "Condition") +
  scale_x_log10() +
  xlab("Unique barcode cumulative frequency") +
  ylab("Average barcode abundance cumulative proportion") +
  theme_bw() + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15))

ggsave("plots/HuH6_cumulative_freq_selection.pdf", width = 5.8, height = 5.8)
ggsave("plots/HuH6_cumulative_freq_selection.png", width = 5.8, height = 5.8)

## Plot the top 2000 barcodes from each sample
## We don't have a POT in this case so select all the UT samples to use as a "POT"
ut.barcodes <- unique(huh6.barcodes$real_bc44[huh6.barcodes$condition == "UT"]) ## 34160 clones
filtered.barcode.df <- huh6.barcodes %>%
  filter(real_bc44 %in% ut.barcodes)

## Define colour palette for plotting 2000 unique clones
set.seed(12)
p2000 = createPalette(2000, c("#ff0000", "#00ff00", "#0000ff"), M = 2000)
p2000 <- as.vector(t(matrix(p2000)))
ut.df <- huh6.barcodes %>%
  filter(condition == "UT")
top2000.ut <- ut.df %>%
  arrange(desc(perc)) %>%
  dplyr::slice(1:2000)
names(p2000) = unique(top2000.ut$real_bc44)

## Plot top 2000 barcodes for each sample
filtered.barcode.df %>%
  mutate(sample = factor(sample, levels = c("UT", "Cisplatin", "Cisplatin R 1week",
                                            "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  ggplot(aes(x = replicate, y = perc, fill = reorder(real_bc44, desc(perc)))) +
  geom_bar(stat = "summary", position = "stack") +
  facet_wrap(~ sample, ncol = 5, scales = "free", labeller = label_wrap_gen(width = 10))+
  scale_fill_manual(values = p2000) +
  ylab("Cellecta Barcode (%)")+
  xlab("")+
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = "none")

ggsave("plots/HuH6_barcode_proportions.pdf", width = 8.7, height = 5.8)
ggsave("plots/HuH6_barcode_proportions.png", width = 8.7, height = 5.8)

## Calculate percentage reduction in each clonal population for each sample compared to POT
reduction.df <- filtered.barcode.df %>%
  group_by(full_sample, sample) %>%
  summarise(barcodes = n_distinct(real_bc44)) %>%
  ungroup() %>%
  mutate(representation = barcodes/barcodes[full_sample == "UT A" | full_sample == "UT B" | full_sample == "UT C"] * 100,
         reduction = 100 - representation)

## Plot reduction in clonal population
reduction.df %>%
  mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                          "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  ggplot(aes(x = sample, y = reduction)) +
  geom_boxplot(colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Experimental Sample")+
  ylab("Cellecta Clonal Reduction (%)")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15))

ggsave("plots/HuH6_reduction_boxplot.png", width = 5.8, height = 5.8)
ggsave("plots/HuH6_reduction_boxplot.png", width = 5.8, height = 5.8)

## Because the entire sample was library prepped proceed without estimating cell frequency
## Find top 75 clones in each sample
top.75 <- filtered.barcode.df %>%
  filter(condition == "UT") %>%
  arrange(desc(perc)) %>%
  distinct(real_bc44) %>%
  head(75)

## Identify clones from POT sample and order these clones based on representation in POT sample
ut.perc.df <- filtered.barcode.df %>%
  filter(condition == "UT") %>%
  arrange(desc(perc)) %>%
  distinct(real_bc44)
ut.perc.barcodes <- unique(ut.perc.df$real_bc44)

## Bubble plot for top 75
filtered.barcode.df %>%
  filter(real_bc44 %in% top.75$real_bc44) %>%
  mutate(sample = factor(sample, levels = c("UT", "Cisplatin", "Cisplatin R 1week",
                                            "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels = ut.perc.barcodes)) %>%
  ggplot( aes(x=replicate, y=desc(real_bc44), size = N, colour=real_bc44)) +
  geom_point(alpha = 0.7) +
  facet_grid(~ sample, scales = "free", labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(values = p2000)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 6), breaks = seq(0, 30000, 5000)) +
  guides(colour = FALSE, size = guide_legend(ncol = 3,
                                             title = "Clone\nFrequency",
                                             title.position = "left")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.45, -0.15),
        plot.margin = unit(c(1, 1, 3, 1), "cm"))

ggsave("plots/HuH6_top75_bubble_plot.pdf", width = 5.8, height = 8.7)
ggsave("plots/HuH6_top75_bubble_plot.png", width = 5.8, height = 8.7)

## Select subset of clones to look at more closely
top.15 <- filtered.barcode.df %>%
  filter(condition == "UT") %>%
  arrange(desc(perc)) %>%
  distinct(real_bc44) %>%
  head(15)

## Bubble plot for top 15
filtered.barcode.df %>%
  filter(real_bc44 %in% top.15$real_bc44) %>%
  mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                          "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels = ut.perc.barcodes)) %>%
  ggplot(aes(x = replicate, y = desc(real_bc44), size = N, colour = real_bc44)) +
  geom_point(alpha = 0.7) +
  facet_grid(~ sample, scales = "free", labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(values = p2000)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 6), breaks = seq(0, 30000, 5000)) +
  guides(colour = FALSE, size = guide_legend(ncol = 3,
                                             title = "Clone\nFrequency",
                                             title.position = "left")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.45, -0.15),
        plot.margin = unit(c(1, 1, 3, 1), "cm"))

ggsave("plots/HuH6_top15_bubble_plot.pdf", width = 5.8, height = 8.7)
ggsave("plots/HuH6_top15_bubble_plot.png", width = 5.8, height = 8.7)

## Isolate and plot clones of interest to visualise change of population overtime
bubble.filtered.df <- filtered.barcode.df %>% 
  filter(real_bc44 %in% top.15$real_bc44)

bubble.filtered.df$timepoint[bubble.filtered.df$sample == "UT"] <- "1"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "Cisplatin"] <- "2"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "Cisplatin R 1week"] <- "3"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "Cisplatin R 3weeks"] <- "4"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "Cisplatin R 3months"] <- "5"

bubble.filtered.df$timepoint <- as.numeric(bubble.filtered.df$timepoint)

## Cisplatin clonal dynamics plots 
ggplot(bubble.filtered.df %>%
         mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                                 "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
         mutate(real_bc44 = factor(real_bc44, levels = ut.perc.barcodes)),
       aes(x = sample, y = N, colour = real_bc44, group = replicate)) +
  geom_line(aes(x = sample, y = N), size = 1.25) +
  geom_point(size = 2)+
  scale_y_continuous(limits = c(0, 30000)) +
  scale_colour_manual(values = p2000)+
  xlab("") +
  facet_wrap(~ real_bc44, nrow = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15),
        strip.text.x = element_blank(),
        legend.position = "none")

ggsave("plots/HuH6_cisplatin_dynamics.pdf", width = 8.7, height = 5.8)
ggsave("plots/HuH6_cisplatin_dynamics.png", width = 8.7, height = 5.8)

ggplot(bubble.filtered.df %>%
         mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                                 "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
         mutate(real_bc44 = factor(real_bc44, levels = ut.perc.barcodes)),
       aes(x = sample, y = log(N), colour = real_bc44, group = replicate)) +
  geom_line(aes(x = sample, y = log(N)), size = 1.25) +
  geom_point(size = 2)+
  scale_y_continuous(limits = c(0, 12)) +
  scale_colour_manual(values = p2000)+
  facet_wrap(~real_bc44, nrow = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 15),
        strip.text.x = element_blank(),
        legend.position = "none")

ggsave("plots/HuH6_cisplatin_dynamics_log.pdf", width = 8.7, height = 5.8)
ggsave("plots/HuH6_cisplatin_dynamics_log.png", width = 8.7, height = 5.8)

## Isolate Recovery 3 barcodes with an abundance of > 0.1%
rec3.df <- filtered.barcode.df %>%
  mutate(sample = factor(sample, levels = c("UT", "Cisplatin", "Cisplatin R 1week",
                                            "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  filter(sample == "Cisplatin R 3months", perc > 0.1) %>%
  arrange(replicate, desc(perc))
rec3.barcodes <- unique(rec3.df$real_bc44)

rec3.cols <- hues::iwanthue(length(rec3.barcodes))
names(rec3.cols) <- rec3.barcodes

filtered.barcode.df %>%
  filter(real_bc44 %in% rec3.barcodes) %>%
  mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                          "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels = rec3.barcodes)) %>%
  ggplot(aes(x = replicate, y = desc(real_bc44), size = N, colour = real_bc44)) +
  geom_point(alpha = 0.7) +
  facet_grid(~ sample, scales = "free", labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(values = rec3.cols)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 6)) +
  guides(colour = FALSE, size = guide_legend(ncol = 2,
                                             title = "Clone\nFrequency",
                                             title.position = "left")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.45, -0.15),
        plot.margin = unit(c(1, 1, 3, 1), "cm"))

ggsave("plots/HuH6_rec_3_bubble_plot.pdf", width = 5.8, height = 8.7)
ggsave("plots/HuH6_rec_3_bubble_plot.png", width = 5.8, height = 8.7)

filtered.barcode.df %>%
  filter(real_bc44 %in% rec3.barcodes) %>%
  mutate(sample = factor(sample, levels=c("UT", "Cisplatin", "Cisplatin R 1week",
                                          "Cisplatin R 3weeks", "Cisplatin R 3months"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels = rec3.barcodes)) %>%
  ggplot(aes(x = replicate, y = desc(real_bc44), size = log(N), colour = real_bc44)) +
  geom_point(alpha = 0.7) +
  facet_grid(~ sample, scales = "free", labeller = label_wrap_gen(width = 10)) +
  scale_colour_manual(values = rec3.cols)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 6)) +
  guides(colour = FALSE, size = guide_legend(ncol = 4,
                                             title = "Log(Clone Frequency)",
                                             title.position = "left")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        legend.position = c(0.45, -0.15),
        plot.margin = unit(c(1, 1, 3, 1), "cm"))

ggsave("plots/HuH6_rec_3_bubble_plot_log.pdf", width = 5.8, height = 8.7)
ggsave("plots/HuH6_rec_3_bubble_plot_log.png", width = 5.8, height = 8.7)

## [ Diversity plots ] ----

huh6.barcodes <- readRDS("data/huh6_barcodes_unfiltered.rds")

huh6.filt <- huh6.barcodes %>%
  select(real_bc44, full_sample, N)

reps <- unique(huh6.filt$full_sample)
barcodes <- unique(huh6.filt$real_bc44)

huh6.mat <- data.frame(matrix(nrow = length(reps), ncol = length(barcodes)))
row.names(huh6.mat) <- reps
colnames(huh6.mat) <- barcodes

for (i in 1:nrow(huh6.filt)) {
  replicate <- huh6.filt$full_sample[i]
  barcode <- huh6.filt$real_bc44[i]
  huh6.mat[replicate, barcode] <- huh6.filt$N[i]
}

huh6.mat[is.na(huh6.mat)] <- 0

saveRDS(huh6.mat, "data/huh6_barcode_mat.rds")

huh6.meta <- unique(huh6.barcodes %>%
                      select(full_sample, recovery, recovery_timepoint) %>%
                      group_by(full_sample))

saveRDS(huh6.meta, "data/huh6_barcode_meta.rds")

## Load data
huh6.data <- readRDS("data/huh6_barcode_mat.rds")
huh6.meta <- readRDS("data/huh6_barcode_meta.rds")

## /////////////////////////////////////////////////////////////////////////////
## Alpha diversity /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Calculate alpha diversity
H <- diversity(huh6.data) ## Shannon's H'

richness <- specnumber(huh6.data) ## Observed richness

evenness <- H/log(richness) ## Pielou's evenness

alpha <- cbind(shannon = H, richness = richness, pielou = evenness, huh6.meta)
head(alpha)

## Plot alpha diversity
plot.shan <- ggplot(alpha, aes(x = full_sample, y = shannon, colour = full_sample)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rep.cols) +
  ylab("Shannon's H'") + 
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

plot.rich <-ggplot(alpha, aes(x = full_sample, y = richness, colour = full_sample)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rep.cols) +
  ylab("Species Richness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

plot.even <- ggplot(alpha, aes(x = full_sample, y = evenness, colour = full_sample)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rep.cols) +
  ylab("Pielou's Evenness") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

legend <- get_legend(plot.even)

plot_grid(plot.shan + theme(legend.position = "none"),
          plot.rich + theme(legend.position = "none"),
          plot.even + theme(legend.position = "none"),ncol = 3)

ggsave("plots/huh6_alpha_diversity.pdf", width = 8.7, height = 5.8)
ggsave("plots/huh6_alpha_diversity.png", width = 8.7, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Beta diversity //////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Calculate pairwise dissimilarity
huh6.mdf <- as.matrix.data.frame(huh6.data)
rownames(huh6.mdf) <- huh6.meta$full_sample

huh6.bray <- vegdist(huh6.mdf, method = "bray")
huh6.bray

huh6.jac <- vegdist(huh6.mdf, method = "jaccard", binary = T)
huh6.jac

#simper(huh6.data, huh6.meta$full_sample, permutations = 999)

## Ordination
## Principal Coordinates Analysis (PCoA)

## Calculate principal coordinates analysis (Bray-Curtis)
pcoa.huh6.bray <- cmdscale(huh6.bray, k = 2, eig = T)

pcoa.huh6.bray.plotting <- as.data.frame(pcoa.huh6.bray$points)
colnames(pcoa.huh6.bray.plotting) <- c("axis_1", "axis_2")
pcoa.huh6.bray.plotting$site <- rownames(pcoa.huh6.bray.plotting)

pcoa.huh6.bray$eig[1]/sum(pcoa.huh6.bray$eig) ## 0.4143492

pcoa.huh6.bray$eig[2]/sum(pcoa.huh6.bray$eig) ## 0.1473642

pcoa.huh6.bray.plot <- ggplot(pcoa.huh6.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rep.cols) +
  theme_bw() +
  xlab("PCoA 1 (41.4%)") +
  ylab("PCoA 2 (14.7%)") +
  annotate(geom = "text", label = "Bray-Curtis", x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

## Repeat process with Jaccard dissimilarity matrix
pcoa.huh6.jac <- cmdscale(huh6.jac, k = 2, eig = T)

pcoa.huh6.jac.plotting <- as.data.frame(pcoa.huh6.jac$points)
colnames(pcoa.huh6.jac.plotting) <- c("axis_1", "axis_2")
pcoa.huh6.jac.plotting$site <- rownames(pcoa.huh6.jac.plotting)

pcoa.huh6.jac$eig[1]/(sum(pcoa.huh6.jac$eig)) ## 0.2625131

pcoa.huh6.jac$eig[2]/(sum(pcoa.huh6.jac$eig)) ## 0.1087771

pcoa.huh6.jac.plot <- ggplot(pcoa.huh6.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_manual(values = rep.cols) +
  theme_bw() + 
  xlab("PCoA 1 (26.3%)") +
  ylab("PCoA 2 (10.9%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

legend <- get_legend(pcoa.huh6.jac.plot)

plot_grid(pcoa.huh6.bray.plot + theme(legend.position = 'none'),
          pcoa.huh6.jac.plot + theme(legend.position = 'none'),
          legend, ncol = 3, rel_widths = c(1,1,0.5))

ggsave("plots/huh6_dissimilarity.pdf", width = 8.7, height = 5.8)
ggsave("plots/huh6_dissimilarity.png", width = 8.7, height = 5.8)

## [ Alluvial plot ] ----

huh6.df <- readRDS("data/huh6_barcodes_unfiltered.rds")
ut.barcodes <- unique(huh6.df$real_bc44[huh6.df$condition == "UT"]) ## 34160 clones
anno.df <- data.frame(Barcode = ut.barcodes,
                      Name = paste0("Barcode", seq_along(ut.barcodes)),
                      stringsAsFactors = FALSE)
filt.df <- huh6.df %>%
  filter(real_bc44 %in% ut.barcodes) %>%
  select(real_bc44, condition, replicate, N, full_sample, sample)
filt.df <- filt.df %>%
  left_join(anno.df, by = c("real_bc44" = "Barcode"))
saveRDS(filt.df, "data/huh6_barcodes_filtered.rds")

rep.list <- split(filt.df, filt.df$replicate)

add_percentages <- function(df) {
  df %>%
    group_by(sample) %>%
    mutate(Percent = 100 * N / sum(N)) %>%
    ungroup()
}
rep.list <- map(rep.list, add_percentages)

process_replicate <- function(df) {
  top.barcodes <- df %>%
    group_by(sample) %>%
    slice_max(order_by = Percent, n = 500, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(sample, Name)
  
  df <- df %>%
    left_join(top.barcodes %>% mutate(top = TRUE), by = c("sample", "Name")) %>%
    mutate(top = ifelse(is.na(top), FALSE, top))
  
  df.top <- df %>% filter(top) %>% select(-top)
  df.other <- df %>%
    filter(!top) %>%
    group_by(sample, condition, replicate) %>%
    summarise(Percent = sum(Percent), .groups = "drop") %>%
    mutate(Name = "Other")
  
  bind_rows(df.top, df.other)
}
process.list <- map(rep.list, process_replicate)
process.list <- map(process.list, ~ .x %>%
                      mutate(sample = factor(sample, levels = c("UT", "Cisplatin",
                                                                "Cisplatin R 1week",
                                                                "Cisplatin R 3weeks",
                                                                "Cisplatin R 3months"))))

# unique.barcodes <- process.list %>%
#   map(~ pull(.x, Name)) %>%
#   unlist() %>%
#   unique()
# barcode.pal <- Polychrome::createPalette(2265, c("#ff0000", "#00ff00", "#0000ff"), M = 2265)
# barcode.pal <- as.vector(t(matrix(barcode.pal)))
# names(barcode.pal) <- unique.barcodes
# barcode.pal["Other"] <- "#DCDCDC"
# scales::show_col(barcode.pal[c("Barcode6423", "Barcode14623", "Barcode1662", "Barcode9640", "Barcode3297", "Barcode5495")])
# saveRDS(barcode.pal, "data/barcode_palette.rds")
barcode.pal <- readRDS("data/barcode_palette.rds")

for (df in 1:length(process.list)) {
  name <- names(process.list)[df]
  plot.df <- process.list[[df]]
  
  barcode.allu <- ggplot(plot.df, 
                         aes(x = sample, 
                             stratum = Name, 
                             alluvium = Name,           
                             y = Percent,
                             fill = Name, 
                             label = Name)) + 
    scale_x_discrete(expand = c(0,0)) + 
    geom_flow() +
    geom_stratum(colour = NA)
  ord.levels <- levels(barcode.allu$Name)
  ord.levels <- as.numeric(gsub("Other", "999", gsub("Barcode", "", ord.levels)))
  
  column_heights <- plot.df %>%
    group_by(sample) %>%
    summarise(total = sum(Percent))
  
  print(barcode.allu + 
          geom_rect(data = column_heights,
                    aes(xmin = as.numeric(factor(sample)) - 0.17,
                        xmax = as.numeric(factor(sample)) + 0.17,
                        ymin = 0, ymax = total),
                    fill = NA, color = "black", size = 0.6, inherit.aes = FALSE) +
          scale_fill_manual(values = c(barcode.pal)) +
          scale_x_discrete(expand = c(0.1, 0)) +
          xlab("") +
          ylab("Percentage Abundance") +
          ggtitle(paste0("Replicate ", name)) +
          theme_bw() + 
          theme(panel.grid = element_blank(),
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1),
                text = element_text(size = 12)))
  ggsave(paste0("plots/figures/huh6_alluvial_rep", name, ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/figures/huh6_alluvial_rep", name, ".png"), width = 8.3, height = 5.8)
}

## [ Common recovery barcodes ] ----

filt.df <- readRDS("data/huh6_barcodes_filtered.rds")
add_percentages <- function(df) {
  df %>%
    group_by(sample) %>%
    mutate(Percent = 100 * N / sum(N)) %>%
    ungroup()
}
filt.df <- add_percentages(filt.df)
filt.df <- filter(filt.df, sample == "Cisplatin R 3months")

## How many barcodes are shared between replicates by Recovery t3
rep.list <- split(filt.df$Name, filt.df$replicate)
all.barcodes <- unique(filt.df$Name)
overlap.matrix <- data.frame(Barcode = all.barcodes,
                             A = as.integer(all.barcodes %in% rep.list[[1]]),
                             B = as.integer(all.barcodes %in% rep.list[[2]]),
                             C = as.integer(all.barcodes %in% rep.list[[3]]))
colnames(overlap.matrix)[-1] <- names(rep.list)

png("plots/huh6_shared_t3_barcodes_upset.png", units = "in", width = 8.7, height = 5.8, res = 200)
pdf("plots/huh6_shared_t3_barcodes_upset.pdf", width = 8.7, height = 5.8)
upset(overlap.matrix,
      sets = names(rep.list),
      order.by = "freq",   
      mainbar.y.label = "Barcode Replicate Intersections",    
      sets.x.label = "Number of Barcodes",   
      text.scale = 1.5,   
      keep.order = TRUE)
dev.off()

## Add overlap information back into dataframe
overlap.long <- overlap.matrix %>%
  pivot_longer(cols = -Barcode, names_to = "Group", values_to = "Present") %>%
  filter(Present == 1) %>%
  group_by(Barcode) %>%
  summarise(Overlap = paste0(sort(Group), collapse = ""), .groups = "drop")

filt.df <- filt.df %>%
  left_join(overlap.long, by = c("Name" = "Barcode"))

ggplot(abc.df, aes(x = replicate, y = Percent, colour = replicate)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  scale_y_log10() +
  facet_wrap(~Name, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("A" = "#9fdff4", "B" = "#66CCEE", "C" = "#408095")) +
  labs(title = "Barcode Percent by Replicate",
       x = "Replicate",
       y = "Percent (log10 scale)") +
  theme_bw() +
  theme(text = element_text(size = 12))
ggsave("plots/huh6_shared_t3_barcodes_percent_point_facet.png", width = 5.8, height = 16.5)
ggsave("plots/huh6_shared_t3_barcodes_percent_point_facet.pdf", width = 5.8, height = 16.5)
