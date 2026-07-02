# ==========================================================================
# OBP Gene Family: Genomic Localisation & Visualisation
# Bombyx mori (silkworm) odorant-binding protein genes
# Course: Biological Omics & Big Data
# ==========================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

# ── Configuration ──────────────────────────────────────────────────────────
# Adjust these paths if your data files are located elsewhere.
DATA_DIR   <- "data"
FIG_DIR    <- "figures"
dir.create(FIG_DIR, showWarnings = FALSE)

GENE_FILE  <- file.path(DATA_DIR, "obp_locations.txt")
CHR_FILE   <- file.path(DATA_DIR, "BMSK_chr_lengths.txt")

# ── 1. Load Data ───────────────────────────────────────────────────────────

# Gene positions: expected columns: GeneName (col 1), Chromosome (col 4)
obp_data <- read.csv(GENE_FILE, stringsAsFactors = FALSE)
obp_data <- obp_data[, c(1, 4)]
colnames(obp_data) <- c("GeneName", "Chromosome")

# Chromosome lengths: two columns (Chr, Length)
chr_len <- read.table(CHR_FILE, header = FALSE,
                      col.names = c("Chr", "Length"))

# ── 2. Parse & Clean ───────────────────────────────────────────────────────

# Extract numeric chromosome IDs
obp_data$Chr_num <- as.numeric(gsub("[^0-9]", "", obp_data$Chromosome))
chr_len$Chr_num  <- as.numeric(gsub("[^0-9]", "", chr_len$Chr))

# Keep only chromosomes that contain at least one OBP gene
obp_data <- obp_data %>% filter(!is.na(Chr_num))
relevant_chrs <- sort(unique(obp_data$Chr_num))

obp_plot_data <- obp_data %>% filter(Chr_num %in% relevant_chrs)
chr_plot_len  <- chr_len  %>% filter(Chr_num %in% relevant_chrs)

obp_plot_data$Chr_fact <- factor(obp_plot_data$Chr_num, levels = relevant_chrs)
chr_plot_len$Chr_fact  <- factor(chr_plot_len$Chr_num,  levels = relevant_chrs)

# ── 3. Compute Plot Positions ──────────────────────────────────────────────
# Genes are evenly spaced along each chromosome for clear visualisation.

obp_plot_data <- obp_plot_data %>%
  left_join(chr_plot_len[, c("Chr_num", "Length")], by = "Chr_num") %>%
  group_by(Chr_num) %>%
  arrange(GeneName) %>%
  mutate(
    Plot_Start = seq(
      from = max(Length) * 0.15,
      to   = max(Length) * 0.85,
      length.out = n()
    )
  ) %>%
  ungroup()

# ── 4. Identify Gene Clusters ──────────────────────────────────────────────
# Chromosomes carrying >= CLUSTER_THRESHOLD genes are annotated as clusters.

CLUSTER_THRESHOLD <- 4

clusters_info <- obp_plot_data %>%
  group_by(Chr_fact) %>%
  summarise(
    n_genes = n(),
    min_y   = min(Plot_Start) / 1e6,
    max_y   = max(Plot_Start) / 1e6,
    mean_y  = mean(Plot_Start) / 1e6
  ) %>%
  filter(n_genes >= CLUSTER_THRESHOLD) %>%
  mutate(Cluster_Name = paste0("Cluster ", row_number(), "\n(", n_genes, " genes)"))

# Isolated genes (chromosomes with fewer than threshold genes) get labels
singles_info <- obp_plot_data %>%
  group_by(Chr_fact) %>%
  filter(n() < CLUSTER_THRESHOLD) %>%
  ungroup()

# ── 5. SCI-Style Chromosome Ideogram ───────────────────────────────────────

p <- ggplot() +

  # Background: chromosome backbones (light grey bars)
  geom_segment(data = chr_plot_len,
               aes(x = Chr_fact, xend = Chr_fact, y = 0, yend = Length / 1e6),
               colour = "#E5E7EB", linewidth = 5, lineend = "round") +

  # Gene positions: red tick marks
  geom_segment(data = obp_plot_data,
               aes(x = as.numeric(Chr_fact) - 0.2,
                   xend = as.numeric(Chr_fact) + 0.2,
                   y = Plot_Start / 1e6, yend = Plot_Start / 1e6),
               colour = "#E41A1C", linewidth = 0.8) +

  # ── Cluster annotations ──
  # Vertical bracket line
  geom_segment(data = clusters_info,
               aes(x = as.numeric(Chr_fact) + 0.35,
                   xend = as.numeric(Chr_fact) + 0.35,
                   y = min_y, yend = max_y),
               colour = "black", linewidth = 0.6) +
  # Upper horizontal tick
  geom_segment(data = clusters_info,
               aes(x = as.numeric(Chr_fact) + 0.25,
                   xend = as.numeric(Chr_fact) + 0.35,
                   y = min_y, yend = min_y),
               colour = "black", linewidth = 0.6) +
  # Lower horizontal tick
  geom_segment(data = clusters_info,
               aes(x = as.numeric(Chr_fact) + 0.25,
                   xend = as.numeric(Chr_fact) + 0.35,
                   y = max_y, yend = max_y),
               colour = "black", linewidth = 0.6) +
  # Cluster label
  geom_text(data = clusters_info,
            aes(x = as.numeric(Chr_fact) + 0.45, y = mean_y,
                label = Cluster_Name),
            hjust = 0, size = 3.5, fontface = "bold",
            colour = "#333333", lineheight = 0.9) +

  # ── Isolated gene labels (ggrepel, italicised) ──
  geom_text_repel(data = singles_info,
                  aes(x = as.numeric(Chr_fact) + 0.2,
                      y = Plot_Start / 1e6, label = GeneName),
                  size = 3.5, fontface = "italic", colour = "black",
                  direction = "y", nudge_x = 0.4, hjust = 0,
                  segment.size = 0.4, segment.colour = "grey50") +

  # ── Coordinate system ──
  scale_y_reverse(
    name   = "Physical Position (Mb)",
    limits = c(max(chr_plot_len$Length) / 1e6, -2),
    expand = c(0.02, 0)
  ) +
  scale_x_discrete(name = "Chromosome") +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    axis.line.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.x      = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y      = element_text(size = 11, colour = "black"),
    axis.title       = element_text(size = 14, face = "bold"),
    plot.margin      = margin(t = 20, r = 80, b = 10, l = 10)
  )

# ── 6. Export ──────────────────────────────────────────────────────────────

print(p)

pdf_out <- file.path(FIG_DIR, "OBP_Genomic_Locations_SCI_Style.pdf")
png_out <- file.path(FIG_DIR, "OBP_Genomic_Locations_SCI_Style.png")
ggsave(pdf_out, plot = p, width = 10, height = 7, dpi = 300)
ggsave(png_out, plot = p, width = 10, height = 7, dpi = 300)

cat(sprintf("Figures saved:\n  %s\n  %s\n", pdf_out, png_out))

# ── Session Info ───────────────────────────────────────────────────────────
cat("\n── Session Info ──\n")
sessionInfo()
