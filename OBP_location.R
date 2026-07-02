# ==========================================
# 1. 加载包
# ==========================================
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

library(ggplot2)
library(dplyr)
library(ggrepel)

# ==========================================
# 2. 读取数据 (交互式选择文件)
# ==========================================
print("请选择你的 OBP 基因数据文件 (CSV/Excel)...")
file_path <- file.choose()
# 兼容读取 CSV
obp_data <- read.csv(file_path, stringsAsFactors = FALSE)

# 提取基因名(第1列)和染色体(第4列)
obp_data <- obp_data[, c(1, 4)]
colnames(obp_data) <- c("GeneName", "Chromosome")

print("请选择你的 染色体长度 TXT 文件 (BMSK_chr_lengths.txt)...")
txt_path <- file.choose()
chr_len <- read.table(txt_path, header=FALSE, col.names=c("Chr", "Length"))

# ==========================================
# 3. 数据清洗与虚拟坐标生成
# ==========================================
obp_data$Chr_num <- as.numeric(gsub("[^0-9]", "", obp_data$Chromosome))
chr_len$Chr_num <- as.numeric(gsub("[^0-9]", "", chr_len$Chr))
obp_data <- obp_data %>% filter(!is.na(Chr_num))

relevant_chrs <- sort(unique(obp_data$Chr_num))
obp_plot_data <- obp_data %>% filter(Chr_num %in% relevant_chrs)
chr_plot_len <- chr_len %>% filter(Chr_num %in% relevant_chrs)

obp_plot_data$Chr_fact <- factor(obp_plot_data$Chr_num, levels = relevant_chrs)
chr_plot_len$Chr_fact <- factor(chr_plot_len$Chr_num, levels = relevant_chrs)

# 生成展示用的均匀坐标
obp_plot_data <- obp_plot_data %>%
  left_join(chr_plot_len[, c("Chr_num", "Length")], by = "Chr_num") %>%
  group_by(Chr_num) %>%
  arrange(GeneName) %>% 
  mutate(
    # 缩小分布范围，让基因簇看起来更紧凑一点
    Plot_Start = seq(from = max(Length) * 0.15, to = max(Length) * 0.85, length.out = n())
  ) %>%
  ungroup()

# ==========================================
# 4. 【核心优化】智能成簇 (Clustering) 计算
# ==========================================
# 定义一个阈值：大于等于 4 个基因被视为一个 Cluster
cluster_threshold <- 4

# 计算哪些染色体形成了 Cluster
clusters_info <- obp_plot_data %>%
  group_by(Chr_fact) %>%
  summarise(
    n_genes = n(),
    min_y = min(Plot_Start) / 1e6,
    max_y = max(Plot_Start) / 1e6,
    mean_y = mean(Plot_Start) / 1e6
  ) %>%
  filter(n_genes >= cluster_threshold) %>%
  mutate(Cluster_Name = paste0("Cluster ", row_number(), "\n(", n_genes, " genes)"))

# 筛选出不成簇的独立基因（用于单独显示名字）
singles_info <- obp_plot_data %>%
  group_by(Chr_fact) %>%
  filter(n() < cluster_threshold) %>%
  ungroup()

# ==========================================
# 5. SCI 级美化绘图
# ==========================================
p <- ggplot() +
  # [底层] 画染色体骨架 (优雅的浅灰)
  geom_segment(data = chr_plot_len,
               aes(x = Chr_fact, xend = Chr_fact, y = 0, yend = Length / 1e6),
               color = "#E5E7EB", linewidth = 5, lineend = "round") +
  
  # [中层] 画所有基因的具体位置 (细红线横切染色体)
  geom_segment(data = obp_plot_data,
               aes(x = as.numeric(Chr_fact) - 0.2, 
                   xend = as.numeric(Chr_fact) + 0.2, 
                   y = Plot_Start / 1e6, yend = Plot_Start / 1e6),
               color = "#E41A1C", linewidth = 0.8) +
  
  # ================== 绘制 Cluster 学术括号 ==================
# 1. 括号的主垂直线
geom_segment(data = clusters_info,
             aes(x = as.numeric(Chr_fact) + 0.35, xend = as.numeric(Chr_fact) + 0.35, 
                 y = min_y, yend = max_y), color = "black", linewidth = 0.6) +
  # 2. 括号顶部的小横线
  geom_segment(data = clusters_info,
               aes(x = as.numeric(Chr_fact) + 0.25, xend = as.numeric(Chr_fact) + 0.35, 
                   y = min_y, yend = min_y), color = "black", linewidth = 0.6) +
  # 3. 括号底部的小横线
  geom_segment(data = clusters_info,
               aes(x = as.numeric(Chr_fact) + 0.25, xend = as.numeric(Chr_fact) + 0.35, 
                   y = max_y, yend = max_y), color = "black", linewidth = 0.6) +
  # 4. 添加 "Cluster X" 文字
  geom_text(data = clusters_info,
            aes(x = as.numeric(Chr_fact) + 0.45, y = mean_y, label = Cluster_Name),
            hjust = 0, size = 3.5, fontface = "bold", color = "#333333", lineheight = 0.9) +
  
  # ================== 绘制独立基因名称 ==================
geom_text_repel(data = singles_info,
                aes(x = as.numeric(Chr_fact) + 0.2, y = Plot_Start / 1e6, label = GeneName),
                size = 3.5, fontface = "italic", color = "black", direction = "y", 
                nudge_x = 0.4, hjust = 0, segment.size = 0.4, segment.color = "grey50") +
  
  # [外围] 坐标系与主题设定
  scale_y_reverse(name = "Physical Position (Mb)", 
                  limits = c(max(chr_plot_len$Length)/1e6, -2), expand = c(0.02, 0)) + 
  scale_x_discrete(name = "Chromosome") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 80, b = 10, l = 10) # 右侧留足80的空白，防止Cluster文字被裁
  )

# 打印并保存高清矢量图 (推荐发文章用 PDF)
print(p)
ggsave("OBP_Genomic_Locations_SCI_Style.pdf", plot = p, width = 10, height = 7, dpi = 300)
ggsave("OBP_Genomic_Locations_SCI_Style.png", plot = p, width = 10, height = 7, dpi = 300)