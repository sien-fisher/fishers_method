# ============================================================
# 连接组 Dissimilarity Matrix 分析
#
# 数据: conn_mat = 113x113 有向FLN矩阵 (行=源区域, 列=目标区域)
#       每个区域的特征向量 = afferent(列) + efferent(行) 拼接 = 226维
#
# 分析流程:
#   1. 构建226维特征向量
#   2. Bray-Curtis dissimilarity matrix (113x113)
#   3. UMAP / t-SNE 降维可视化
#   4. Hierarchical clustering
#   5. (接口预留) cell type composition dissimilarity matrix
#      + 跨模态 Mantel test 相关
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(tibble)
library(vegan)       # Bray-Curtis + Mantel test
library(umap)        # UMAP
library(Rtsne)       # t-SNE
library(dendextend)  # 美化dendrogram
library(compositions)

# 命名空间锁定，防止MASS/stats冲突
select  <- dplyr::select
filter  <- dplyr::filter
mutate  <- dplyr::mutate

# ============================================================
# 输入数据
# conn_mat: 113x113 数值矩阵
#           行名 = 源区域(source), 列名 = 目标区域(target)
#           值 = FLN (已经去掉自连接对角线，或保留均可，下面会处理)
# ============================================================
# conn_mat <- readRDS("your_conn_mat.rds")   # ← 替换为你的实际读入方式
# conn_mat <- read.csv("your_conn_mat.csv", row.names = 1) %>% as.matrix()

# ============================================================
# 步骤 1: 构建226维特征向量 (afferent || efferent)
# ============================================================
# 对角线设为0 (不考虑自连接)
diag(conn_mat) <- 0

# afferent: 某区域接收来自其他区域的输入 -> 取该区域所在的列向量
# efferent: 某区域向其他区域的输出       -> 取该区域所在的行向量
regions <- rownames(conn_mat)
stopifnot(identical(regions, colnames(conn_mat)))

feat_mat <- cbind(
  t(conn_mat),      # afferent: 每列转置后 -> 行=region, 列=source
  conn_mat          # efferent: 每行        -> 行=region, 列=target
)
# 列名加前缀避免afferent/efferent列名重复
colnames(feat_mat) <- c(
  paste0("aff_", regions),
  paste0("eff_", regions)
)
cat("特征矩阵维度:", nrow(feat_mat), "regions x", ncol(feat_mat), "features\n")

# ============================================================
# 步骤 2: CLR变换 (组成数据处理, 和FLN逻辑一致)
# ============================================================
feat_raw <- feat_mat   # 保留原始FLN值版本

# CLR版本
feat_clr <- feat_mat
if (any(feat_clr == 0)) {
  feat_clr[feat_clr == 0] <- 1e-6
  feat_clr <- t(prop.table(t(feat_clr), margin = 1))  # 按行归一化
}
feat_clr <- as.matrix(compositions::clr(feat_clr))

cat("原始FLN版本 & CLR版本特征矩阵均已准备好\n")

# ============================================================
# 步骤 3: Bray-Curtis Dissimilarity Matrix
# ============================================================
# vegdist 期望: 行=样本(region), 列=特征
dist_bc_raw <- vegdist(feat_raw, method = "bray")
dist_bc_clr <- vegdist(feat_clr, method = "euclidean")  # CLR后用欧氏距离等价于Aitchison距离

dist_mat_raw <- as.matrix(dist_bc_raw)
dist_mat_clr <- as.matrix(dist_bc_clr)

cat("Bray-Curtis dissimilarity matrix (FLN raw):", dim(dist_mat_raw), "\n")
cat("Aitchison distance matrix (CLR):           ", dim(dist_mat_clr), "\n")

# ============================================================
# 步骤 4: 区域分组标签 (用于后续着色)
# ============================================================
# 根据你的区域命名体系定义分组，修改这里的逻辑即可
get_region_group <- function(regions) {
  case_when(
    grepl("V1", regions)          ~ "V1",
    grepl("V2", regions)          ~ "V2",
    grepl("V3", regions)          ~ "V3",
    grepl("V4|DLP", regions)      ~ "V4",
    grepl("MT|FST|MST", regions)  ~ "Dorsal Stream",
    grepl("TEO|TE|IT|Inf", regions) ~ "Ventral Stream",
    grepl("LIP|VIP|MIP|7|PG", regions) ~ "Parietal",
    grepl("F[1-7]|PMC|M1", regions)    ~ "Frontal Motor",
    grepl("PFC|46|9|10|11|12|13|14|47", regions) ~ "PFC",
    TRUE                          ~ "Other"
  )
}

region_df <- data.frame(
  Region = regions,
  Group  = get_region_group(regions),
  stringsAsFactors = FALSE
)

group_colors <- c(
  "V1"             = "#D9001B",
  "V2"             = "#F07C74",
  "V3"             = "#FAAD6C",
  "V4"             = "#F5D76E",
  "Dorsal Stream"  = "#76C376",
  "Ventral Stream" = "#4A90C5",
  "Parietal"       = "#2472B5",
  "Frontal Motor"  = "#AD1F80",
  "PFC"            = "#9B59B6",
  "Other"          = "#AAAAAA"
)

# ============================================================
# 步骤 5: UMAP (基于dissimilarity matrix)
# ============================================================
run_umap <- function(dist_mat, label, seed = 42) {
  set.seed(seed)
  umap_config <- umap.defaults
  umap_config$input <- "dist"
  umap_config$n_neighbors <- min(15, nrow(dist_mat) - 1)
  umap_config$min_dist    <- 0.3
  umap_config$metric      <- "precomputed"

  umap_res <- umap(dist_mat, config = umap_config)
  df <- data.frame(
    Region = rownames(dist_mat),
    UMAP1  = umap_res$layout[, 1],
    UMAP2  = umap_res$layout[, 2]
  ) %>%
    left_join(region_df, by = "Region") %>%
    mutate(Method = label)
  df
}

umap_raw <- run_umap(dist_mat_raw, "FLN Bray-Curtis")
umap_clr <- run_umap(dist_mat_clr, "CLR Aitchison")

plot_umap <- function(df, title) {
  centroids <- df %>%
    group_by(Group) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

  ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Group, fill = Group), shape = 21, size = 3, alpha = 0.85, stroke = 0.3) +
    geom_text_repel(aes(label = Region), size = 2.5, max.overlaps = 30,
                     segment.alpha = 0.3, segment.size = 0.3) +
    scale_color_manual(values = group_colors, name = "Region Group") +
    scale_fill_manual(values = group_colors, name = "Region Group") +
    theme_void(base_family = "sans") +
    theme(
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      plot.margin   = margin(10, 10, 10, 10)
    ) +
    labs(title = title)
}

p_umap_raw <- plot_umap(umap_raw, "UMAP: FLN Bray-Curtis dissimilarity")
p_umap_clr <- plot_umap(umap_clr, "UMAP: CLR Aitchison distance")
p_umap_compare <- p_umap_raw + p_umap_clr + plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_umap_compare)

# ============================================================
# 步骤 6: t-SNE (基于dissimilarity matrix)
# ============================================================
run_tsne <- function(dist_mat, label, seed = 42) {
  set.seed(seed)
  tsne_res <- Rtsne(
    dist_mat,
    is_distance  = TRUE,
    perplexity   = min(30, floor((nrow(dist_mat) - 1) / 3)),
    max_iter     = 1000,
    check_duplicates = FALSE
  )
  data.frame(
    Region = rownames(dist_mat),
    tSNE1  = tsne_res$Y[, 1],
    tSNE2  = tsne_res$Y[, 2]
  ) %>%
    left_join(region_df, by = "Region") %>%
    mutate(Method = label)
}

tsne_raw <- run_tsne(dist_mat_raw, "FLN Bray-Curtis")
tsne_clr <- run_tsne(dist_mat_clr, "CLR Aitchison")

plot_tsne <- function(df, title) {
  ggplot(df, aes(tSNE1, tSNE2)) +
    geom_point(aes(color = Group, fill = Group), shape = 21, size = 3, alpha = 0.85, stroke = 0.3) +
    geom_text_repel(aes(label = Region), size = 2.5, max.overlaps = 30,
                     segment.alpha = 0.3, segment.size = 0.3) +
    scale_color_manual(values = group_colors, name = "Region Group") +
    scale_fill_manual(values = group_colors, name = "Region Group") +
    theme_void(base_family = "sans") +
    theme(
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      plot.margin   = margin(10, 10, 10, 10)
    ) +
    labs(title = title)
}

p_tsne_raw <- plot_tsne(tsne_raw, "t-SNE: FLN Bray-Curtis dissimilarity")
p_tsne_clr <- plot_tsne(tsne_clr, "t-SNE: CLR Aitchison distance")
p_tsne_compare <- p_tsne_raw + p_tsne_clr + plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_tsne_compare)

# ============================================================
# 步骤 7: Hierarchical Clustering
# ============================================================
run_hclust <- function(dist_obj, label, method = "ward.D2") {
  hc <- hclust(dist_obj, method = method)
  hc$labels <- regions
  list(hc = hc, label = label)
}

hc_raw <- run_hclust(dist_bc_raw, "FLN Bray-Curtis")
hc_clr <- run_hclust(dist_bc_clr, "CLR Aitchison")

plot_hclust <- function(hc_obj, region_df, title, n_cut = 8) {
  dend <- as.dendrogram(hc_obj$hc)

  # 按分组给叶片着色
  group_map <- setNames(region_df$Group, region_df$Region)
  label_colors <- group_colors[group_map[labels(dend)]]
  label_colors[is.na(label_colors)] <- "#AAAAAA"

  dend <- dend %>%
    set("labels_cex", 0.55) %>%
    set("labels_colors", label_colors)

  # 切割成 n_cut 个cluster，给树枝着色
  dend <- color_branches(dend, k = n_cut, col = scales::hue_pal()(n_cut))

  cairo_pdf(file.path(out_dir, paste0("Hclust_", gsub(" ", "_", hc_obj$label), ".pdf")),
            width = 20, height = 8)
  par(mar = c(8, 4, 3, 2))
  plot(dend, main = title, ylab = "Distance")
  legend("topright", legend = names(group_colors),
         fill = group_colors, cex = 0.6, ncol = 2, title = "Region Group")
  dev.off()

  png(file.path(out_dir, paste0("Hclust_", gsub(" ", "_", hc_obj$label), ".png")),
      width = 2400, height = 960, res = 150)
  par(mar = c(8, 4, 3, 2))
  plot(dend, main = title, ylab = "Distance")
  legend("topright", legend = names(group_colors),
         fill = group_colors, cex = 0.6, ncol = 2, title = "Region Group")
  dev.off()

  # 返回cluster membership (供后续分析用)
  data.frame(
    Region  = regions,
    Cluster = cutree(hc_obj$hc, k = n_cut)[regions]
  ) %>% left_join(region_df, by = "Region")
}

out_dir <- "Dissimilarity_Analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cluster_raw <- plot_hclust(hc_raw, region_df, "Hierarchical Clustering: FLN Bray-Curtis")
cluster_clr <- plot_hclust(hc_clr, region_df, "Hierarchical Clustering: CLR Aitchison")

cat("\n[FLN] Cluster membership:\n"); print(table(cluster_raw$Cluster, cluster_raw$Group))
cat("\n[CLR] Cluster membership:\n"); print(table(cluster_clr$Cluster, cluster_clr$Group))

# ============================================================
# 步骤 8: 保存
# ============================================================
write.csv(dist_mat_raw, file.path(out_dir, "Dissimilarity_matrix_FLN_BrayCurtis.csv"))
write.csv(dist_mat_clr, file.path(out_dir, "Dissimilarity_matrix_CLR_Aitchison.csv"))
write.csv(bind_rows(umap_raw, umap_clr), file.path(out_dir, "UMAP_coords.csv"), row.names = FALSE)
write.csv(bind_rows(tsne_raw, tsne_clr), file.path(out_dir, "tSNE_coords.csv"), row.names = FALSE)
write.csv(cluster_raw, file.path(out_dir, "Hclust_membership_FLN.csv"), row.names = FALSE)
write.csv(cluster_clr, file.path(out_dir, "Hclust_membership_CLR.csv"), row.names = FALSE)

ggsave(file.path(out_dir, "UMAP_compare.pdf"), p_umap_compare, width = 18, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, "UMAP_compare.png"), p_umap_compare, width = 18, height = 8, dpi = 300)
ggsave(file.path(out_dir, "tSNE_compare.pdf"), p_tsne_compare, width = 18, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, "tSNE_compare.png"), p_tsne_compare, width = 18, height = 8, dpi = 300)

# ============================================================
# 步骤 9: 跨模态相关 (Mantel test)
# ============================================================
# 等你提供cell type composition矩阵后，取消下方注释运行
# comp_mat: region x cell_type 矩阵 (行名需与 conn_mat 行名一致)

# --- 接口预留 ---
# shared_regions <- intersect(rownames(comp_mat), regions)
# cat("\n共享区域数:", length(shared_regions), "\n")
#
# # 按共享区域对齐
# comp_sub <- comp_mat[shared_regions, ]
# dist_conn_sub <- as.dist(dist_mat_raw[shared_regions, shared_regions])
#
# # cell type composition 的 Bray-Curtis dissimilarity
# dist_comp_bc <- vegdist(comp_sub, method = "bray")
#
# # Mantel test: 两个距离矩阵之间的相关性
# # H0: connectivity距离 与 composition距离 无关
# set.seed(42)
# mantel_res <- mantel(dist_conn_sub, dist_comp_bc, method = "spearman",
#                       permutations = 999, na.rm = TRUE)
# cat("\n========== Mantel test: Connectivity vs Composition ==========\n")
# print(mantel_res)
#
# # 散点图: 两个距离矩阵的元素级相关
# dist_conn_vec <- as.vector(dist_conn_sub)
# dist_comp_vec <- as.vector(dist_comp_bc)
# valid <- !is.na(dist_conn_vec) & !is.na(dist_comp_vec)
# cor_val <- cor(dist_conn_vec[valid], dist_comp_vec[valid], method = "spearman")
#
# p_mantel <- ggplot(data.frame(Connectivity = dist_conn_vec[valid],
#                                Composition  = dist_comp_vec[valid]),
#                     aes(Connectivity, Composition)) +
#   geom_point(alpha = 0.2, size = 0.8, color = "#2472B5") +
#   geom_smooth(method = "lm", color = "#D9001B", linewidth = 0.8, se = TRUE) +
#   annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
#             label = sprintf("Mantel r = %.3f, p = %.4f (999 perms)",
#                              mantel_res$statistic, mantel_res$signif),
#             size = 4, fontface = "italic") +
#   theme_bw(base_size = 12) +
#   labs(title = "Connectivity vs Composition dissimilarity",
#        x = "FLN Bray-Curtis distance", y = "Composition Bray-Curtis distance")
#
# print(p_mantel)
# ggsave(file.path(out_dir, "Mantel_Connectivity_vs_Composition.pdf"), p_mantel,
#        width = 7, height = 6, device = cairo_pdf)
# ggsave(file.path(out_dir, "Mantel_Connectivity_vs_Composition.png"), p_mantel,
#        width = 7, height = 6, dpi = 300)

cat("\n全部分析完成，输出目录:", out_dir, "\n")
cat("Cell type composition接口已预留在步骤9，提供数据后取消注释即可运行\n")
