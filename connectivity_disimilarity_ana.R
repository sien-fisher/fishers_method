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
library(uwot)        # UMAP (支持预计算距离矩阵，替代umap包)
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

rownames(conn_mat) <- gsub("/",".",rownames(conn_mat))
colnames(conn_mat) <- gsub("X","",colnames(conn_mat))
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
# 精确命名映射表 (直接查找，不用正则)
# 未匹配到的区域自动归为 "Other"
region_category <- c(
  # Somatosensory
  "1" = "Somatosensory", "2" = "Somatosensory", "3" = "Somatosensory",
  "SII" = "Somatosensory",
  # Prefrontal
  "10" = "Prefrontal", "11" = "Prefrontal", "12" = "Prefrontal",
  "13" = "Prefrontal", "14" = "Prefrontal", "32" = "Prefrontal",
  "44" = "Prefrontal", "45A" = "Prefrontal", "45B" = "Prefrontal",
  "46d" = "Prefrontal", "46v" = "Prefrontal", "8B" = "Prefrontal",
  "8l" = "Prefrontal", "8m" = "Prefrontal", "8r" = "Prefrontal",
  "9" = "Prefrontal", "9.46d" = "Prefrontal", "9.46v" = "Prefrontal",
  # Cingulate
  "23" = "Cingulate", "24a" = "Cingulate", "24b" = "Cingulate","v23" = "Cingulate", 
  "24c" = "Cingulate", "24d" = "Cingulate", "25" = "Cingulate",
  # Paralimbic
  "29.30" = "Paralimbic", "ENTO" = "Paralimbic", "POLE" = "Paralimbic",
  # Association
  "31" = "Association", "5" = "Association", "7A" = "Association",
  "7B" = "Association", "7m" = "Association", "7op" = "Association",
  "AIP" = "Association", "IPa" = "Association",
  # Auditory
  "CORE" = "Auditory", "LB" = "Auditory", "MB" = "Auditory",
  "PBc" = "Auditory", "PBr" = "Auditory", "Pir" = "Auditory",
  # Motor
  "F1" = "Motor", "F2" = "Motor", "F3" = "Motor", "F4" = "Motor",
  "F5" = "Motor", "F6" = "Motor", "F7" = "Motor", "ProM" = "Motor",
  # Gustatory
  "Gu" = "Gustatory",
  # Insula
  "INSULA" = "Insula", "OPAI" = "Insula", "OPRO" = "Insula", "Pi" = "Insula",
  # Dorsal Stream
  "DP" = "Dorsal Stream", "FST" = "Dorsal Stream", "LIP" = "Dorsal Stream",
  "MIP" = "Dorsal Stream", "MST" = "Dorsal Stream", "MTc" = "Dorsal Stream",
  "MTp" = "Dorsal Stream", "PGa" = "Dorsal Stream", "PIP" = "Dorsal Stream",
  "STPc" = "Dorsal Stream", "STPi" = "Dorsal Stream", "STPr" = "Dorsal Stream",
  "TPt" = "Dorsal Stream", "V3LF" = "Dorsal Stream", "V6" = "Dorsal Stream",
  "V6A" = "Dorsal Stream", "VIP" = "Dorsal Stream", "V3A" = "Dorsal Stream",
  # Ventral Stream
  "PERI" = "Ventral Stream", "TEa.ma" = "Ventral Stream", "TEa.mp" = "Ventral Stream",
  "TEad" = "Ventral Stream", "TEav" = "Ventral Stream", "TEO" = "Ventral Stream",
  "TEOm" = "Ventral Stream", "TEpd" = "Ventral Stream", "TEpv" = "Ventral Stream",
  "TH.TF" = "Ventral Stream", "V4Ac" = "Ventral Stream", "V4Apc" = "Ventral Stream",
  "V4t" = "Ventral Stream",
  # Retinotopic Subdivisions
  "Pro.St." = "Retinotopic", "V1c" = "Retinotopic",
  "V1cLF"  = "Retinotopic", "V1cUF"  = "Retinotopic",
  "V1fpLF" = "Retinotopic", "V1fpUF" = "Retinotopic",
  "V1pcLF" = "Retinotopic", "V1pcUF" = "Retinotopic",
  "V2c"    = "Retinotopic",
  "V2cLF"  = "Retinotopic", "V2cUF"  = "Retinotopic",
  "V2fpLF" = "Retinotopic", "V2fpUF" = "Retinotopic",
  "V2pcLF" = "Retinotopic", "V2pcUF" = "Retinotopic",
  "V3c"    = "Retinotopic",
  "V3cLF"  = "Retinotopic", "V3cUF"  = "Retinotopic",
  "V3fpLF" = "Retinotopic", "V3fpUF" = "Retinotopic",
  "V3pcLF" = "Retinotopic", "V3pcUF" = "Retinotopic",
  "V4c"    = "Retinotopic",
  "V4cLF"  = "Retinotopic", "V4cUF"  = "Retinotopic",
  "V4fpLF" = "Retinotopic", "V4fpUF" = "Retinotopic","V4ApcUF" = "Retinotopic",
  "V4pcLF" = "Retinotopic", "V4pcUF" = "Retinotopic",
  "DLP"    = "Retinotopic", "DLPLF"  = "Retinotopic", "DLPUF"  = "Retinotopic"
)
 
# 精确查找分组，未命中的标为 "Other"，并打印出来便于检查
region_group_vec <- region_category[regions]
unmatched <- regions[is.na(region_group_vec)]
if (length(unmatched) > 0) {
  cat("以下区域未在 region_category 中找到匹配，归为 'Other':\n")
  print(unmatched)
}
region_group_vec[is.na(region_group_vec)] <- "Other"
 
region_df <- data.frame(
  Region = regions,
  Group  = unname(region_group_vec),
  stringsAsFactors = FALSE
)
 
cat("\n区域分组分布:\n")
print(sort(table(region_df$Group), decreasing = TRUE))
 
group_colors <- c(
  "Retinotopic"   = "#D9001B",
  "Dorsal Stream" = "#76C376",
  "Ventral Stream"= "#4A90C5",
  "Prefrontal"    = "#9B59B6",
  "Motor"         = "#AD1F80",
  "Cingulate"     = "#F07C74",
  "Association"   = "#2472B5",
  "Somatosensory" = "#FAAD6C",
  "Auditory"      = "#F5D76E",
  "Insula"        = "#1ABC9C",
  "Paralimbic"    = "#E67E22",
  "Gustatory"     = "#95A5A6",
  "Other"         = "#AAAAAA"
)
 

### umap 
library(umap)
run_umap <- function(dist_mat, label, seed = 42) {
  set.seed(seed)
  mat <- as.matrix(dist_mat)
  
  cfg <- umap.defaults
  cfg$input       <- "dist"
  cfg$n_neighbors <- min(15, nrow(mat) - 1)
  cfg$min_dist    <- 0.3
  cfg$n_epochs    <- 500
  cfg$random_state <- seed
  
  umap_res <- umap(mat, config = cfg)
  
  data.frame(
    Region = rownames(mat),
    UMAP1  = umap_res$layout[, 1],
    UMAP2  = umap_res$layout[, 2]
  ) %>%
    left_join(region_df, by = "Region") %>%
    mutate(Method = label)
}
plot_umap <- function(df, title) {
  ggplot(df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = Group, fill = Group), shape = 21, size = 3, alpha = 0.85, stroke = 0.3) +
    geom_text_repel(aes(label = Region), size = 2.5, max.overlaps = 30,
                    segment.alpha = 0.3, segment.size = 0.3) +
    scale_color_manual(values = group_colors, name = "Region Group") +
    scale_fill_manual(values = group_colors, name = "Region Group") +
    theme_void(base_family = "sans") +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      plot.margin     = margin(10, 10, 10, 10)
    ) +
    labs(title = title)
}
umap_raw <- run_umap(dist_mat_raw, "FLN Bray-Curtis")
umap_clr <- run_umap(dist_mat_clr, "CLR Aitchison")
p_umap_raw <- plot_umap(umap_raw, "UMAP: FLN Bray-Curtis dissimilarity")
p_umap_clr <- plot_umap(umap_clr, "UMAP: CLR Aitchison distance")
p_umap_compare <- p_umap_raw + p_umap_clr +
  plot_layout(guides = "collect") & theme(legend.position = "right")





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
