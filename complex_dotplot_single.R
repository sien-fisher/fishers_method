# complex_dotplot_single: 基础逻辑实现
# 展示单个基因在各 cluster 中按 group 分组的 dotplot
# x轴: group (如 orig.ident), y轴: celltype/cluster, facet: 无（或自选）

complex_dotplot_single <- function(
  seu_obj,
  feature,
  groups,           # meta.data列名，用于x轴分组（如 "orig.ident"）
  celltypes = NULL, # 筛选的cluster，NULL则用全部levels
  split.by = NULL,  # 可选，额外facet变量
  vmin = NULL,
  vmax = NULL,
  color.palette = NULL,
  dot.scale = 10    # 点最大尺寸
) {
  library(Seurat)
  library(ggplot2)
  library(dplyr)

  # ── 1. 确定要展示的 celltypes ──────────────────────────────
  all_ct <- levels(seu_obj)  # Seurat默认identity levels
  if (!is.null(celltypes)) {
    all_ct <- intersect(celltypes, all_ct)
  }

  # ── 2. 提取元数据 + 基因表达量 ────────────────────────────
  meta <- seu_obj@meta.data
  expr <- FetchData(seu_obj, vars = c(feature, groups, if (!is.null(split.by)) split.by))
  expr$celltype <- Idents(seu_obj)
  expr <- expr[expr$celltype %in% all_ct, ]

  # ── 3. 分组计算 avg.exp 和 pct.exp ────────────────────────
  group_vars <- c("celltype", groups, if (!is.null(split.by)) split.by)

  dot_data <- expr %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      avg.exp  = mean(expm1(.data[[feature]])),   # 还原log1p后取均值
      pct.exp  = mean(.data[[feature]] > 0) * 100,
      .groups  = "drop"
    )

  # scale avg.exp to avg.exp.scaled（按基因在所有格中归一化）
  exp_range <- range(dot_data$avg.exp, na.rm = TRUE)
  if (!is.null(vmin)) exp_range[1] <- vmin
  if (!is.null(vmax)) exp_range[2] <- vmax
  dot_data$avg.exp.scaled <- (dot_data$avg.exp - exp_range[1]) /
                              (exp_range[2] - exp_range[1])
  dot_data$avg.exp.scaled <- pmin(pmax(dot_data$avg.exp.scaled, 0), 1)

  # ── 4. 因子化保持顺序 ─────────────────────────────────────
  dot_data$celltype <- factor(dot_data$celltype, levels = all_ct)
  dot_data$groups   <- dot_data[[groups]]

  # ── 5. 配色 ───────────────────────────────────────────────
color.palette <- colorRampPalette(
  c("#89C4E1", "#B8D9E8", "#E8C8DC", "#F4A7C3")
)(255)

  # ── 6. 画图 ───────────────────────────────────────────────
  p <- ggplot(dot_data, aes(x = groups, y = celltype)) +
    geom_point(aes(color = avg.exp.scaled, size = pct.exp), alpha = 0.9) +
    scale_color_gradientn(
      colours  = color.palette,
      name     = "Avg Exp\n(scaled)",
      limits   = c(0, 1)
    ) +
    scale_size_continuous(
      range = c(0, dot.scale),
      name  = "Pct Exp (%)"
    ) +
    labs(
      title = feature,
      x     = groups,
      y     = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major  = element_line(color = "grey92"),
      panel.grid.minor  = element_blank(),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y       = element_text(size = 12),
      plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position   = "right",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      strip.text        = element_text(size = 12, face = "bold"),
      strip.background  = element_rect(fill = "grey90", color = "black")
    ) +
    guides(
      size  = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black")),
      color = guide_colorbar(barwidth = 1, barheight = 6)
    )

  # ── 7. 可选 split.by facet ────────────────────────────────
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~ .data[[split.by]], ncol = length(unique(dot_data[[split.by]])))
  }

  return(p)
}


# ─── 多基因版本 ─────────────────────────────────────────────────────────────
complex_dotplot_multiple <- function(
  seu_obj,
  features,
  groups,
  celltypes    = NULL,
  split.by     = NULL,
  color.palette = NULL,
  vmin         = NULL,
  vmax         = NULL,
  dot.scale    = 10,
  strip.color  = NULL   # 每个celltype facet的顶部背景色，长度=celltype数
) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  # ── 1. 逐基因收集数据 ─────────────────────────────────────
  all_ct <- levels(seu_obj)
  if (!is.null(celltypes)) all_ct <- intersect(celltypes, all_ct)

  plot_list <- lapply(features, function(feat) {
    p_tmp <- complex_dotplot_single(
      seu_obj   = seu_obj,
      feature   = feat,
      groups    = groups,
      celltypes = celltypes,
      split.by  = split.by,
      vmin      = vmin,
      vmax      = vmax,
      color.palette = color.palette,
      dot.scale = dot.scale
    )
    df        <- p_tmp$data
    df$gene   <- feat
    df
  })

  all_data <- do.call(rbind, plot_list)
  all_data$gene     <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = all_ct)

  # ── 2. 配色 ───────────────────────────────────────────────
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(
      c("grey80", "lemonchiffon1", "indianred1", "darkred")
    )(255)
  }

  # ── 3. 主图：y=gene, x=group, facet=celltype ──────────────
  p <- ggplot(all_data, aes(x = groups, y = gene)) +
    geom_tile(fill = "white", color = "white") +          # 保持背景均匀
    geom_point(aes(color = avg.exp.scaled, size = pct.exp), alpha = 0.9) +
    scale_color_gradientn(
      colours = color.palette,
      name    = "Avg Exp\n(scaled)",
      limits  = c(0, 1)
    ) +
    scale_size_continuous(
      range = c(0, dot.scale),
      name  = "Pct Exp (%)"
    ) +
    facet_wrap(~ celltype, ncol = length(all_ct)) +
    labs(x = "", y = "") +
    theme_bw(base_size = 12) +
    theme(
      panel.background  = element_rect(fill = "white", colour = "black"),
      panel.grid.major  = element_line(color = "grey92"),
      panel.grid.minor  = element_blank(),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y       = element_text(size = 12),
      legend.position   = "right",
      legend.title      = element_text(size = 10),
      strip.text        = element_text(size = 13, face = "bold", colour = "black"),
      strip.background  = element_rect(fill = "grey88", colour = "black")
    ) +
    guides(
      size  = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black")),
      color = guide_colorbar(barwidth = 1, barheight = 6)
    )

  # ── 4. 自定义 strip 背景色 ────────────────────────────────
  if (!is.null(strip.color)) {
    library(grid)
    library(gtable)
    gb   <- ggplot_build(p)
    gt   <- ggplot_gtable(gb)
    stri <- grep("strip-t", gt$layout$name)
    # 循环赋色（strip.color长度需 >= celltype数）
    cols <- rep_len(strip.color, length(stri))
    for (i in seq_along(stri)) {
      j <- which(grepl("rect", gt$grobs[[stri[i]]]$grobs[[1]]$childrenOrder))
      gt$grobs[[stri[i]]]$grobs[[1]]$children[[j]]$gp$fill <- cols[i]
    }
    grid::grid.draw(gt)
    invisible(gt)
  } else {
    print(p)
    invisible(p)
  }
}


# ─── 使用示例 ───────────────────────────────────────────────────────────────
# 单基因（x=orig.ident, 各cluster行，可选split.by="Phase"）
Idents(seu) <- seu$refined_annotation
celltype <- levels(factor(seu$refined_annotation))

nothit_genes <- c("BASP1-AS1", "C1orf94", "PTGER3", "HOXA7", "LINC00229")

for (feat in nothit_genes) {
  message(sprintf("PLOT  '%s'", feat))

  p <- tryCatch(
    complex_dotplot_single(
      seu_obj   = seu,
      feature   = feat,
      groups    = "age",
      celltypes = celltype
    ),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )

  if (is.null(p)) next

  fname <- paste0(gsub("[^A-Za-z0-9_\\-]", "_", feat), "_vanco.pdf")
  ggsave(fname, plot = p, width = 10, height = 6, device = cairo_pdf)
}
                  
    
                  
# 多基因（y=gene, x=orig.ident, facet=celltype）
# png("dotplot_multiple.png", width=14, height=6, units="in", res=150)
# complex_dotplot_multiple(
#   seu_obj   = iri.integrated,
#   features  = c("Havcr1","Slc34a1","Lrp2"),
#   groups    = "orig.ident",
#   celltypes = c("PTS1","PTS2","PTS3","NewPT1","NewPT2")
# )
# dev.off()
