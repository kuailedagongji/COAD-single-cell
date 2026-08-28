
library(CytoTRACE2)
library(tidyverse)
library(Seurat)
options(timeout = 10000)
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
#remotes::install_local("digitalcytometry-cytotrace2-6e8041f.tar.gz",subdir = "cytotrace2_r", upgrade = F,dependencies = T)
sce2 <- readRDS("上皮评分后.rds")
Idents(sce2)=sce2$geneSet.Type
#sce2 <-clean
future::plan("multisession", workers = 4)
sce2@meta.data$CB <- rownames(sce2@meta.data)
#sce2$celltype<-Idents(sce2)
sample_CB <- sce2@meta.data %>% 
  group_by(geneSet.Type) %>% 
  sample_frac(0.8)
#options(future.globals.maxSize = 5 * 1024^3)
sce3 <- subset(sce2,CB %in% sample_CB$CB) 
rm(sce2)
sce3
# An object of class Seurat 
#######输入seurat 对象###########
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)
cytotrace2_result_sce
save(cytotrace2_result_sce,file = "cytotrace2.rda")
load("cytotrace2.rda")
# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = sce3@meta.data$geneSet.Type) %>% 
  set_rownames(., colnames(sce3))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)


# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno
p5 <- plots$Phenotype_UMAP
library(patchwork)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
p1
p2
p3
p4
p5
pdf(file = "CytoTrace_umap1.pdf", width = 5, height = 4)
FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 1.5) + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16), 
        plot.title = element_text(size = 16, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)

dev.off()

# 假设 p5 是一个标准的 DimPlot 或 ggplot 对象
# 1. 定义你 FeaturePlot 使用的渐变色列表（作为颜色库）
# 这是一个连续渐变的颜色序列
rainbow_colors <-c("#FEE08B", "#F46D43", "#9E0142", "#E6F598", 
                   "#66C2A5", "#5E4FA2")
#FEE08B
# 2. 从你的 Seurat 对象中获取细胞类型的数量和名称
cell_types <- levels(cytotrace2_result_sce$geneSet.Type)
num_cell_types <- length(cell_types)

# 3. 为每种细胞类型选择颜色
# 我们可以使用 scales::viridis_pal() 或 RColorBrewer::brewer.pal() 
# 或者直接从 rainbow_colors 列表中采样

# 假设你只有 6 种或更少的细胞类型，我们直接使用 rainbow_colors 作为离散颜色
# 如果细胞类型 > 6，请考虑使用其他离散调色板，如 scales::hue_pal() 或 RColorBrewer::brewer.pal(8, "Dark2")
if (num_cell_types <= length(rainbow_colors)) {
  custom_colors <- setNames(rainbow_colors[1:num_cell_types], cell_types)
} else {
  # 如果细胞类型太多，我们使用一个高对比度的离散调色板
  custom_colors <- scales::hue_pal()(num_cell_types) 
  custom_colors <- setNames(custom_colors, cell_types)
}


p5_themed_colored <- p5 + 
  
  # --- 应用自定义离散颜色 ---
  scale_color_manual(values = custom_colors) + # 使用 scale_color_manual 设置离散颜色
  
  # --- 应用与 FeaturePlot 相同的主题和标签 ---
  xlab("UMAP1") + ylab("UMAP2") +
  ggtitle("Cell Phenotype") + 
  
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16,
                                  face = "bold", hjust = 0.5,
                                  margin = margin(b = 20))) +
  
  theme(aspect.ratio = 1)


# 打印查看最终效果
pdf(file = "CytoTrace_umap2.pdf", width = 4, height = 4)
p5_themed_colored
dev.off()



library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# 1. 定义参数
Y_MIN <- 0.0
Y_MAX <- 0.4
VIOLIN_ALPHA <- 0.85

# 2. 准备数据
plot_data <- cytotrace2_result_sce@meta.data

# 3. 创建精致的颜色方案
n_colors <- length(unique(plot_data$geneSet.Type))
color_palette <- c("#5E4FA2","#FEE08B","#9E0142")

# 4. 构建精致的小提琴箱线图
p_violin_box <- ggplot(plot_data, 
                       aes(x = geneSet.Type, y = CytoTRACE2_Score)) +
  
  # --- 小提琴图层 ---
  geom_violin(aes(fill = geneSet.Type), 
              alpha = VIOLIN_ALPHA,
              color = "white",
              linewidth = 0.8,
              trim = TRUE,
              scale = "width"
  ) +
  
  # --- 箱线图层 ---
  geom_boxplot(width = 0.15,
               outlier.shape = NA,
               color = "gray25",
               size = 0.5,
               fill = "white",
               alpha = 0.9
  ) +
  
  # --- 颜色和填充 ---
  scale_fill_manual(values = color_palette) +
  
  # --- 坐标轴和标签 ---
  labs(x = "", 
       y = "Predicted ordering\nby CytoTRACE") +  # 使用\n换行符
  
  coord_cartesian(ylim = c(Y_MIN, Y_MAX)) +
  
  # --- 精细主题调整 ---
  theme_pubr(base_size = 12) +
  theme(
    # X轴文字 - 改为纯黑色
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, 
                               color = "black", face = "plain"),
    # Y轴文字 - 改为纯黑色
    axis.text.y = element_text(size = 16, color = "black"),
    
    # Y轴标题 - 居中显示
    axis.title.y = element_text(size = 16, color = "black", face = "plain",
                                margin = margin(r = 10),
                                hjust = 0.5, vjust = 0.5),  # 水平和垂直居中
    
    # X轴和Y轴线 - 改为纯黑色并加粗
    axis.line = element_line(color = "black", linewidth = 0.8),  # 坐标轴线加黑
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 刻度线加黑
    
    # 图例调整
    legend.title = element_blank(),
    legend.text = element_text(size = 18, color = "black"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    
    # 面板和背景
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    
    # 整体绘图区域
    plot.background = element_rect(fill = "white", color = NA),
    
    # 移除所有标题相关的边距
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  ) +
  
  # 确保y轴从0开始
  expand_limits(y = 0)

# 显示图形
print(p_violin_box)

# 可选：保存高质量图片
ggsave("cytotrace_ordering_plot.pdf", p_violin_box, 
       width = 6, height = 4, dpi = 300, bg = "white")











library(ggstatsplot)
options(timeout = 10000)
remotes::install_github("IndrajeetPatil/ggstatsplot", force = TRUE)
install.packages("ggstatsplot")
library(ggplot2)
library(dplyr)
cytotrace2_result_sce@meta.data$CytoTRACE2_Score
# 1. 提取数据并创建数据框
# 假设两个Seurat对象的细胞顺序相同
plot_data <- data.frame(
  celltype = cytotrace2_result_sce@meta.data[["CytoTRACE2_Score"]],
  AUC = sce3@meta.data[["STING"]]
)

plot_data <- data.frame(
  celltype = cytotrace2_result_sce@meta.data[["CytoTRACE2_Relative"]],
  AUC = cytotrace2_result_sce@meta.data[["CytoTRACE2_Score"]]
)
# 2. 移除可能存在的NA值
plot_data <- plot_data %>% 
  filter(!is.na(celltype) & !is.na(AUC))

write.csv(plot_data,"plot_data2.csv")
# 3. 创建统计散点图
ggstatsplot::ggscatterstats(
  data = plot_data,
  x = AUC,                    # x轴为AUC值
  y = celltype,               # y轴为细胞类型（分类变量）
  xlab = "AUC Score",         # x轴标签
  ylab = "Cell Type",         # y轴标签
  title = "Relationship between AUC Score and Cell Type", # 图表标题
  messages = FALSE            # 不显示额外消息
)