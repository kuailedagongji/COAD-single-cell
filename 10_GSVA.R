library(Seurat)
library(GSVA)
library(clusterProfiler)
library(pheatmap)

# 加载预处理过的单细胞数据对象
yjsl1<-readRDS("af3_注释后.rds")
#load("Macrophage.rda")
# 查看不同细胞类型的分布情况
Idents(yjsl1) <- yjsl1$group
table(Idents(yjsl1))

{
# 计算每种细胞类型的平均基因表达量
exp = AverageExpression(yjsl1)[[1]]  # 使用 AverageExpression 函数计算每个细胞群的平均基因表达矩阵
# 将结果转换为矩阵格式
exp = as.matrix(exp)  
# 去掉所有表达量为 0 的基因
exp = exp[rowSums(exp) > 0, ]     
}
#看一下数据情况
#exp[1:4, 1:4]                      

# 读取基因集数据库文件
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
{
genes = read.gmt("D:/R/h.all.v2023.1.Hs.symbols.gmt")[,c(2,1)]
# 将基因与基因集的配对信息转换为一个列表
list = unstack(genes)
}

# Score <- gsva(exprData = exp, 
#               gset.idx.list = list, 
#               method = "gsva", 
#               kcdf = "Gaussian", 
#               max.diff = TRUE)


#BiocManager::install("GSVA",force = TRUE)
# 2. 重新运行你的代码
# 注意：你的变量 'exp' 应该是表达矩阵，'list' 应该是基因集列表
#Score = GSVA::gsva(GSVA::gsvaParam(exp, list, maxDiff = TRUE))
Score = GSVA::gsva(expr = exp, 
                   gset.idx.list = list, 
                   mx.diff=TRUE)

write.csv(Score,"Score1.csv")
Score<-read.csv("Score1.csv",row.names = 1)

# 1. 加载必要的包
library(pheatmap)

# 2. 读取数据 (假设您的数据框叫 Score)
# 如果是从文件读取，请取消下面这行的注释:
# Score <- read.csv("Score.csv", row.names = 1)

# 3. 设置自定义颜色 (蓝-白-红)
my_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# 4. 绘图并直接导出
# 使用 pdf() 或 png() 并不是必须的，pheatmap 自带 filename 参数可以直接保存
pheatmap(Score,
         # --- 核心设置 ---
         scale = "row",               # 按行标准化 (Z-score)
         color = my_colors,           # 自定义颜色
         
         # --- 布局与聚类 ---
         cluster_rows = TRUE,         # 是否聚类行
         cluster_cols = F,         # 是否聚类列
         angle_col = 45,              # 列名旋转 45 度 (防重叠)
         
         # --- 字体优化 (这里是您最需要的) ---
         fontsize = 11,               # 基础字体大小
         fontsize_row = 11,            # 行名(通路)字体大小，太挤就把这个调小
         fontsize_col = 14,           # 列名(细胞)字体大小，可以调大
         
         # --- 格子美化 ---
         border_color = "white",      # 格子边框颜色 (白色看起来更精致)
         cellwidth = NA,              # 自动调整，如果想固定宽度可设为数值(如 15)
         cellheight = NA,             # 自动调整，如果想固定高度可设为数值(如 10)
         
         # --- 导出设置 ---
         filename = "Heatmap.pdf",  # 直接保存为 PDF (矢量图，无限放大不糊)
         width = 12,                   # 图片宽度 (英寸)
         height = 10                  # 图片高度 (英寸)，根据通路数量适当拉长
)


library(tidyr)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# 1. 你的彩色调色板 (保持不变)
rainbow_colors <- c("#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", 
                    "#F46D43", "#9E0142")

# 2. 将 Score 矩阵转换为长格式数据框 (保持不变)
score_long_data <- sting_filtered_data %>%
  as.data.frame() %>%
  rownames_to_column(var = "Pathway") %>% 
  pivot_longer(
    cols = -Pathway, 
    names_to = "Sample_Group", 
    values_to = "GSVA_Score" 
  )

# ==============================================================================
# 3. 对通路名称进行换行处理 (修改核心：每 2 个单词换一次行)
# ==============================================================================
score_long_data$Pathway_wrapped <- sapply(score_long_data$Pathway, function(x) {
  # 1. 按下划线分割成单词
  words <- str_split(x, "_")[[1]]
  
  # 2. 创建分组索引：这里改成除以 2
  # 结果示例: 1, 1, 2, 2, 3, 3...
  groups <- ceiling(seq_along(words) / 2)
  
  # 3. 将每一组的单词用 "_" 重新连接
  lines <- tapply(words, groups, paste, collapse = "_")
  
  # 4. 将各组用换行符 "\n" 连接起来
  final_string <- paste(lines, collapse = "\n")
  
  return(final_string)
})

# 4. 转换因子 (保持不变)
score_long_data$Pathway_wrapped <- factor(score_long_data$Pathway_wrapped, 
                                          levels = unique(score_long_data$Pathway_wrapped))
score_long_data$Sample_Group <- factor(score_long_data$Sample_Group, 
                                       levels = unique(score_long_data$Sample_Group))

# 5. 创建气泡图
p_bubble <- ggplot(score_long_data, 
                   aes(x = Pathway_wrapped, y = Sample_Group, 
                       color = GSVA_Score, 
                       size = GSVA_Score)) +
  
  geom_point() +
  
  scale_color_gradientn(
    colors = rainbow_colors,
    name = "GSVA Score\n(High <--> Low)"
  ) +
  
  scale_size_continuous(range = c(1, 8)) + 
  
  theme_minimal() + 
  
  theme(
    # --- X 轴字体设置 ---
    # lineheight = 0.8: 行间距，防止多行文字太散
    # size = 12: 字号，如果觉得字太大挤在一起，可以改小成 10 或 11
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 12, color = "black", face = "bold",
                               lineheight = 0.8), 
    
    # Y 轴字体设置
    axis.text.y = element_text(angle = 0, 
                               size = 16, color = "black", face = "bold"),
    
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    axis.title = element_blank(),
    
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) 

# 6. 保存为PDF
# 由于变成了2个单词一行，标签会变得很高，建议稍微增加 height
pdf(file = "GSVA1_wide_2words.pdf", width = 16, height = 4) 
print(p_bubble)
dev.off()




library(ggplot2)
library(dplyr)
library(tibble)


library(tidyverse)

library(ggplot2)



# --- 1. 数据准备 (同前) ---

df_plot <- Score %>%
  
  as.data.frame() %>%
  
  rownames_to_column("Pathway") %>%
  
  mutate(diff = Normal - keloid) %>%
  
  arrange(diff)



# 美化通路名称 (去掉前缀和下划线)

df_plot$Pathway <- df_plot$Pathway %>%
  
  gsub("^KEGG_", "", .) %>%
  
  gsub("_", " ", .)



# 锁定顺序

df_plot$Pathway <- factor(df_plot$Pathway, levels = df_plot$Pathway)



# --- 2. 绘图 ---



p <- ggplot(df_plot, aes(x = Pathway, y = diff)) +
  
  
  
  # A. 0刻度参考线 (最先画，作为背景)
  
  # 用浅灰色虚线，不抢眼
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.8) +
  
  
  
  # B. 棒棒糖的“棒” (带渐变色)
  
  # 这里的 color 映射到 diff 数值，实现线条颜色的渐变
  
  geom_segment(aes(x = Pathway, xend = Pathway, y = 0, yend = diff, color = diff), 
               
               linewidth = 1) +  # 线条粗细设为 1，既不细得看不见，也不显粗糙
  
  
  
  # C. 棒棒糖的“糖” (带渐变色)
  
  # shape=21 允许填充颜色。外圈用白色 (color="white") 勾勒，增加精致感
  
  geom_point(aes(fill = diff), shape = 21, size = 5, color = "white", stroke = 0.6) +
  
  
  
  # D. 坐标轴翻转
  
  coord_flip() +
  
  
  
  # E. 颜色渐变设置 (关键)
  
  # 使用 gradient2 设置三点颜色：低(负值)=冷色, 中(0)=白色, 高(正值)=暖色
  
  # 这里的颜色推荐：#3288BD (深蓝) -> White -> #D53E4F (深红)，这是经典的 Spectral 配色
  
  scale_color_gradient2(low = "#3288BD", mid = "white", high = "#D53E4F", midpoint = 0) +
  
  scale_fill_gradient2(low = "#3288BD", mid = "white", high = "#D53E4F", midpoint = 0) +
  
  
  
  # F. 标签
  
  labs(x = NULL, y = "GSVA Score Difference") +
  
  
  
  # G. 主题去框与极简风 (The Cleanest Look)
  
  theme_classic() + # 基础去框
  
  
  
  theme(
    
    # 1. 去除多余线条
    
    axis.line.y = element_blank(), # 去掉 Y 轴的那根竖线 (只保留左侧文字)
    
    axis.line.x = element_line(color = "black", linewidth = 0.5), # 保留 X 轴底线
    
    axis.ticks.y = element_blank(), # 去掉 Y 轴刻度短线
    
    
    
    # 2. 字体调整
    
    text = element_text(color = "black", family = "sans"),
    
    axis.text.y = element_text(size = 11, color = "black"), # 通路名称
    
    axis.text.x = element_text(size = 10, color = "black"),
    
    
    
    # 3. 图例处理
    
    # 如果觉得颜色条图例多余（因为X轴已经有数值了），可以去掉。
    
    # 这里我保留了颜色条，但把它做得很细
    
    legend.position = "right", 
    
    legend.title = element_blank(), # 去掉图例标题
    
    legend.text = element_text(size = 8),
    
    legend.key.height = unit(1.5, "cm"), # 图例拉长一点好看
    
    legend.key.width = unit(0.3, "cm")
    
  )



# --- 3. 预览与保存 ---

print(p)



ggsave(filename = "Gradient_Lollipop.pdf", plot = p, width = 7, height = 7)
