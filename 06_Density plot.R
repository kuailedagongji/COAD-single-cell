library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(GSEABase)
library(harmony)
library(plyr) # 用于重命名



#yjsl1=readRDS("小胶质评分后.rds")
yjsl1=yjsl
# 1. 确保 ggplot2 已加载
library(ggplot2)
Q1_STING <- quantile(yjsl@meta.data$STING, probs = 0.25, na.rm = TRUE)
Q3_STING <- quantile(yjsl@meta.data$STING, probs = 0.75, na.rm = TRUE)

# 2. 准备绘图（p1 包含了基础图层）
p1 <- ggplot(yjsl1@meta.data, aes(x = STING)) +
  
  # 绘制灰色的“密度”柱状图
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 50,
                 fill = "lightgrey", 
                 color = "grey",
                 alpha = 0.7,           # 增加透明度使外观更柔和
                 linewidth = 0.3) +     # 细化边框线条
  
  # 绘制粉色的“密度”曲线
  geom_density(color = "magenta", 
               linetype = "dashed", 
               linewidth = 1,         # 使用linewidth代替size
               alpha = 0.8) +           # 增加透明度
  
  # 添加 Q1 和 Q3 的红色垂直分割线
  geom_vline(xintercept = c(Q1_STING, Q3_STING), 
             color = "red", 
             linetype = "dashed", 
             linewidth = 1,           # 使用linewidth代替size
             alpha = 0.8) +             # 增加透明度
  
  # 设置标题和坐标轴标签
  labs(title = "Score Distribution",
       x = "Score",
       y = "Density") +
  
  # 自定义x轴范围（请根据您的数据调整范围）
  coord_cartesian(xlim = c(0, 6)) +     # 示例范围，请根据实际情况调整
  
  # 使用经典主题并自定义
  theme_classic() +
  
  # 自定义主题细节
  theme(
    # 所有文本设为黑色
    text = element_text(color = "black"),
    
    # 坐标轴标题 - 加大但不加粗
    axis.title.x = element_text(size = 18, face = "plain", color = "black"),
    axis.title.y = element_text(size = 18, face = "plain", color = "black"),
    
    # 坐标轴刻度标签 - 加大但不加粗
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    
    # 标题样式
    plot.title = element_text(size = 16, face = "plain", color = "black", hjust = 0.5),
    
    # 坐标轴线
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

# 打印绘图
print(p1)
ggsave("Sorce分布密度图.pdf", p1, 
       width = 7, height = 4, dpi = 300, bg = "white")
