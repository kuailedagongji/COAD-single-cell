########## 展示每个通路与每个细胞群的通讯关系 ##########
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# 方法1：使用netAnalysis_signalingRole_heatmap - 显示每个细胞群在每个通路中的信号角色
pdf("通路-细胞群通讯热图1.pdf", width = 12, height = 8)
netAnalysis_signalingRole_heatmap(cellchat, 
                                  pattern = "all",
                                  width = 12, 
                                  height = 8,
                                  color.heatmap = "RdBu",
                                  font.size = 8)
dev.off()

# 方法2：使用netAnalysis_signalingRole_heatmap分别显示outgoing和incoming信号
pdf("细胞群通讯热图1.pdf", width = 5, height = 5)
netAnalysis_signalingRole_heatmap(cellchat, 
                                  pattern = "outgoing",
                                  width = 8, 
                                  height = 8,
                                  color.heatmap = "GnBu",
                                  font.size = 8)
dev.off()

pdf("细胞群通讯热图2.pdf", width =5 , height = 5)
netAnalysis_signalingRole_heatmap(cellchat, 
                                  pattern = "incoming",
                                  width = 8, 
                                  height = 8,
                                  color.heatmap = "GnBu",
                                  font.size = 8)
dev.off()


# 方法3：自定义热图 - 展示每个通路在每个细胞群中的总信号强度
# 计算每个细胞群在每个通路中的信号强度
signaling_heatmap_data <- cellchat@netP[["pathways"]] %>% 
  map_dfr(function(pathway){
    # 获取该通路的信号矩阵
    prob_matrix <- cellchat@netP$prob[,,pathway]
    # 计算每个细胞群的总信号强度（发送+接收）
    rowSums(prob_matrix) + colSums(prob_matrix)
  }, .id = "pathway")

# 整理数据为热图格式
heatmap_matrix <- signaling_heatmap_data %>%
  pivot_longer(cols = -pathway, names_to = "celltype", values_to = "strength") %>%
  pivot_wider(names_from = celltype, values_from = strength) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# 绘制热图
pdf("通路-细胞群通讯热图4_自定义.pdf", width = 10, height = 8)
pheatmap::pheatmap(heatmap_matrix,
                   color = colorRampPalette(c("white", "yellow", "red"))(100),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   scale = "row",  # 按行标准化
                   main = "Pathway-Celltype Communication Strength",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   angle_col = 45)
dev.off()

# 方法5：使用dot plot展示通路活性（推荐，更直观）
pdf("通路-细胞群通讯点图.pdf", width = 14, height = 10)
# 计算每个细胞群在每个通路中的相对活性
pathway_activity <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


pdf("通路-细胞群通讯点图.pdf", width =4.5, height = 3.5)
# 绘制点图展示通路活性
netAnalysis_signalingRole_scatter(cellchat, 
                                  signaling = cellchat@netP$pathways,
                                  font.size = 11,
                                  dot.size = c(2, 6))
dev.off()

# 方法6：通路活性热图 - 展示每个细胞群对每个通路的贡献
pdf("通路活性热图.pdf", width = 12, height = 10)
netVisual_heatmap(cellchat, 
                  signaling = cellchat@netP$pathways,
                  color.heatmap = "Reds",
                  width = 12,
                  height = 10,
                  font.size = 8,
                  font.size.title = 10)
dev.off()

# 方法7：重点展示特定细胞群相关的通路
# 选择您感兴趣的细胞群
interesting_celltypes <- c("T cells", "Macrophages", "Endothelial cells") # 请替换为您的实际细胞类型

# 提取与这些细胞群相关的通路
celltype_specific_pathways <- identifyEnrichedPaths(cellchat, 
                                                    celltypes = interesting_celltypes)

# 绘制特定细胞群相关的通路热图
if(length(celltype_specific_pathways) > 0){
  pdf("特定细胞群相关通路热图.pdf", width = 10, height = 8)
  netVisual_heatmap(cellchat, 
                    signaling = celltype_specific_pathways,
                    color.heatmap = "Purples",
                    width = 10,
                    height = 8)
  dev.off()
}

# 方法8：使用rankNet函数排序并展示最重要的通路
pdf("通路重要性排序图.pdf", width = 12, height = 8)
rankNet(cellchat, 
        mode = "comparison", 
        stacked = T, 
        font.size = 8,
        sources.use = 1:length(levels(cellchat@idents)), # 所有细胞群
        targets.use = 1:length(levels(cellchat@idents))) # 所有细胞群
dev.off()

print("通路-细胞群通讯关系分析完成！")








########## 批量可视化贝壳图（所有通路） ##########
########## 批量生成细胞群间通讯网络图（基础版） ##########

# 创建保存文件夹
if(!dir.exists("细胞群通讯网络图")) dir.create("细胞群通讯网络图")

# 获取通信矩阵和细胞群大小
mat <- cellchat@net$count  # 使用count或weight根据需要选择
groupSize <- as.numeric(table(cellchat@idents))

# 批量生成每个细胞群作为信号源的网络图
pdf("细胞群通讯网络图/所有细胞群通讯网络图1.pdf", width = 6, height = 6)
par(mfrow = c(3, 4), xpd = TRUE, mar = c(2,2,3,2)) # 调整布局根据细胞群数量

for (i in 1:nrow(mat)) {
  # 创建只包含当前细胞群发出信号的矩阵
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  # 绘制网络图
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   edge.weight.max = max(mat), 
                   title.name = paste0(rownames(mat)[i]),
                   vertex.size = 13,
                   vertex.label.cex = 0.6,
                   edge.width.max = 8)
}
dev.off()

cat("已生成所有细胞群通讯网络图组合PDF\n")
