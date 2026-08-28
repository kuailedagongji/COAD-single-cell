library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
#devtools::install_github("smorabit/hdWGCNA", ref = "dev")
#BiocManager::install("enrichR")

#options(timeout = 10000)
{
  theme_set(theme_cowplot())
  
  set.seed(12345)
  
  enableWGCNAThreads(nThreads = 7)
}

yjsl=readRDS("上皮评分后.rds")
#load('yjsl_Type.rda')
#可视化一下看看
yjsl$geneSet.Type
DimPlot(yjsl, group.by='geneSet.Type', label=TRUE) 

seurat_obj<-yjsl

# 清除yjsl对象
rm(yjsl)
# 强制垃圾回收，释放内存
gc()


seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 5000 # 指定要查找的高变基因数量
)

# 2. 重新运行 SetupForWGCNA
# 它会自动使用上面步骤中存储在 VariableFeatures(seurat_obj) 里的 5000 个基因
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "variable",  # <--- 这样设置是正确的
  wgcna_name = "tutorial"
  # <--- top_n_genes 参数已移除
)


seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("geneSet.Type", "seurat_clusters"), 
  reduction = 'harmony', 
  k = 25,  # KNN参数
  max_shared = 10, # 两个metacell可共享的最大细胞数
  ident.group = 'geneSet.Type' 
)
options(future.globals.maxSize = 4000 * 1024^2)
#rm(yjsl1)
{
  # 标准化metacell表达矩阵
  seurat_obj <- NormalizeMetacells(seurat_obj)
}

#创建表达矩阵用于WGCNA
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "High", #挑选感兴趣的细胞类型
  group.by='geneSet.Type', 
  assay = 'RNA',
  slot = 'data' 
)

#选取软阈值
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'## 此外，还可选择“unsigned”或“signed hybrid”参数
)




plot_list <- PlotSoftPowers(seurat_obj)
pdf(file = "hdWGCNA1.pdf", width = 8, height = 6.5)
P1=wrap_plots(plot_list, ncol=2)
P1
dev.off()
{
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)
}
# 构造共表达网络

{
  seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power = 9,        # 你选定的软阈值
    setDatExpr = FALSE,
    overwrite_tom = TRUE,
    tom_name = 'High',
    
    # --- 关键修改参数 ---
    mergeCutHeight = 0.25,  # 【核心】合并阈值：0.25 代表合并相似度 > 0.75 的模块 (默认是 0.15)
    minModuleSize = 100,    # 【辅助】最小基因数：由默认的 50 提高到 100，强制小模块合并或丢弃
    deepSplit = 2           # 【微调】灵敏度：范围 0-4，值越小模块越少且越大 (默认通常是 2，可尝试改 1)
  )
}

pdf(file = "hdWGCNA2.pdf", width = 7, height = 4.5)
P2=PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
P2
#灰色模块表示无法组成任何共表达模块的基因，这些基因不纳入后续的分析
dev.off()
# 计算模块特征值ME: Module Eigengenes，模块特征基因（值）
#计算所有的MEs，比较耗时

#seurat_obj=readRDS("hdWGCNA_object1.rds")
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="seurat_clusters"# 根据其去批次
) #时间较长


# 计算模块连接性
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'geneSet.Type', 
  group_name = 'High'# 感兴趣的细胞类型的kME
)


pdf(file = "hdWGCNA3.pdf", width = 8, height = 4.5)

# 可视化每个模块中的基因的kME
P3=PlotKMEs(seurat_obj, ncol=4, n_hubs = 10)
P3
dev.off()
{
  # 获取模块分配表
  modules <- GetModules(seurat_obj)%>%subset(module!="grey")
}
#获取那些与每个模块高度连接的基因（hub genes）
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)

write.csv(hub_df,"hub_df基因.csv")

#保存一下hdWGCNA结果
saveRDS(seurat_obj, file = 'hdWGCNA_object1.rds')
#seurat_obj<-read_rds('hdWGCNA_object.rds')


#两种评分方式
#一种是seurat自带的AddModuleScore评分，另一种是UCell评分。
#利用hub genes进行打分
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'#AddModuleScore
)


###可视化
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs',# 可选择MEs、hMEs、scores、average
  order=TRUE 
)
pdf(file = "hdWGCNA4.pdf", width = 7, height = 6)
p4=wrap_plots(plot_list, ncol=3)
p4
dev.off()
# 模块相关性
pdf(file = "hdWGCNA5.pdf", width = 4.5, height = 5.5)
{
  ModuleCorrelogram(seurat_obj)
}
dev.off()
#下面两种图很重要，直接把模块与亚群，也就是表型相联系，
{
  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
}



# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'geneSet.Type')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p6 <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
pdf(file = "hdWGCNA6.pdf", width = 4.5, height = 5.5)
p6
dev.off()


p <- VlnPlot(
  seurat_obj,
  features = 'Blue',##换成自己的模块
  group.by = 'cellType_detailed',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')


# 确保 modules 和 hub_df 已经存在（根据您的脚本是存在的）
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = module_colors <- c("yellow", "brown", "turquoise", "blue", "black", "red")
)
# 1. 定义模块颜色（确保与 ResetModuleNames 一致）
module_colors <- c("yellow", "brown", "turquoise", "blue", "black", "red","green","pink","magenta","greenyellow","tan","purple")

# 2. 提取每个模块的前 50 个 Hub 基因
top_genes_df <- hub_df %>%
  filter(module != "grey") %>% # 排除灰色模块
  group_by(module) %>%
  slice_max(order_by = kME, n = 20, with_ties = FALSE) %>% # 选取kME最高的50个基因
  ungroup() %>%
  mutate(gene_id = as.character(gene_name)) # 确保基因名是字符

# 获取最终的基因列表
gene_list <- unique(top_genes_df$gene_id)
# 节点数据 (Nodes)
nodes <- top_genes_df %>%
  select(gene_id, module) %>%
  distinct() # 确保每个基因只出现一次
edges_list <- list()

for (mod in module_colors) {
  # 找出当前模块的所有基因
  mod_genes <- nodes %>% filter(module == mod) %>% pull(gene_id)
  
  if (length(mod_genes) > 1) {
    # 创建所有基因对的组合 (全连接)
    edges_mod <- t(combn(mod_genes, 2))
    
    # 转换为数据框
    edges_mod_df <- as.data.frame(edges_mod)
    colnames(edges_mod_df) <- c("from", "to")
    edges_mod_df$module <- mod # 标记这条边属于哪个模块
    
    edges_list[[mod]] <- edges_mod_df
  }
}

# 合并所有模块的边
edges <- bind_rows(edges_list)

library(ggraph)
library(igraph)
library(Cairo)

# 定义所有模块的颜色列表
module_colors <- c("yellow", "brown", "turquoise", "blue", "black", "red","green","pink","magenta","greenyellow","tan","purple")

# 确保颜色映射使用小写（与 module_colors 一致）
module_color_map <- c(
  "yellow" = "yellow", 
  "brown" = "brown", 
  "turquoise" = "turquoise", 
  "blue" = "blue", 
  "black" = "black", 
  "red" = "red",
  "green" = "green",
  "pink" = "pink",
  "magenta" = "magenta",
  "greenyellow"="greenyellow",
  "tan"="tan",
  "purple"="purple"
)

for (mod_name in module_colors) {
  
  # 1. 过滤当前模块的节点和边
  current_nodes <- nodes %>% filter(module == mod_name)
  current_edges <- edges %>% filter(module == mod_name)
  
  # 检查是否有足够的边来构建图
  if (nrow(current_edges) == 0) {
    cat(paste("跳过模块", mod_name, ": 边数不足。\n"))
    next
  }
  
  # 2. 从节点和边数据创建 igraph 对象
  g_mod <- graph_from_data_frame(d = current_edges, vertices = current_nodes, directed = FALSE)
  
  # 生成文件名
  file_name <- paste0(mod_name, "_HubGene_NetworkPlot.pdf")
  
  # 打开 PDF 文件
  cairo_pdf(filename = file_name, 
            width = 4, height = 4,
            family = "sans")
  
  # 4. 绘制网络图
  P_mod_network <- ggraph(g_mod, layout = 'fr') + 
    # 绘制边 (Edge)
    geom_edge_fan(
      alpha = 0.1, 
      color = "grey50"
    ) +
    # 绘制点 (Node)
    geom_node_point(
      aes(color = module), 
      size = 3.5 # 适当增大节点大小
    ) +
    # 绘制标签
    geom_node_text(
      aes(label = name), 
      repel = TRUE,
      size = 3.5, # 适当增大标签大小
      bg.color = "white" # 给标签加白色背景，增强可读性
    ) +
    # 设置颜色，只使用当前模块的颜色
    scale_color_manual(values = module_color_map[mod_name]) +
    # 设置标题
    labs(title = paste0("Top 20 Hub Genes Network: ", mod_name, " Module")) +
    # 应用主题并移除轴元素
    theme_cowplot() +
    theme_graph() +
    # 移除颜色图例，因为只有一个颜色
    guides(color = "none") 
  
  print(P_mod_network)
  
  # 关闭 PDF 文件
  dev.off()
  
  cat(paste("✅ 模块", mod_name, "的网络图已生成，文件名为:", file_name, "\n"))
}



