rm(list = ls())  
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
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
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)

yjsl<-read_rds("af5_注释后.rds")
table(yjsl$final_annotation)
library(Seurat)
library(dplyr)

# 1. 安全备份：先把原本的注释备份一列，万一分错了还能找回
# 如果已经备份过可以跳过这一步
yjsl$final_annotation_bak <- yjsl$final_annotation

# 2. 提取 Malignant 细胞的 ID
# 确保你的 meta.data 里列名是 final_annotation，且其中有 "Malignant"
malignant_cells <- WhichCells(yjsl, expression = final_annotation == "Malignant")

# 3. 获取 MACC1 在这些细胞中的表达量
# 默认使用 data slot (归一化后的数据)，这比 counts 更适合做高低分组
# 如果你的 MACC1 表达量极低，可能需要检查一下
macc1_expr <- FetchData(yjsl, vars = "MACC1", cells = malignant_cells)

# 4. 计算中位数并定义分组
# 注意：在单细胞数据中，很多基因的中位数可能是 0。
# 如果 median 是 0，那么 "High" 就是所有表达 > 0 的细胞，"Low" 就是所有表达 = 0 的细胞。
median_val <- median(macc1_expr$MACC1)
print(paste0("Malignant 群体中 MACC1 的中位表达量为: ", median_val))

# 创建一个新的分组向量
# 如果表达量 > 中位数 -> High_MACC1
# 否则 (<= 中位数) -> Low_MACC1
new_labels <- ifelse(macc1_expr$MACC1 > median_val, "High_MACC1", "Low_MACC1")

# 将行名（细胞ID）赋予这个新向量，确保一一对应
names(new_labels) <- rownames(macc1_expr)

# 5. 更新 Seurat 对象的元数据
# 先将 final_annotation 转为字符型 (如果是 factor 因子型，直接加新水平可能会报错)
yjsl$final_annotation <- as.character(yjsl$final_annotation)

# 仅替换 Malignant 细胞的注释
yjsl$final_annotation[names(new_labels)] <- new_labels

# (可选) 将其转回 factor，如果你有固定的顺序需求可以在这里指定 levels
yjsl$final_annotation <- factor(yjsl$final_annotation)

# 6. 检查结果
print("替换后的细胞群统计：")
table(yjsl$final_annotation)

# 7. 更新 Idents (设为当前默认身份)
Idents(yjsl) <- "final_annotation"
#setwd("D:/BaiduNetdiskDownload/!免费的午餐单细胞代码分析/5、细胞通讯//")  
#load("yjsl_Type.rda")

head(yjsl@meta.data)
table(yjsl$final_annotation)
yjsl$final_annotation

N=length(colnames(yjsl))/10# 获取总细胞数的10%

{
  N=round(N)# 对细胞数量进行四舍五入
  # 从数据集中随机抽取N个细胞
  yjsl<-yjsl[,sample(x=colnames(yjsl),size = N,replace=F)]
}
# #创建cellchat对象
cellchat = createCellChat(object = yjsl,
                          group.by = "final_annotation")#通过 group.by 定义分组


#展示以下现在的细胞分组
levels(cellchat@idents) 
# #细胞亚群各组数量
# group <- as.numeric(table(cellchat@idents))

#设置配体受体交互数据库 
CellChatDB <- CellChatDB.human #如果是老鼠的话使用内置“CellChatDB.mouse”数据

showDatabaseCategory(CellChatDB)
#人的数据包括61.8%的旁分泌/自分泌信号互作、
#21.7%的细胞外基质(ECM)受体互作
#16.5%的细胞-细胞通讯互作


unique(CellChatDB$interaction$annotation)

#如果想用全部的用于cellchat分析，不进行subsetDB，直接指定cellchat@DB <- CellChatDB 即可。
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# # set the used database in the object
# cellchat@DB <- CellChatDB.use
# use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use 


######################对表达数据进行预处理######################
rm(yjsl)
##This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)
#根据配置设置
if (exists("FUN") && is.function(FUN) && object.size(FUN) > 100000) {
  message("Found and removing the giant FUN object (", format(object.size(FUN), units = "auto"), ")...")
  rm(FUN)
}

# 1. 设置并行计划
future::plan("multisession", workers = 4)
future::plan("sequential")
options(future.globals.maxSize = 5 * 1024^3)
#devtools::install_github('immunogenomics/presto')
{
  # 识别过表达基因(很慢)
  cellchat <- identifyOverExpressedGenes(cellchat)
  # 识别过表达配体受体对
  cellchat <- identifyOverExpressedInteractions(cellchat)
}
# {


{
  # 直接进行细胞通讯推断（不需要 projectData）
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  
  # 过滤低可信度的通讯
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  # 计算通讯网络
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 聚合网络
  cellchat <- aggregateNet(cellchat)
}




#建议保存一下
save(cellchat,file = "cellchat1.rda")
load("cellchat.rda")


##提取 保存结果
df.net <- subsetCommunication(cellchat)
head(df.net)
write.csv(df.net,"df.net.csv")

#df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
##查看细胞通讯分群###
levels(cellchat@idents)

########3种方式建立网络########
##通过序号
df.net1 <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
head(df.net1)
##通过名字
df.net2 <- subsetCommunication(cellchat, sources.use = c("T cells"), targets.use = c("Mast cells" ,"Mast cells")) 
head(df.net2)
##通过通路
df.net3 <- subsetCommunication(cellchat, signaling = c("EGF"))
head(df.net3)


############计算cell-cell communication#####
#使用computeCommunProbPathway计算每个信号通路的所有配体-受体相互作用的通信结果，结存存放在net 和 netP中 
cellchat <- computeCommunProbPathway(cellchat)
##使用aggregateNet计算细胞类型间整合的细胞通讯结果
cellchat <- aggregateNet(cellchat)


##########快乐的可视化过程##########

##########1、celltype之间通讯结果##########
##########通讯次数（左）通讯强度(右)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")


# 导出高质量PDF
pdf("细胞通讯1图.pdf", width = 6, height = 6)
#par(mfrow = c(1,2), xpd=TRUE, mar = c(2,2,4,2))

netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Number of interactions",
                 vertex.size = 18,
                 vertex.label.cex = 1.3,
                 edge.width.max = 8)
dev.off()
pdf("细胞通讯2图.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge = F, 
                 title.name = "Interaction weights/strength",
                 vertex.size = 18,
                 vertex.label.cex = 1.3,
                 
                 edge.width.max = 8)

dev.off()
#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。
#发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。
#右图：互作的概率或者强度值（强度就是概率值相加）

##换一种可视化方式
p3 <- netVisual_heatmap(cellchat)
p3
p4 <- netVisual_heatmap(cellchat, measure = "weight")
p4
p3 + p4

##########分别展示
#根据celltype的个数，灵活调整mfrow = c(2,3) 参数
mat <- cellchat@net$weight

###另一种
{
  mat <- cellchat@net$count
  par(mfrow = c(2,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

##########单个信号通路可视化

#展示当前有哪些通路结果
cellchat@netP$pathways
#选择感兴趣的进行展示
pathways.show <- c("TGFb") 

#查看当前的celltype顺序
levels(cellchat@idents) 
#通过vertex.receiver指定target 的细胞类型
##左边是你选定的细胞群，右边是剩下的细胞群
vertex.receiver = c(1,2) 

##层次图
netVisual_aggregate(cellchat, signaling = "TGFb", 
                    vertex.receiver = vertex.receiver,layout="hierarchy")
#在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。

##圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="TGFb", layout = "circle")

##和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="TGFb", layout = "chord", vertex.size = groupSize)

##热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "TGFb", color.heatmap = "Reds")

##########气泡图##############################
#通过sources.use 和 targets.use指定定受体-配体
levels(cellchat@idents)
pdf("细胞群通讯图3.pdf", width =5.5, height = 7)
netVisual_bubble(cellchat, sources.use = c(2,6,7), 
                 targets.use = c(1,3,4,5,8,9,10,11,12), remove.isolate = FALSE)

dev.off()


# 1. 备份原始的 p-value 数据（以防后续分析还需要用到原始 P 值）
original_pval <- cellchat@net$pval

# 2. 提取原始 P 值并进行多重假设检验校正
# 注意：cellchat@net$pval 是一个三维数组，p.adjust 需要对向量操作，因此需要先打平再还原
pvals_vector <- as.vector(cellchat@net$pval)

# 这里以 "fdr" (即 "BH") 方法为例，你也可以换成 "bonferroni" 等其他方法
pvals_adj_vector <- p.adjust(pvals_vector, method = "fdr") 

# 3. 还原维度结构并覆盖回 cellchat 对象
pvals_adj_array <- array(pvals_adj_vector, 
                         dim = dim(original_pval), 
                         dimnames = dimnames(original_pval))
cellchat@net$pval <- pvals_adj_array

# 4. 重新绘制气泡图 (此时图例中的 p-value 代表的即为 FDR 校正后的 p-value)
# 加载 ggplot2（如果之前没加载的话）
library(ggplot2)

# 4. 重新绘制气泡图并修改图例标题
library(ggplot2)

pdf("细胞群通讯图3_FDR校正.pdf", width = 5.5, height = 7)

p <- netVisual_bubble(cellchat, 
                      sources.use = c(2,6,7), 
                      targets.use = c(1,3,4,5,8,9,10,11,12), 
                      remove.isolate = FALSE)

# 【关键修改】：使用 guides 强制重写图例标题
p <- p + guides(size = guide_legend(title = "FDR")) 

# （可选）如果你想把颜色图例也一起修改，可以这样写：
# p <- p + guides(
#   size = guide_legend(title = "FDR"),
#   color = guide_colorbar(title = "Comm. Prob.")
# )

print(p)
dev.off()
# 5. (可选) 绘图完成后，如果需要进行其他基于原始 P 值的分析，可以将数据恢复
# cellchat@net$pval <- original_pval
##指定受体-配体细胞类型且指定通路
cellchat@netP$pathways 
netVisual_bubble(cellchat, sources.use = c(3,5), targets.use = c(1,2,4,6), 
                 signaling = c("TGFb","SPP1"), remove.isolate = FALSE)


#某条信号通路（如SPP1）的所有基因在细胞群中的表达情况展示
pdf("细胞群通讯图4.pdf", width =5, height = 4)
plotGeneExpression(cellchat, signaling = "TGFb")
dev.off()
