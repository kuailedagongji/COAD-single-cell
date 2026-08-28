# =============================================================================
# R包自动安装脚本
# 功能：自动检查并安装单细胞数据分析所需的各类R包
# 包括CRAN、Bioconductor和GitHub三个来源的包
# =============================================================================

# -----------------------------------------------------------------------------
# 第一部分：从CRAN安装包
# CRAN (Comprehensive R Archive Network) - R官方包仓库
# -----------------------------------------------------------------------------

# 定义需要从CRAN安装的包列表
cran.packages <- c(
  "aplot",          # 图形组合和排版
  "BiocManager",    # Bioconductor包管理器
  "data.table",     # 高效数据处理
  "devtools",       # 开发工具，用于安装GitHub包
  "doParallel",     # 并行计算支持
  "doRNG",          # 可重复的并行随机数生成
  "dplyr",          # 数据操作和转换
  "ggfun",          # ggplot2功能扩展
  "gghalves",       # 半半图（如半小提琴图+半箱线图）
  "ggplot2",        # 图形语法绘图系统
  "ggplotify",      # 将图形转换为grob对象
  "ggridges",       # 山脊图（密度分布图）
  "ggsci",          # 科学期刊配色方案
  "irlba",          # 大规模矩阵的快速SVD计算
  "magrittr",       # 管道操作符 %>%
  "Matrix",         # 稀疏矩阵处理
  "msigdbr",        # MSigDB基因集数据库接口
  "pagoda2",        # 单细胞RNA-seq分析工具
  "pointr",         # 指针操作工具
  "purrr",          # 函数式编程工具
  "RcppML",         # 矩阵分解机器学习算法
  "readr",          # 快速数据读取
  "reshape2",       # 数据重塑（宽转长等）
  "reticulate",     # Python与R接口
  "rlang",          # 高级R元编程
  "RMTstat",        # 随机矩阵理论统计
  "RobustRankAggreg", # 稳健排名聚合算法
  "roxygen2",       # 文档生成工具
  "Seurat",         # 单细胞数据分析主流工具
  "SeuratObject",   # Seurat对象类定义
  "stringr",        # 字符串处理
  "tibble",         # 现代数据框
  "tidyr",          # 数据整理
  "tidyselect",     # 选择变量辅助函数
  "tidytree",       # 树结构数据处理
  "VAM"             # 基因集变异分析
)

# 循环检查并安装CRAN包
for (i in cran.packages) {
  # 检查包是否已安装，quietly=TRUE表示不显示警告信息
  if (!requireNamespace(i, quietly = TRUE)) {
    # 如果未安装，则安装该包
    install.packages(i, ask = F, update = F)
  }
}

# -----------------------------------------------------------------------------
# 第二部分：从Bioconductor安装包
# Bioconductor - 生物信息学专用包仓库
# -----------------------------------------------------------------------------

# 定义需要从Bioconductor安装的包列表
bioconductor.packages <- c(
  "AUCell",               # 基于AUC的细胞类型鉴定
  "BiocParallel",         # Bioconductor并行计算
  "ComplexHeatmap",       # 复杂热图绘制
  "decoupleR",            # 网络生物学分析
  "fgsea",                # 快速基因集富集分析
  "ggtree",               # 系统发育树可视化
  "GSEABase",             # 基因集富集分析基础
  "GSVA",                 # 基因集变异分析
  "Nebulosa",             # 单细胞数据密度可视化
  "scde",                 # 单细胞差异表达
  "singscore",            # 单样本基因集评分
  "SummarizedExperiment", # 基因组实验数据容器
  "UCell",                # 单细胞基因集评分
  "viper",                # 蛋白活性推断
  "sparseMatrixStats"     # 稀疏矩阵统计
)

# 循环检查并安装Bioconductor包
for (i in bioconductor.packages) {
  # 检查包是否已安装
  if (!requireNamespace(i, quietly = TRUE)) {
    # 注意：这里原代码有误，应该使用BiocManager::install()
    # 修正后的代码应该是：
     BiocManager::install(i, ask = F, update = F)
    #install.packages(i, ask = F, update = F)  # 这行需要修正
  }
}

# -----------------------------------------------------------------------------
# 第三部分：从GitHub安装包
# GitHub - 开发中的最新版本包
# -----------------------------------------------------------------------------

# 检查并安装irGSEA包（单细胞基因集富集分析工具）
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  # 从GitHub安装开发版本
  devtools::install_github("chuiqin/irGSEA", force = T)
}


library("irGSEA")
# =============================================================================
# 使用说明：
# 1. 首次运行此脚本将安装所有缺失的包
# 2. 如果某个包安装失败，可以单独安装该包
# 3. Bioconductor部分需要先安装BiocManager包
# 4. 建议定期更新包以获得最新功能
# =============================================================================

# 注意：在实际使用中，Bioconductor包的安装应该使用：
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(i, ask = FALSE, update = FALSE)




#install.packages("rlang", type = "binary")
#BiocManager::install("SeuratData")
#remotes::install_github("satijalab/seurat-data")

# options(timeout = 10000)
# devtools::install_github("samuel-marsh/scCustomize")

rm(list=ls())
library(irGSEA)
library(Seurat)
library(SeuratData)
library(RcppML)
#load("Macrophage.rda")
yjsl=readRDS("af3_注释后.rds")
#Macrophage
sc_dataset<-yjsl
rm(yjsl)

#Check
sc_dataset$final_annotation 
UMAP_celltype<-DimPlot(sc_dataset,reduction="umap",
                         group.by="final_annotation",label=T);UMAP_celltype
Idents(sc_dataset)<-sc_dataset$final_annotation
scCustomize::DimPlot_scCustom(sc_dataset,figure_plot=TRUE)




sc_dataset<-SeuratObject::UpdateSeuratObject(sc_dataset)
sc_dataset2<-CreateSeuratObject(counts=CreateAssay5Object(GetAssayData(sc_dataset,
                                                                           assay="RNA",
                                                                           slot="counts")),
                                  meta.data=sc_dataset[[]])
STING<- data.table::fread("script/基因集.txt")
{
  gene <- STING$Gene
  ##格式需要一个list
  genes=list(gene) 
}
rm(sc_dataset)
sc_dataset2<-NormalizeData(sc_dataset2)
options(future.globals.maxSize=100000*1024^5)
sc_dataset2<-irGSEA.score(object=sc_dataset2,assay="RNA",
                            slot="data",seeds=123,
                            #ncores=1,
                            min.cells=3,min.feature=0,
                            custom=T,geneset=genes,msigdb=F,
                            species="Homo sapiens",
                            category="H",
                            subcategory=NULL,
                            geneid="symbol",
                            method=c("AUCell","UCell","singscore",
                                     "JASMINE","VAM","AddModuleScore","scSE"),
                            aucell.MaxRank=NULL,
                            ucell.MaxRank=NULL,
                            kcdf='Gaussian')
saveRDS(sc_dataset2,"sc_dataset2-1.rds")


#setwd("D:/R/任务/ZX分析/肾缺血再灌注/单细胞/")   

#sc_dataset2<-readRDS("sc_dataset2.rds")
result.dge<-irGSEA.integrate(object=sc_dataset2,
                               group.by="newType",
                               method=c("AUCell","UCell","singscore",
                                        "JASMINE","VAM","AddModuleScore","scSE"))



Idents(sc_dataset2)<-sc_dataset2$final_annotation
vlnplot<-irGSEA.vlnplot(object=sc_dataset2,
                          method=c("AUCell","UCell","singscore",
                                   "JASMINE","VAM","AddModuleScore","scSE"),
                          show.geneset="geneset1")
pdf(file = "打分方法汇总1.pdf", width = 8, height = 4)
vlnplot

dev.off()
table(sc_dataset2$newType)
upset.plot<-irGSEA.upset(object=result.dge,
                           method="AUCell")
upset.plot


ridgeplot<-irGSEA.ridgeplot(object=sc_dataset2,
                              method="AUCell",
                              show.geneset="geneset1")

ridgeplot
dev.off()
densityheatmap<-irGSEA.densityheatmap(object=sc_dataset2,
                                        method="AUCell",
                                        show.geneset="geneset1")
pdf(file = "多方法打分密度分布热图.pdf", width = 6, height = 5)
densityheatmap
dev.off()
#sc_dataset2@assays$AUCell@scale.data

median(sc_dataset2@assays$AUCell@scale.data)


library(Seurat)
library(dplyr)
score_columns
# 从metadata中提取所有评分列
score_columns <- c("singscore_score", "ssgsea_score", "VAM_score",
                     "scSE_score")

# 确保这些列都存在
available_scores <- score_columns[score_columns %in% colnames(sc_dataset2@meta.data)]
print(paste("可用的评分方法:", length(available_scores)))
print(available_scores)

# 提取原始评分矩阵
raw_scores <- sc_dataset2@meta.data[, available_scores, drop = FALSE]

# 检查数据
print("原始评分统计:")
print(summary(raw_scores))

# 步骤1: Z-score标准化（使均值为0，标准差为1）
normalized_scores <- scale(raw_scores)

# 步骤2: Min-Max归一化到[0,1]范围
min_max_normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

final_scores <- apply(normalized_scores, 2, min_max_normalize)

# 步骤3: 按行求和得到总评分
combined_score <- rowSums(final_scores, na.rm = TRUE)

# 将总评分添加到metadata
sc_dataset2@meta.data$combined_gene_set_score <- combined_score

# 检查结果
print("组合评分统计:")
print(summary(combined_score))





library(ggridges)
library(ggplot2)
yjsl$final_annotation

# 创建浅色系调色板
light_colors <- c("#CD5C5C", "#4169E1", "#32CD32", "#FFD700", "#8A2BE2",
                 "#DAA520", "#20B2AA", "#FF8C00", "#4682B4", "#D2B48C",
                 "#9370DB", "#40E0D0", "#90EE90", "#DB7093", "#9370DB")
# 修改山脊图 - 只按cellType分组，使用浅色系
ridgeplot_simple <- ggplot(sc_dataset2@meta.data, 
                           aes(x = combined_gene_set_score, 
                               y = final_annotation, 
                               fill = final_annotation)) +  # 按cellType填充颜色
  geom_density_ridges(
    alpha = 0.9,
    scale = 1.5,
    rel_min_height = 0.01,
    jittered_points = FALSE,
    color = "black",  # 改为黑色边框
    size = 0.9        # 稍微加大边框线粗细
  ) +
  scale_fill_manual(
    values = light_colors[1:length(unique(sc_dataset2@meta.data$cellType))],
    name = "Cell Type"
  ) +
  theme_ridges() +
  labs(
    x = "Score",
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(size = 18, face = "plain"),  # 加大x轴标签字体
    axis.text.y = element_text(size = 18, face = "plain"),  # 加大y轴标签字体
    axis.title.x = element_text(size = 18, face = "plain"), # 加大x轴标题字体
    axis.title.y = element_text(size = 18, face = "plain"), # 加大y轴标题字体
    legend.position = "none",  # 隐藏图例，因为y轴已经显示了细胞类型
    
    # 添加坐标轴线
    axis.line = element_line(color = "black", size = 0.8),  # 坐标轴线
    axis.ticks = element_line(color = "black", size = 0.6), # 坐标轴刻度线
    panel.grid.major = element_line(color = "white", size = 0.3)  # 保持网格线但调浅
  )

print(ridgeplot_simple)

# 保存图形
ggsave("ridgeplot_simple1.pdf", ridgeplot_simple, 
       width = 5, height = 4, dpi = 300)


saveRDS(yjsl,"yjsl评分后.rds")



scores_to_transfer <- sc_dataset2$IRI1
names(scores_to_transfer) <- colnames(sc_dataset2)

# 3. 将评分添加到 yjsl 中
# AddMetaData 会自动匹配名字，如果 yjsl 中有 sc_dataset2 没有的细胞，会填入 NA
yjsl <- AddMetaData(object = yjsl1, 
                    metadata = scores_to_transfer, 
                    col.name = "IR_Score")

# 4. (可选) 检查是否转移成功
head(yjsl$combined_gene_set_score)

# 5. (可选) 可视化检查
FeaturePlot(yjsl, features = "IR_Score")
