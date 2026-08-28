# # 设置 Bioconductor 版本为 3.8
# BiocManager::install(version = "3.16")
# BiocManager::install("monocle")
# # 尝试安装 Monocle 2.32.0
# install.packages("https://mirrors.westlake.edu.cn/CRAN/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz",repos = NULL,type="source")
#devtools::install_local("monocle1050/monocle.zip",upgrade = F)
{
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
}
# package.version("monocle")
# options(timeout = 10000)
#package.version("monocle")

# 加载已处理的数据
#load("Macrophage.rda")
yjsl=readRDS("上皮评分后.rds")


data <- as(as.matrix(yjsl@assays$RNA$counts), 'sparseMatrix')# 转换为稀疏矩阵格式
pd <- new('AnnotatedDataFrame', data = yjsl@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# 构建CellDataSet对象
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

######数据预处理######
#估计size factor和离散度，为后续的差异分析做铺垫
##如果用原文件，这一步时间会很长
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)


###细胞过滤
# 检测在至少0.1表达量以上的基因，这些基因会被保留用于后续分析
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 输出一些基因特征信息
print(head(fData(monocle_cds)))

# 筛选出关键的高变异基因用于后续分析
{
HSMM=monocle_cds# 复制数据集
disp_table <- dispersionTable(HSMM)# 获取每个基因的离散度
# 筛选出平均表达量大于0.1且离散度高于模型拟合值的基因
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
}
#黑点表示的就是筛选出来用于后续分析的差异基因
plot_ordering_genes(HSMM)

# 基于关键基因进行降维分析，使用DDRTree算法，降低维度到2个主成分
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
?reduceDimension

##沿着时间轨迹排序细胞（monocle 2.26.0版本这一步会报错）
# 对细胞进行拟时序排序，以时间轨迹排序细胞
HSMM <- orderCells(HSMM)


# --- 第一步：定义一个简单的函数来寻找起点 State ---
get_root_state <- function(cds, target_group_col, target_value){
  # 统计每个 State 中，目标细胞群（这里是 "Low"）的数量分布
  state_table <- table(pData(cds)$State, pData(cds)[[target_group_col]])
  
  # 打印分布表，方便你自己肉眼核对一下
  print("各 State 中细胞类型的分布情况：")
  print(state_table)
  
  # 找到 target_value (即 "Low") 细胞数量最多的那个 State
  # 注意：这假设 "Low" 细胞主要聚集在某一个树梢上
  target_counts <- state_table[, target_value]
  root_state <- as.numeric(names(which.max(target_counts)))
  
  return(root_state)
}

# --- 第二步：计算并应用起点 ---

# 1. 自动计算 "Low" 对应的 State
# 请确保 "Low" 这个拼写与你数据里的完全一致（大小写敏感）
my_root_state <- get_root_state(HSMM, "geneSet.Type", "Low")

print(paste0("检测到 'Low' 细胞主要集中在 State ", my_root_state, "，将以此为起点重新排序。"))

# 2. 指定 root_state 重新运行 orderCells
HSMM <- orderCells(HSMM, root_state = my_root_state)
###好不容易跑出来，最好保存一下
save(HSMM,file = "HSMM.rda")
                                                                                                                                                                                                                                                                                                                                                                                                                                                       load("HSMM.rda")
table(HSMM$geneSet.Type)
#########下面就是快乐的可视化做图#########
### 分别按照细胞状态、时间、细胞类型、细胞来源上色
pdf(file = "拟时序2-1.pdf", width = 4, height = 3)
plot_cell_trajectory(HSMM, color_by = "geneSet.Type")
dev.off()

plot_cell_trajectory(HSMM, color_by = "group")

pdf(file = "拟时序1.pdf", width = 3, height = 3)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()
HSMM$CytoTRACE2_Relative

plot_cell_trajectory(HSMM, color_by = "CytoTRACE2_Relative")

# 将结果按细胞类型分面展示
plot_cell_trajectory(HSMM, color_by = "cellType") +
  facet_wrap(~cellType, nrow = 2)


# 关注特定基因（CTSL和PDLIM2），并展示这些基因在不同细胞类型中的表达情况
blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("CTSL", "PDLIM2")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "cellType",
                  min_expr = 0.1)

{# 筛选在至少10个细胞中表达的基因
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,] # 保留高表达基因

}
# 提取并关注特定核心基因，并绘制这些基因在拟时序中的表达变化
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CTSL", "PDLIM2")))
cds_subset <- HSMM_filtered[my_genes,]

plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")

plot_genes_in_pseudotime(cds_subset, color_by =  "State")
plot_genes_in_pseudotime(cds_subset, color_by =  "cellType")


# 绘制小提琴图，展示基因CTSL和PDLIM2在不同细胞类型中的表达分布
genes <- c("CTSL", "PDLIM2")
plot_genes_violin(HSMM[genes,], grouping = "cellType", color_by = "cellType")


# 可视化多个基因在拟时序中的表达变化
marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("Abcc3", "Acsbg1", "Acta2", "Actb", "Actg1", "Actn1", "Adam12", 
                                                        "Adam8", "Adamts1", "Adamts2", "Aebp1", "Akap12", "Akr1b8", 
                                                        "Aldh1a2", "Alox5ap", "Ankrd1", "Anxa1", "Anxa2", "Anxa3", 
                                                        "Anxa4", "Anxa5", "Apoe", "Arhgdib", "Arl4c", "Arpc1b", "Arpc2", 
                                                        "Atp1a1", "Atp1b1", "Atp2b1", "Axl", "B2m", "B4galnt1", "Basp1", 
                                                        "Bcam", "Bcl3", "Bgn", "Bst2", "Btg1", "C3", "C5ar1", "C7", 
                                                        "Capg", "Ccdc80", "Ccl12", "Ccl2", "Ccl6", "Ccl7", "Ccl9", 
                                                        "Cd14", "Cd151", "Cd24a", "Cd44", "Cd59a", "Cd63", "Cd68", 
                                                        "Cd74", "Cd81", "Cd9", "Cd93", "Cfh", "Cfp", "Ckb", "Cks2", 
                                                        "Clca1", "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Clu", "Cmtm3", 
                                                        "Cnn2", "Cnn3", "Col15a1", "Col18a1", "Col1a1", "Col1a2", 
                                                        "Col3a1", "Col4a1", "Col4a2", "Col5a1", "Col5a2", "Col6a1", 
                                                        "Col6a2", "Colec12", "Cotl1", "Cox5a", "Crip1", "Crlf1", 
                                                        "Cryab", "Csf1", "Cst3", "Cstb", "Ctsd", "Ctsb", "Ctsl", "Ctsz", 
                                                        "Ctsc", "Ctsk", "Cxcl1", "Cxcl14", "Cxcl16", "Cxcl17", "Cxcl2", 
                                                        "Cyba", "Cygb", "Cystm1", "Dab2", "Dcn", "Ddr1", "Des", "Dpysl3", 
                                                        "Dstn", "Edn1", "Efemp1", "Efemp2", "Egr1", "Egr2", "Elf3", 
                                                        "Emilin1", "Emp1", "Emp2", "Emp3", "Eno1", "Epcam", "Ermp1", 
                                                        "Ezr", "F13a1", "Fbln1", "Fbln2", "Fbn1", "Fcer1g", "Fcgr3", 
                                                        "Fcgrt", "Fermt2", "Fgg", "Fhl2", "Fkbp10", "Fkbp1a", "Flna", 
                                                        "Fn1", "Fos", "Fstl1", "Fth1", "Ftl1", "Fxyd5", "Gale", "Gas6", 
                                                        "Gatm", "Gdf15", "Gfpt2", "Gja1", "Gkn3", "Glipr2", "Gnai2", 
                                                        "Gnas", "Gng11", "Gpnmb", "Gprc5a", "Gpx1", "Gpx3", "Grn", 
                                                        "Gsn", "H2-Aa", "H2-Ab1", "H2-D1", "H2-Eb1", "H2-K1", "H2-Q7", 
                                                        "Hmbox1", "Hoxb9", "Hsp90ab1", "Hspb1", "Hspb8", "Ifi27l2a", 
                                                        "Ifitm2", "Ifitm3", "Igfbp7", "Il18", "Il33", "Inhba", "Islr", 
                                                        "Itga3", "Itga5", "Itga6", "Itgal", "Itgam", "Itgb1", "Itgb2", 
                                                        "Itgb4", "Itgb5", "Jun", "Junb", "Klf4", "Krt18", "Krt19", 
                                                        "Krt20", "Krt7", "Krt8", "Lamp1", "Lamp2", "Lap3", "Lcn2", 
                                                        "Lgals1", "Lgals3", "Lgals3bp", "Lif", "Lilrb4", "Loxl1", 
                                                        "Loxl2", "Lpl", "Lrp1", "Lsp1", "Lyz2", "Mal", "Marcks", 
                                                        "Marcksl1", "Met", "Mfap4", "Mgp", "Mgst1", "Micall2", "Mif", 
                                                        "Mmp14", "Mmp19", "Mmp2", "Mmp23", "Mmp3", "Mmp7", "Mrc2", 
                                                        "Ms4a6c", "Ms4a6d", "Msn", "Muc1", "Mxra7", "Mxra8", "Myl12a", 
                                                        "Myl6", "Myof", "Nbl1", "Ndrg1", "Nfkbiz", "Ngf", "Nid1", "Nid2", 
                                                        "Npc2", "Oit1", "P4hb", "Pdgfrb", "Pdlim4", "Pea15a", "Pf4", 
                                                        "Pfn1", "Pgk1", "Phlda3", "Pkp3", "Plac8", "Plet1", "Plin2", 
                                                        "Pltp", "Ppl", "Prdx1", "Prdx5", "Prdx6", "Prss23", "Ptn", 
                                                        "Ptrf", "Rab25", "Rac1", "Rcn3", "Rgs16", "Rplp1", "S100a10", 
                                                        "S100a11", "S100a13", "S100a4", "S100a6", "Sbno2", "Scn1b", 
                                                        "Scx", "Sdc1", "Sdc4", "Sepn1", "Serpine1", "Serpine2", 
                                                        "Serping1", "Sfn", "Sh3bgrl3", "Slc3a2", "Smoc2", "Socs3", 
                                                        "Sod2", "Sox4", "Sox9", "Sparc", "Spint2", "Spns2", "Spon1", 
                                                        "Spp1", "Sprr1a", "Sprr2f", "Src", "Srxn1", "Stmn1", "S100a16", 
                                                        "Tacstd2", "Tagln", "Tagln2", "Tcn2", "Tenc1", "Tgfb1", "Tgfb2", 
                                                        "Tgfbi", "Thbs1", "Thy1", "Timp1", "Timp2", "Timp3", "Tm4sf1", 
                                                        "Tmem119", "Tmem158", "Tmem176a", "Tmem176b", "Tmsb10", "Tmsb4x", 
                                                        "Tnfaip2", "Tnfaip3", "Tnfrsf12a", "Tnfrsf1b", "Tpm1", "Tpm4", 
                                                        "Tspan8", "Tubb5", "Tyrobp", "Ubb", "Ucp2", "Vcan", "Vim", "Vwf", 
                                                        "Wfdc2", "Wnt4", "Ywhaz")
))
{# 差异基因测试，确定哪些基因在拟时序中的表达变化显著
diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1)) # 提取显著基因名
}
# 绘制拟时序热图，展示基因在不同时间点的表达模式
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 6,  # 聚类成4组
                        cores = 1,  # 使用1个核心计算
                        show_rownames = T)  # 显示基因名


