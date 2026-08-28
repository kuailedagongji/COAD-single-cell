library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel) 
library(magrittr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)

yjsl=readRDS("af3_注释后.rds")

#h获取差异基因，耗时间！
markers <- FindAllMarkers(yjsl, only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 0)


save(markers,file = "markers.rda")
load("markers.rda")

#做图

p1=jjVolcano(diffData = markers)
p1


#标记自己的目标基因

mygene <- c('LTB','CD79B','CCR7')

{
p2=jjVolcano(diffData = markers,
          myMarkers = mygene)
}
p2

#修改点颜色
p3=jjVolcano(diffData = markers,
          aesCol = c('aquamarine3','lightgoldenrod3'))

p3

pdf(file = "差异火山图.pdf", width = 14, height = 4)
##更加高级的
{
p4=markerVolcano(markers = markers,
              topn = 5,
              labelCol = ggsci::pal_npg()(11))
}

p4
dev.off()
