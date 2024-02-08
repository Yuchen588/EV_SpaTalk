setwd("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_SpaTalkdb")
#write.csv(CellChatDB_expanded.filter, file="CellChatDB_expanded.filter.csv")

CellChatDB_expanded.filter <- read.csv(file="CellChatDB_expanded.filter.csv", header = T, row.names = 1, as.is = T)

database <- database.new <- cbind.data.frame(Pair.Name=paste(CellChatDB_expanded.filter$ligand.new.id, CellChatDB_expanded.filter$receptor.new.id, sep = "_"),
                                             Ligand.ApprovedSymbol=CellChatDB_expanded.filter$ligand.new.id,
                                             Ligand.Name=CellChatDB_expanded.filter$ligand.secreted_type,
                                             Receptor.ApprovedSymbol=CellChatDB_expanded.filter$receptor.new.id,
                                             Receptor.Name=CellChatDB_expanded.filter$receptor.secreted_type,
                                             Classification=CellChatDB_expanded.filter$annotation)
database.new <- database.new[database$Classification != "Non-protein Signaling",]

load("/home/lyc/sc/3CA/rawdata/meta.program/tumor/results/modules.RData")


save(db_mac, file="/home/lyc/sc/3CA/rawdata/meta.program/tumor/results/db_mac.Rdata")
load("/home/lyc/sc/3CA/rawdata/meta.program/tumor/results/db_mac.Rdata")
#required package
#library(iTALK)
library(Seurat)
#library(iTALK)
library(tidyverse)
library(future)
library(nnls)
library(parallel)
library(Matrix)
library(parallel)
library(Matrix)
library(VennDiagram)
library(ggplot2)
library(rstatix)

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(RobustRankAggreg)
library(gridExtra)
library(circlize)




#创建EV_spatalk_object
setClass("EV_spatalk",
         slots = c(database = "data.frame",
                   all.cell.type = "character",
                   s.cell.type = "character",
                   Sender.celltype.id = "character",
                   Receiver.celltype.id = "character",
                   Sender.spot.id = "character",
                   Receiver.spot.id = "character",
                   Distance.corr.LR.results = "data.frame",
                   EVrelease.corr.LR.results = "data.frame",
                   inter.LR.results = "list",
                   all.can.LR.id="character",
                   st.seurat.obj = "Seurat",
                   sc.seurat.obj = "Seurat",
                   all.spot.module.score = "data.frame",
                   distance.results = "data.frame",
                   neighborhood.results = "list",
                   metadata = "data.frame",
                   interaction_df = "data.frame",
                   interaction_df_raw.indensity = "data.frame",
                   modules="array",
                   db_mac="list",
                   db_tcell="list",
                   nnls_bin="character",
                   nearest.spot.list = "list",
                   all.nitch.LR.talk.list="list",
                   all.spot.lr_interaction="list",
                   sender.spot.lr_interaction="list",
                   spot.LR.freq.results="data.frame",
                   LR.in.spot.frequency_table="data.frame",
                   LR.in.spot.RRA.results="data.frame",
                   EV_spatalk_stat_results="data.frame"))


EV_spatalk_object <- new("EV_spatalk")
EV_spatalk_object@database <- database.new
EV_spatalk_object@modules <- modules
EV_spatalk_object@db_mac <- db_mac
EV_spatalk_object@st.seurat.obj <- st
EV_spatalk_object@sc.seurat.obj <- sc

saveRDS(EV_spatalk_object, file="/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/source.data/EV_spatalk_object.rds")



#load the EV_spatalk_object
####
###
##
#demo
output.file.path <- "/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/results"

EV_spatalk_object <- readRDS(file = "/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/source.data/EV_spatalk_object.rds")

st <- readRDS("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/rawdata/st.rds")
sc <- readRDS("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/rawdata/sc.rds")

source("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/Rcode/functions_public.R")
source("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/Rcode/functions_main.R")

cell.types <- unique(sc$pop)
select.celltype <- cell.types
print(cell.types)
#other parameters need to be setted
nrand = 2 # interation times, 2 option is indicate 10^2
nbin = 5
pval_thresh=0.05
mc.cores = ncores = 10
method = 'rand'
set.seeds = 1

custom_colors <- c("Malignant" = "#4B96B5",
                   "Epi" = "#9CED5F",
                   "T_cells" = "#5201BC",
                   "Macrophage" = "#DC153A",
                   "B_cell" = "#C785C8",
                   "Fibro" = "#EFE633",
                   "Endothelial_cells"="#FFB27D",
                   "NK_cell"="#AAAEEB",
                   "DC"="#8F5861", "Neutrophils"="#BF5B16")


EV.spatalk.results <- st_deco_anno(st = st, sc = sc, nrand = 2, EV_spatalk_object=EV_spatalk_object, nbin = 5, pval_thresh=0.05, mc.cores=30, sc.anno.id="pop", set.seeds = 1, select.celltype = select.celltype)
saveRDS(EV.spatalk.results, file="/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/results/EV.spatalk.results.rds")

#EV mediated spatiallt ccc
EV.spatalk.results <- readRDS("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/results/EV.spatalk.results.rds")

s.cell.type <- c("Malignant", "T_cells")#这里可以自定义细胞类型
mc.cores=30
comm_list= c("Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling")
datatype='mean count'
#find_niche_LR是求每个spot中的LR互作情况

EV.spatalk.results <- find_niche_LR(EV_spatalk_object=EV.spatalk.results, prox="inverse", mc.cores=30, s.cell.type = c("Malignant", "T_cells"), comm_list=comm_list, datatype='mean count', method = "pseudocount")
saveRDS(EV.spatalk.results, file="/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/results/EV.spatalk.results.rds")

#find.inter.LR寻找每个niche中overlap的LR
EV.spatalk.results <- find.inter.LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30)

#find_EV_spatalk_LR寻找有统计学意义的LR并和EV分泌以及空间距离相关联
EV.spatalk.results <- find_EV_spatalk_LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30, seeds = 1)


EV.spatalk.results@EV_spatalk_stat_results
EV_spatalk_stat_results <- EV.spatalk.results@EV_spatalk_stat_results

#可视化
#所有候选LR的robust情况，即EV.spatalk.results中701个LR的
##
# 筛选出p值小于0.05



################
#########
#####
#空间ST展示LR
EV.spatalk.results <- add_interaction_score(EV_spatalk_object = EV.spatalk.results, mc.cores=30)
saveRDS(EV.spatalk.results, file="/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/results/EV.spatalk.results.rds")

#到此位置分析结束
#Idents(st) <- "SCT"
#SpatialFeaturePlot(st, features = "EV_release", image.alpha = 0.6)

# st now contains the interaction data in its metadata
st.sender <- subset(st, cells=EV.spatalk.results@Sender.spot.id)
past.id <- paste0("Sender:", s.cell.type[1], " | ", "Receiver:", s.cell.type[2])
past.id



#LR indensity空间分布
LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="CD86_CTLA4")#这是EVspatalk的作图function

s.LR.pair <- c("LGALS9_HAVCR2", "PVR_TIGIT", "CD274_PDCD1", "HLA-DRB1_CD4")
p1 <- LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="LGALS9_HAVCR2")#这是EVspatalk的作图function
p2 <- LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="PVR_TIGIT")#这是EVspatalk的作图function
p3 <- LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="CD274_PDCD1")#这是EVspatalk的作图function
p4 <- LR_spatial_indensity_plot(EV_spatalk_object=EV.spatalk.results, s.LR.pair="HLA-DRB1_CD4")#这是EVspatalk的作图function
p1+p2+p3+p4

###
##
#
#饼图展示
LR_pie.plot(EV_spatalk_object=EV.spatalk.results)

#韦恩图
LR_venn.plot(EV_spatalk_object=EV.spatalk.results)
dev.off()

#circlize图
set.seed(2)
#select.LR.id <- sample(LR.freq.res$LR, 90)
select.LR.id <- c(EV.spatalk.results@inter.LR.results$common_positive_LR, EV.spatalk.results@inter.LR.results$common_negative_LR)
LR_circos.plot(EV_spatalk_object = EV.spatalk.results, select.LR.id=select.LR.id, specific.LR.id="CD96_PVR")
circos.clear()



#创建R包
library(devtools)
create_package("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/Rpackage/EVSpaTalk")



#
devtools::build(path = "/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/Rpackage/EVSpaTalk")


#完成document撰写
setwd("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/Rpackage/EVSpaTalk")
devtools::document()


#测试成功
devtools::install_github("Yuchen588/EV_SpaTalk@master", force = T)
library(EVSpaTalk)
remove.packages("EVSpaTalk")

#读取R包
# 加载所需的包
library(Seurat)
library(tidyverse)
library(future)
library(nnls)
library(parallel)
library(Matrix)
library(VennDiagram)
library(ggplot2)
library(rstatix)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(RobustRankAggreg)
library(gridExtra)
library(circlize)
#创建create_EV.SpaTalk.obj对象
if ("EV_spatalk" %in% getClassDef()) {
  removeClass("EV_spatalk")
}
create_EV.SpaTalk.obj()

EV_spatalk_object <- readRDS(file = "/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/source.data/EV_spatalk_object.rds")
st <- readRDS("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/rawdata/st.rds")
sc <- readRDS("/home/lyc/sc/3CA/rawdata/meta.program/tumor/NG_pancaner_ST/Analysis/EV-spatalk/EV_spatalk_demo/rawdata/sc.rds")

cell.types <- unique(sc$pop)
select.celltype <- cell.types
print(cell.types)

#main process
EV.spatalk.results <- st_deco_anno(st = st, sc = sc, nrand = 2, EV_spatalk_object=EV_spatalk_object, nbin = 5, pval_thresh=0.05, mc.cores=30, sc.anno.id="pop", set.seeds = 1, select.celltype = select.celltype)
#
s.cell.type <- c("Malignant", "T_cells")
mc.cores=30
comm_list= c("Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling")
datatype='mean count'
EV.spatalk.results <- find_niche_LR(EV_spatalk_object=EV.spatalk.results, prox="inverse", mc.cores=30, s.cell.type = c("Malignant", "T_cells"), comm_list=comm_list, datatype='mean count', method = "pseudocount")
#
EV.spatalk.results <- find.inter.LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30)
#
EV.spatalk.results <- find_EV_spatalk_LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30, seeds = 1)
#
EV.spatalk.results <- add_interaction_score(EV_spatalk_object = EV.spatalk.results, mc.cores=30)
