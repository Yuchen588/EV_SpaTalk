
#' Decompose cell types in spatial transcriptomics data
#'
#' Identify the cellular composition for single-cell or spot-based spatial transcriptomics data with non-negative regression.
#'
#' @param EV_spatalk_object An S4 object reconstructed for EV_SpaTalk. This object must correctly include both spatial transcriptomics (st) and single-cell RNA-seq (sc) data entries.
#' @param sc.anno.id A character containing the cell type of the reference single-cell RNA-seq data.
#' @param select.celltype A character the all selected cell types for ST spot deconvolution.
#' @return An updated EV_spatalk_object containing the results of the decomposition analysis, including cell type distributions and spatial distances between identified cell types.
#' @export
#' @examples
#' EV.spatalk.results <- st_deco_anno(st = st, sc = sc, nrand = 2, EV_spatalk_object=EV_spatalk_object, nbin = 5, pval_thresh=0.05, mc.cores=30, sc.anno.id="pop", set.seeds = 1, select.celltype = select.celltype)
st_deco_anno <- function(st = st.data, sc = sc.data, EV_spatalk_object=EV_spatalk_object, nrand = 2, nbin = 5, pval_thresh=0.05, mc.cores=30, sc.anno.id="pop", set.seeds = 1, select.celltype = select.celltype, pred_bin_cutoff=0.9, coef_bin_cutoff=2){
  ncores = mc.cores
  EV_spatalk_object@all.cell.type <- select.celltype
  modules <- EV_spatalk_object@modules#
  ma <- names(modules)
  print(ma)
  st@assays$SCT <- NULL
  st@assays$Spatial <- NULL
  st@assays$integrated <- NULL


  sc@assays$SCT <- NULL
  sc@assays$Spatial <- NULL
  sc@assays$integrated <- NULL

  DefaultAssay(st) <- "RNA"
  DefaultAssay(sc) <- "RNA"

  print("Processing SCTransform on sc and st data")

  st = SCTransform(st, return.only.var.genes = FALSE, vst.flavor = "v1")
  sc = SCTransform(sc, return.only.var.genes = FALSE, vst.flavor = "v1")

  DefaultAssay(st) <- "SCT"
  DefaultAssay(sc) <- "SCT"

  set.seed(set.seeds)

  cell.types <- unique(sc@meta.data[,sc.anno.id])

  columns_to_remove <- c('subclone', 'sumsq', 'topcor', 'Cycle', 'Stress', 'Interferon',
                         'Hypoxia', 'pEMT', 'Squamous', 'Squamous2', 'Glandular', 'Glandular2',
                         'Cilium', 'Metal', 'Undetermined', 'Undetermined2', 'Glandular1',
                         'Squamous1', 'AC', 'OPC', 'NPC', 'Mesenchymal', 'Ductal', 'Luminal',
                         'Keratinocyte', 'Basal', names(modules))

  # 删除这些列
  ST.data <- st
  ST.data@meta.data <- ST.data@meta.data[,!colnames(ST.data@meta.data) %in% columns_to_remove]


  ST.data.mscore <- ST.data
  #这里计算module score
  #options(Seurat.object.assay.version = 'v4')

  set.seed(set.seeds)
  modules_rand = MakeRand(ST.data.mscore, db = modules, nrand = nrand, nbin = nbin)
  ini = matrix(0,nrow = ncol(ST.data.mscore), ncol = length(modules))
  rownames(ini) = colnames(ST.data.mscore)
  colnames(ini) = names(modules)
  ST.data.mscore@meta.data[,names(modules)] = as.data.frame(ini)

  for (m in names(modules)){
    tryCatch(expr = {
      ST.data.mscore = GeneToEnrichment(ST.data.mscore, db = modules[m], method = 'rand', db_rand = modules_rand[m])
    }, error = function(e){c()})
  }

  scores = ST.data.mscore@meta.data[,names(modules)]
  scores.1 <- scores
  EV_spatalk_object@all.spot.module.score <- scores.1#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值每个spot的module score
  scores[is.na(scores)] <- 0#NA部分替换为0

  frequency = colMeans(scores > 0.5, na.rm = TRUE)
  ST.data.mscore$state = apply(ST.data.mscore@meta.data[,names(modules)], 1, function(x){
    top = which.max(x)
    return(names(x)[top])
  })
  ST.data.mscore$state = factor(ST.data.mscore$state, levels = names(modules))



  ST.data.mscore@meta.data[is.na(ST.data.mscore@meta.data)] <- 0#NA部分替换为0

  ma_bin = paste0(ma, '_bin')
  scores = ST.data.mscore@meta.data[,ma]
  scores_bin = scores
  scores_bin[] = as.numeric(scores_bin > 0.5)
  ST.data.mscore = AddMetaData(ST.data.mscore, metadata = scores_bin, col.name = ma_bin)


  nmf <- names(modules)

  st <- ST.data.mscore
  st = AddMetaData(st, metadata = st@images$image@coordinates[,c('row','col')], col.name = c('row','col'))
  st$axis = st$row
  st$axis = (st$axis - min(st$axis))/(max(st$axis) - min(st$axis))

  #
  print("Processing intergration of sc and st data")

  srt <- sc

  options(future.globals.maxSize = 5000*1024^2)

  srt@assays$integrated <- NULL
  srt@assays$Spatial <- NULL

  st@assays$integrated <- NULL
  st@assays$Spatial <- NULL

  common_genes <- intersect(rownames(srt), rownames(st))

  Idents(srt) <- "RNA"
  Idents(st) <- "RNA"
  DefaultAssay(srt) <- "RNA"
  DefaultAssay(st) <- "RNA"

  srt@assays$SCT <- NULL
  st@assays$SCT <- NULL

  srt@images$image <- NULL
  srt@images$image.1 <- NULL

  st@images$image <- NULL
  st@images$image.1 <- NULL

  srt <- srt[common_genes,]
  st <- st[common_genes,]

  #st@images$image <- ST.data.mscore@images$image #重新导入image
  srt = SCTransform(srt, return.only.var.genes = FALSE, vst.flavor = "v1") %>% RunPCA() #%>% RunUMAP(dims = 1:20)
  st = SCTransform(st, return.only.var.genes = FALSE, vst.flavor = "v1") %>% RunPCA() #%>% RunUMAP(dims = 1:20)

  object.list = list('SC' = srt, 'ST' = st)


  #
  print("Take the intersection of variable genes between st and sc data")
  genes.use = intersect(VariableFeatures(st), VariableFeatures(srt))

  print("Processing SCT integration")
  object.list = PrepSCTIntegration(object.list = object.list,
                                   anchor.features = genes.use,
                                   verbose = FALSE)


  anchors = FindTransferAnchors(reference = object.list$SC, query = object.list$ST,
                                normalization.method = 'SCT',
                                features = genes.use,
                                verbose = FALSE, reduction = "rpca")



  predictions = TransferData(anchorset = anchors, refdata = srt@meta.data[,sc.anno.id])
  predictions = predictions[,!colnames(predictions) %in% c('predicted.id','prediction.score.max')]
  colnames(predictions) = gsub('prediction.score','pred',colnames(predictions))
  pred = colnames(predictions)
  st = AddMetaData(st, metadata = predictions, col.name = pred)#这步增加了基于单细胞的st spot预测概率
  print("The integration of sc and st finished")



  print("Binarizing predictions")
  predictions_bin = predictions
  predictions_bin[] = as.numeric(predictions_bin > pred_bin_cutoff)#这步cutoff是0.9导致很多细胞被认为是肿瘤细胞, 但是没问题后面有nnls
  pred_bin = paste0(colnames(predictions_bin), '_bin')
  st = AddMetaData(st, predictions_bin, col.name = pred_bin)

  print("Deconvoluting st spot from paired scRNA-seq data")
  genes.use = intersect(VariableFeatures(st), VariableFeatures(srt))

  Idents(srt) <- sc.anno.id
  prof = AverageExpression(srt, assay = 'SCT', layer = 'data')$SCT[genes.use,]
  data = as.matrix(GetAssayData(st, assay = 'SCT', layer = 'data'))[genes.use,]


  print("nnls regression")

  coef = t(apply(data, 2, function(y){
    coef(nnls(as.matrix(prof), y))
  }))
  colnames(coef) = colnames(prof)
  nnls = colnames(coef)[colSums(coef > 0) >= 0]
  prof = prof[,nnls]
  coef = t(apply(data, 2, function(y){
    coef(nnls(as.matrix(prof), y))
  }))

  nnls <- gsub("-", "_", nnls, fixed = T)

  colnames(coef) = nnls
  st = AddMetaData(st, coef, col.name = nnls)

  print("Processing nnls coefficient scaling")
  #
  colnames(coef) <- gsub("-", "_", colnames(coef), fixed = T)#
  #x="T_cells"
  coef_scaled = sapply(nnls, function(x){
    vec = coef[,x]
    y = paste0('pred.',x)
    spots.use = colnames(st)[st@meta.data[,y] == 0]
    #spots.use = colnames(st)[st@meta.data[,y] == 0 & rowSums(predictions_bin) == 1]
    if (length(spots.use) < 5){
      spots.use = colnames(st)[order(st@meta.data[,y])[1:5]]
    }

    nrand=1#原来是2， 100次迭代
    if(nrand==2){nrand=nrand-1}

    data_rand = Reduce(cbind, lapply(1:10^nrand, function(i){
      t(apply(data[,spots.use], 1, function(expr){
        sample(expr, length(expr), replace = FALSE)
      }))
    }))
    coef_rand = t(apply(data_rand, 2, function(y){
      coef(nnls(as.matrix(prof), y))
    }))
    colnames(coef_rand) = nnls
    if (sd(coef_rand[,x]) == 0){
      return((coef[,x] - mean(coef_rand[,x]))/min(coef[coef[,x] > 0,x]))
    } else {
      return((coef[,x] - mean(coef_rand[,x]))/sd(coef_rand[,x]))
    }
  })#>到这是对nnls的coef进行标准化
  nnls_scaled = paste0(nnls, '_scaled')#nnls_scale是coef_scaled 进行nnls的coefficient进行标准化
  st = AddMetaData(st, coef_scaled, col.name = nnls_scaled)

  print("nnls binarizing")
  coef_bin = coef_scaled
  coef_bin[] = (coef_bin > coef_bin_cutoff)#scale的nnls系数相加大于2可以认为是不同类型的细胞
  nnls_bin = paste0(nnls, '_bin')

  EV_spatalk_object@nnls_bin <- nnls_bin#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值nnls_bin
  st = AddMetaData(st, coef_bin, col.name = nnls_bin)
  nnls_scaled = paste0(nnls, '_scaled')
  nnls_bin = paste0(nnls, '_bin')

  #这块就先不加macro和T细胞的
  st@images<- ST.data@images
  print("add neighboor")
  DefaultAssay(st) <- "SCT"
  nei = FindSTNeighbors(st, d_min = 0, d_max = 1.5)#FindSTNeighbors要从seurat_fucntions_public.R里面找
  EV_spatalk_object@neighborhood.results <- nei#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值neighboorhood

  #先定义每个spot周围的spotID, 在metadata中找到这些spot，看每种细胞类型在的占比即mean值。
  coef_nei = sapply(nnls_bin, function(x){

    #spot = nei[1]

    sapply(nei, function(spot){
      y = st@meta.data[spot,x]
      y[y < 0] = 0
      mean(y, na.rm = TRUE)
    })
  })
  colnames(coef_nei) = nnls
  nnls_nei = paste0(nnls, '_nei')
  st = AddMetaData(st, coef_nei, col.name = nnls_nei)

  print("Add cell distance information from ST")
  coord = st@images$image@coordinates[,c('imagerow','imagecol')]

  prox = 'inverse'
  distances = as.matrix(dist(coord))
  distances = distances/min(distances[distances > 0])#矩阵标准化
  distance.all.spot <- 1/(1+distances)
  EV_spatalk_object@distance.results <- as.data.frame(distance.all.spot)#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值distance



  #
  coef_dist = sapply(nnls_bin, function(x){
    w = colnames(st)[as.logical(st@meta.data[,x])]
    w = w[!is.na(w)]
    if (length(w) == 1){
      mi = as.numeric(distances[,w])
    } else {
      mi = apply(distances[,w], 1, function(x){min(as.numeric(x), na.rm = TRUE)})
    }
    return(mi)
  })
  if (prox == 'inverse'){
    coef_dist = 1/(1+coef_dist)
  }
  if (prox == 'opposite'){
    coef_dist = -coef_dist
  }
  colnames(coef_dist) = nnls
  nnls_dist = paste0(nnls, '_dist')#这边同样nnls_dist包括所有
  st = AddMetaData(st, coef_dist, col.name = nnls_dist)


  print("Cell categorization|Classify all cells into normal, both and malignant")

  # Categories
  st$cat = apply(st@meta.data[,nnls_bin], 1, function(x){
    if (x['Malignant_bin']){
      if (sum(x) == 1){
        return('Malignant')
      } else {
        return('Both')
      }
    } else {
      if (sum(x) == 0){
        return(NA)
      } else {
        return('Normal')
      }
    }
  })
  cats = c('Malignant','Both','Normal')
  st$cat = factor(st$cat, levels = cats)

  #add new seurat slot
  #slot(EV_spatalk_object, "st.seurat.obj") <- st  # 替换为您的实际ST Seurat对象
  #slot(EV_spatalk_object, "sc.seurat.obj") <- sc  # 替换为您的实际SC Seurat对象

  #add your gene signature of cell states

  EV_spatalk_object@st.seurat.obj <- st#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值最新的st
  EV_spatalk_object@sc.seurat.obj <- sc#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值最新的sc

  EV.spatalk.results <- EV_spatalk_object
  return(EV.spatalk.results)
}


#' Screening the sender-receiver niche and ligand-receptor interactome
#'
#' This function identifies the spatially nearest sender-receiver niche for EV-mediated cell-cell interactions. It composes EV-related ligand-receptor (LR) pairs and analyzes their patterns based on spatial proximity.
#'
#' @param EV_spatalk_object An S4 object reconstructed for EV_SpaTalk, which should include both spatial transcriptomics data and the reference single-cell RNA-seq data.
#' @param prox The pattern to identify the spatial distance. "inverse" indicates that the closer the distance, the greater the value is considered in the interaction strength.
#' @param s.cell.type A character vector specifying the sender (first cell type) and the receiver (second cell type) in the interaction.
#' @param comm_list The list of ligand-receptor groups used for identifying EV-mediated interactomes.
#' @param datatype The type of data to use for calculating the intensity of LR interactions, defaulting to the 'data' slot from SCTransform normalization.
#' @param method The method used to determine the strength of each LR interaction. The default "pseudocount" method involves a log transformation after summing the values of ligand and receptor.
#' @return An updated EV_spatalk_ object with niche identifier and all LR interaction results.
#' @export
#' @examples
#' EV.spatalk.results <- find_niche_LR(EV_spatalk_object=EV.spatalk.results, prox="inverse", mc.cores=30, s.cell.type = c("Malignant", "T_cells"), comm_list=comm_list, datatype='mean count', method = "pseudocount")
find_niche_LR <- function(EV_spatalk_object=EV.spatalk.results, prox="inverse", mc.cores=30, s.cell.type = c("Malignant", "T_cells"), comm_list=comm_list, datatype='mean count', method="pseudocount"){
  numCores = mc.cores
  st <- EV_spatalk_object@st.seurat.obj
  EV_spatalk_object@s.cell.type <- s.cell.type#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值s.cell.type

  nnls_bin <- EV_spatalk_object@nnls_bin

  coord = EV_spatalk_object@st.seurat.obj@images$image@coordinates[,c('imagerow','imagecol')]
  distances = as.matrix(dist(coord))
  distances = distances/min(distances[distances > 0])#矩阵标准化

  coef_dist = sapply(nnls_bin, function(x){
    w = colnames(st)[as.logical(EV_spatalk_object@st.seurat.obj@meta.data[,x])]
    w = w[!is.na(w)]
    if (length(w) == 1){
      mi = as.numeric(distances[,w])
      mi.ID = w

    } else {
      mi = apply(distances[,w], 1, function(x){min(as.numeric(x), na.rm = TRUE)})
      mi.ID = apply(distances[,w], 1, function(x){names(which.min(x))})
    }
    return(mi)
    #return(mi.ID)
  })


  if (prox == 'inverse'){
    coef_dist = 1/(1+coef_dist)
  }
  colnames(coef_dist) = nnls_bin


  print("finding nearest spot id")

  coef_id = sapply(nnls_bin, function(x){
    w = colnames(st)[as.logical(EV_spatalk_object@st.seurat.obj@meta.data[,x])]
    w = w[!is.na(w)]
    if (length(w) == 1){
      mi = as.numeric(distances[,w])
      mi.ID = w

    } else {
      mi = apply(distances[,w], 1, function(x){min(as.numeric(x), na.rm = TRUE)})
      mi.ID = apply(distances[,w], 1, function(x){names(which.min(x))})
    }
    #return(mi)
    return(mi.ID)
  })


  colnames(coef_dist) <- gsub("_bin", "", colnames(coef_dist), fixed = T)
  colnames(coef_id) <- gsub("_bin", "", colnames(coef_id), fixed = T)

  list.results <- list()
  list.results[["cell_dist"]] <- as.data.frame(coef_dist)
  list.results[["cell_id"]] <- as.data.frame(coef_id)
  #list.results#是包含最近距离和最近细胞的两个矩阵，可以导出了
  EV_spatalk_object@nearest.spot.list <- list.results #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值最近距离和最近细胞

  #
  #s.cell.type <- c("Malignant", "T_cells")#这里可以自定义细胞类型
  sender.receiver.bi.res <- EV_spatalk_object@st.seurat.obj@meta.data[,c(paste(s.cell.type, "_bin", sep = ""))]
  sender.cell.id <- rownames(sender.receiver.bi.res[sender.receiver.bi.res[,1]==1,])#sender的spotID
  receiver.cell.id <- rownames(sender.receiver.bi.res[sender.receiver.bi.res[,2]==1,])#receiver的spotID

  EV_spatalk_object@Sender.spot.id <- sender.cell.id#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值所有sender spotid
  EV_spatalk_object@Receiver.spot.id <- receiver.cell.id#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值所有Receiver.spot.id

  #s.cell.type <- s.cell.type
  print(paste0("Sender:", s.cell.type[1], " | ", "Receiver:", s.cell.type[2]))

  nitch.id <- rownames(list.results$cell_id)#所有spot的id


  all.nitch.LR.talk.list <- list()

  # 使用mclapply进行并行计算

  #s.cell.type <- nnls_bin
  all.nitch.LR.talk.list <- mclapply(nitch.id, function(i) {
    nitch.anno <- cbind.data.frame(cell.id=as.character(list.results$cell_id[i,][,s.cell.type]),
                                   anno=names(list.results$cell_id[i,][,s.cell.type]))
    nitch.anno$cb.cell.id <- paste(nitch.anno$cell.id, nitch.anno$anno, sep = "|")

    iTalk_data <- as.data.frame(t(st@assays$SCT$data[,nitch.anno$cell.id]))
    rownames(iTalk_data) <- nitch.anno$cb.cell.id
    iTalk_data$cell_type <- nitch.anno$cb.cell.id

    unique(iTalk_data$cell_type)
    #unique(iTalk_data$compare_group)
    highly_exprs_genes <- rawParse2(iTalk_data, top_genes=200, stats="mean")#这里可以改成10000
    # 通讯类型
    #comm_list<-c('growth factor','other','cytokine','checkpoint')#这是所有的可以自己设置
    #comm_list<-c('checkpoint')

    cell_types <- unique(iTalk_data$cell_type)
    #cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)

    #comm_type = comm_list[1]
    #database <- database[database$Classification == comm_type,]#这边需要换成cellchat的database


    #comm_list= c("Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling")
    #comm_list= comm_list[c(1, 2, 4)]
    iTalk_res <- NULL
    for(comm_type in comm_list){
      res_cat <- FindLR2(highly_exprs_genes, datatype='mean count', comm_type=comm_type)#这里不对随时转换为FindLR2
      iTalk_res <- rbind(iTalk_res, res_cat)
    }


    if(method=="weighted.sum"){
      iTalk_res$multiply <- log2(iTalk_res$cell_from_mean_exprs+1) + log2(iTalk_res$cell_to_mean_exprs+1)#这步骤是weighted sum
    }


    #iTalk_res$multiply <- log2(iTalk_res$cell_from_mean_exprs+1) + log2(iTalk_res$cell_to_mean_exprs+1)#这步骤是weighted sum

    #pesudocount
    # 应用函数计算互作强度
    if(method=="pseudocount"){
      iTalk_res$multiply <- calculate_interaction_strength(iTalk_res$cell_from_mean_exprs, iTalk_res$cell_to_mean_exprs)#pseudocount method
    }



    iTalk_res <- iTalk_res[order(iTalk_res$multiply, decreasing = T),]
    # 返回iTalk_res
    return(iTalk_res)
  }, mc.cores = numCores)
  names(all.nitch.LR.talk.list) <- nitch.id
  #
  EV_spatalk_object@all.nitch.LR.talk.list <- all.nitch.LR.talk.list#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>在确定好sender和receiver之后，赋值所有spot的受配体表达情况
  return(EV_spatalk_object)
}

#' Identify overlapping LR interactomes across niches
#'
#' Analyzes the frequency of ligand-receptor (LR) pairs enriched in spatially adjacent spots to identify those that are consistently present across different niches. This function facilitates the parallel analysis of LR interactomes and their interaction intensities.
#'
#' @param EV_spatalk_object An updated S4 object that has been prepared for EV_SpaTalk analysis.
#' @param mc.cores The number of cores to use for parallel computation. Adjust according to your server's capabilities; use 1 for a single-core PC.
#' @return Returns the EV_spatalk_object updated with identified overlapping LR pairs and their interaction intensities across all examined niches.
#' @export
#' @examples
#' EV.spatalk.results <- find.inter.LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30)
find.inter.LR <- function(EV_spatalk_object=EV.spatalk.results, mc.cores=30){
  s.cell.type <- EV_spatalk_object@s.cell.type
  sender.receiver.bi.res <- EV_spatalk_object@st.seurat.obj@meta.data[,c(paste(s.cell.type, "_bin", sep = ""))]
  Sender.spot.id <- EV_spatalk_object@Sender.spot.id#sender的spotID
  Receiver.spot.id <- EV_spatalk_object@Receiver.spot.id#receiver的spotID

  all.nitch.LR.talk.list <- EV_spatalk_object@all.nitch.LR.talk.list

  lr_interaction_list <- list()
  # （1）遍历所有spot的LR交互情况并提取multiply值; #所有的spot从里面只挑选Malignant_bin|T_cells_bin的受配体
  direction = paste(s.cell.type[1], s.cell.type[2], sep = "|")
  print(paste("The direction of your EV-mediated cell-cell interation pattern is", direction, sep = ": "))
  #c("Malignant_bin|T_cells_bin")#方要设定好
  #lr_interaction_list 这步lr_interaction_list 可以选择不同的spot
  lr_interaction_list <- mclapply(names(all.nitch.LR.talk.list), function(spot_id) {
    lr_data <- all.nitch.LR.talk.list[[spot_id]]
    lr_data$interaction <- paste(lr_data$ligand, lr_data$receptor, sep = "_")

    # 提取cell_from和cell_to中“|”之后的字符
    extract_after_pipe <- function(string) {
      sub(".*\\|", "", string)
    }

    # 合并提取的字符作为方向
    lr_data$direction <- paste(
      sapply(lr_data$cell_from, extract_after_pipe),
      sapply(lr_data$cell_to, extract_after_pipe),
      sep = "|"
    )
    lr_data <- lr_data[lr_data$direction==direction,]
    return(lr_data[, c("interaction", "direction", "multiply")])
  }, mc.cores = mc.cores)

  #所有的spot从里面只挑选Malignant_bin|T_cells_bin的受配体

  names(lr_interaction_list) <- names(EV_spatalk_object@all.nitch.LR.talk.list)

  EV_spatalk_object@all.spot.lr_interaction <- lr_interaction_list #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值所有spot 每个sender 的LR 的sender>receiver情况
  #pick the malignnat as send cell
  EV_spatalk_object@sender.spot.lr_interaction <- lr_interaction_list[Sender.spot.id]#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>赋值sender spot 每个sender 的LR 的sender>receiver情况
  return(EV_spatalk_object)
}

#' Statistical analysis of LR pairs for EV secretion and spatial distance.
#'
#' This function employs the Robust Rank Aggregation (RRA) algorithm to assess the enrichment significance of ligand-receptor (LR) pairs. It also identifies LR pairs that significantly correlate with the EV release capacity of the sender cell and the spatial proximity between sender and receiver cell niches.
#'
#' @param EV_spatalk_object An updated S4 object that has been prepared for EV_SpaTalk analysis.
#' @param mc.cores The number of cores to utilize for parallel processing. This should be adjusted according to the capabilities of your computational environment; set to 1 for single-core machines.
#' @return An updated EV_spatalk_object including the significance analysis of all LR pairs, with additional annotations for their relevance to EV release and spatial distance correlation.
#' @export
#' @examples
#' EV.spatalk.results <- find_EV_spatalk_LR(EV_spatalk_object=EV.spatalk.results, mc.cores=30, seeds = 1)
find_EV_spatalk_LR <- function(EV_spatalk_object=EV.spatalk.results, mc.cores=30, seeds=1){
  #lr_interaction_list <- EV.spatalk.results@sender.spot.lr_interaction
  #sender.cell.id <- EV.spatalk.results@Sender.spot.id
  lr_interaction_list <- EV_spatalk_object@sender.spot.lr_interaction
  all_lr <- unique(unlist(lapply(lr_interaction_list, function(df) df$interaction)))

  list.results <- EV_spatalk_object@nearest.spot.list
  sender.cell.id <- EV_spatalk_object@Sender.spot.id

  # 初始化结果列表
  correlation_results <- list()
  # 遍历所有LR对
  for (lr in all_lr) {
    # 提取所有spot中特定LR的multiply值
    multiply_values <- get_lr_multiply(lr, lr_interaction_list)

    # 提取所有spot的T_cells_bin距离
    distances <- list.results$cell_dist[sender.cell.id, s.cell.type[2]]

    # 计算multiply值与距离的相关性
    cor_test <- cor.test(multiply_values, distances, method = "pearson")

    # 保存相关性结果
    correlation_results[[lr]] <- tidy(cor_test)
  }

  # 将结果转换为数据框
  distance.correlation_results_df <- bind_rows(correlation_results, .id = "LR")

  EV_spatalk_object@Distance.corr.LR.results <- as.data.frame(distance.correlation_results_df)

  #后面再加上个和EV_release的相关性，和差异分析的结果。
  correlation_results <- list()
  # 遍历所有LR对
  for (lr in all_lr) {
    # 提取所有spot中特定LR的multiply值
    multiply_values <- get_lr_multiply(lr, lr_interaction_list)

    # 提取所有spot的T_cells_bin距离
    EV.release.score <- EV_spatalk_object@st.seurat.obj@meta.data[sender.cell.id,]$EV_release

    # 计算multiply值与距离的相关性
    cor_test <- cor.test(multiply_values, EV.release.score, method = "pearson")

    # 保存相关性结果
    correlation_results[[lr]] <- tidy(cor_test)
  }
  # 将结果转换为数据框
  EVrelease.correlation_results_df <- bind_rows(correlation_results, .id = "LR")
  EV_spatalk_object@EVrelease.corr.LR.results <- as.data.frame(EVrelease.correlation_results_df)

  # 筛选出p值小于0.05
  significant_distance <- subset(distance.correlation_results_df, p.value < 0.05)
  significant_EVrelease <- subset(EVrelease.correlation_results_df, p.value < 0.05)

  # 确定正负相关性
  positive_distance <- subset(significant_distance, estimate > 0)
  negative_distance <- subset(significant_distance, estimate < 0)

  positive_EVrelease <- subset(significant_EVrelease, estimate > 0)
  negative_EVrelease <- subset(significant_EVrelease, estimate < 0)

  # 找出两者都是正相关或都是负相关的LR
  common_positive_LR <- intersect(positive_distance$LR, positive_EVrelease$LR)#正相关
  common_negative_LR <- intersect(negative_distance$LR, negative_EVrelease$LR)#负相关

  #
  cat("Number of LRs with positive correlation in both:", length(common_positive_LR), "\n")
  cat("Number of LRs with negative correlation in both:", length(common_negative_LR), "\n")

  inter.LR.results <- list()
  inter.LR.results[["common_positive_LR"]] <- common_positive_LR
  inter.LR.results[["common_negative_LR"]] <- common_negative_LR
  EV_spatalk_object@inter.LR.results <- inter.LR.results

  #求所有spot中存在的LR的个数
  #每一个spot中出现的>0的LR的个数
  set.seed(seeds)
  LR.list <- EV_spatalk_object@sender.spot.lr_interaction
  interaction_frequencies <- numeric(length(LR.list))
  # 遍历LR.list中的每个元素
  for (i in seq_along(LR.list)) {
    # 获取当前data.frame
    current_df <- LR.list[[i]]

    # 计算multiply大于0的行数
    interaction_frequencies[i] <- sum(current_df$multiply > 0)
  }
  # 创建一个data.frame来存储结果
  result_df <- data.frame(Spot = seq_along(LR.list), Frequency = interaction_frequencies)
  # 打印结果
  #print(result_df)
  # 创建一个函数来计算每个data frame中multiply大于0的次数
  count_interactions <- function(df) {
    sum(df$multiply > 0)
  }
  # 应用这个函数到LR.list中的每个元素
  interaction_counts <- sapply(LR.list, count_interactions)
  # 创建一个data frame来存储每个spot和对应的interaction次数
  spot.LR.freq.results <- cbind.data.frame(Spot = names(interaction_counts), Frequency = interaction_counts)
  EV_spatalk_object@spot.LR.freq.results <- spot.LR.freq.results


  #求所有LR在spot中的额频率和RRA #aggregateRanks {RobustRankAggreg}
  #每一LR在所有spot中的freq
  #>
  lr_interactions <- unique(LR.list[[1]]$interaction)
  # Initialize a named vector to store the frequency of non-zero 'multiply' values for each LR
  lr_frequencies <- setNames(rep(0, length(lr_interactions)), lr_interactions)

  # Function to increment frequency count for non-zero 'multiply' values, counting each interaction only once per spot

  # Calculate frequencies across all spots
  for (spot_id in names(LR.list)) {
    lr_frequencies <- increment_frequency(LR.list[[spot_id]], lr_frequencies)
  }

  # Create a data.frame for the frequencies
  LR.in.spot.frequency_table <- data.frame(LR = names(lr_frequencies), Freq.num = lr_frequencies, Frequency=lr_frequencies/length(LR.list))

  # Print the result
  #print(LR.in.spot.frequency_table)
  EV_spatalk_object@LR.in.spot.frequency_table <- as.data.frame(LR.in.spot.frequency_table)

  #>
  #RAA方法根据这些LR出现频率计算pvalue
  #首先将LR.list替换为计算RRA的list 命名为RRA.list, RRA.list的结构为list名字是每个spot的名字，里面的条目是interaction（就是LR的名字, 前面已经根据表达量排序好了）
  #然后应用RRA算法，对基因进行整合排序

  set.seed(seeds)
  RRA.list <- list()
  # Extract and store the sorted interactions for each spot
  for (spot_id in names(LR.list)) {
    # Assuming that higher 'multiply' values are better and have been pre-sorted in descending order
    RRA.list[[spot_id]] <- LR.list[[spot_id]]$interaction
  }
  LR.in.spot.RRA.results <- aggregateRanks(RRA.list, method="RRA", full = T, N=length(RRA.list))
  EV_spatalk_object@LR.in.spot.RRA.results <- as.data.frame(LR.in.spot.RRA.results)

  #整理EV_spatalk_stat_results 的统计结果，包括RRA，dist.corr，EV.release.corr
  Distance.corr.LR.results <- EV_spatalk_object@Distance.corr.LR.results
  rownames(Distance.corr.LR.results) <- Distance.corr.LR.results$LR
  colnames(Distance.corr.LR.results) <- paste0("dist_", colnames(Distance.corr.LR.results))
  Distance.corr.LR.results <- Distance.corr.LR.results[,c("dist_estimate", "dist_statistic", "dist_p.value")]

  EVrelease.corr.LR.results <- EV_spatalk_object@EVrelease.corr.LR.results
  rownames(EVrelease.corr.LR.results) <- EVrelease.corr.LR.results$LR
  colnames(EVrelease.corr.LR.results) <- paste0("EV.release_", colnames(EVrelease.corr.LR.results))
  EVrelease.corr.LR.results <- EVrelease.corr.LR.results[,c("EV.release_estimate", "EV.release_statistic", "EV.release_p.value")]

  colnames(LR.in.spot.RRA.results) <- paste0("RRA_", colnames(LR.in.spot.RRA.results))
  colnames(LR.in.spot.RRA.results)[1] <- "LR_pairs_ID"

  inter.rowname <- intersect(rownames(Distance.corr.LR.results), rownames(EVrelease.corr.LR.results))
  inter.rowname <- intersect(inter.rowname, rownames(LR.in.spot.RRA.results))

  EV_spatalk_stat_results <- cbind.data.frame(LR.in.spot.RRA.results[inter.rowname,], Distance.corr.LR.results[inter.rowname,], EVrelease.corr.LR.results[inter.rowname,])
  EV_spatalk_object@EV_spatalk_stat_results <- as.data.frame(EV_spatalk_stat_results)
  return(EV_spatalk_object)
}

#' Normalize and integrate interaction intensity of LR pairs
#'
#' Applies Min-Max normalization to the interaction intensity of each ligand-receptor (LR) pair, scaling the range to 0-1. Integrates the normalized results into the spatial transcriptomics (ST) Seurat object contained within the EV_spatalk_object.
#'
#' @param EV_spatalk_object An updated S4 object that has been prepared for EV_SpaTalk analysis.
#' @param mc.cores The number of cores to utilize for parallel processing. This should be adjusted according to the capabilities of your computational environment; set to 1 for single-core machines.
#' @return The EV_spatalk_object with added normalized interaction intensity scores for all LR pairs, incorporated into the metadata of the ST Seurat object.
#' @export
#' @examples
#' EV.spatalk.results <- add_interaction_score(EV_spatalk_object = EV.spatalk.results, mc.cores=30)
add_interaction_score <- function(EV_spatalk_object=EV.spatalk.results, mc.cores=30){
  #选择一个LR，可视化展示
  #例如"CD86_CTLA4"
  distance.correlation_results_df <- EV_spatalk_object@Distance.corr.LR.results
  EVrelease.correlation_results_df <- EV_spatalk_object@EVrelease.corr.LR.results
  common_LR.id <- c(EVrelease.correlation_results_df$LR, distance.correlation_results_df$LR)#这边是所有候选到的LR，没有取pvalue
  common_LR.id <- unique(common_LR.id)
  #提取显著common_positive_LR的LR
  #PVR_TIGIT
  #lr_interaction_list.raw$`AAACACCAATAACTGC-1`

  # Assuming `common_LR.id` is a character vector of the common ligand-receptor IDs
  # And `lr_interaction_list.raw` is your list of data frames

  # Initialize an empty matrix with rows as spot IDs and columns as common LR ids
  lr_interaction_list.raw <- EV_spatalk_object@all.spot.lr_interaction
  interaction_matrix <- matrix(NA, nrow = length(lr_interaction_list.raw), ncol = length(common_LR.id))
  rownames(interaction_matrix) <- names(lr_interaction_list.raw)
  colnames(interaction_matrix) <- common_LR.id

  # 获取核心数，设置要使用的核心数
  #no_cores <- detectCores() - 1

  # 使用mclapply并行处理每个spot
  results <- mclapply(names(lr_interaction_list.raw), process_spot, lr_list = lr_interaction_list.raw, lr_ids = common_LR.id, mc.cores = mc.cores)

  # 将结果转换为矩阵
  interaction_matrix <- do.call(rbind, results)

  # 转换为数据框
  interaction_df_raw.indensity <- as.data.frame(interaction_matrix)
  rownames(interaction_df_raw.indensity) <- names(lr_interaction_list.raw)
  EV_spatalk_object@interaction_df_raw.indensity <- interaction_df_raw.indensity

  interaction_df <- as.data.frame(interaction_matrix)
  interaction_df <- apply(interaction_df, 2, function(x){(x-min(x))/max(x)})

  # 如果需要，可以将rownames设置为spot IDs
  rownames(interaction_df) <- names(lr_interaction_list.raw)
  interaction_df[interaction_df=="NaN"] <- 0

  # Merge this data frame with the metadata of the Seurat object
  # Ensure that the rownames of your Seurat metadata match the spot IDs
  EV_spatalk_object@st.seurat.obj <- AddMetaData(EV_spatalk_object@st.seurat.obj, metadata = interaction_df)

  EV_spatalk_object@interaction_df =  as.data.frame(interaction_df)#增加每个spot的all.spotLR_interaction_score
  return(EV_spatalk_object)
}



