#EV_spatalk other functions

#' Make equivalent random modules
#'
# @export
MakeRand = function(srt, db, assay = "SCT", nrand = 3, nbin = 5)
  {
  #if (is.null(assay)){
  #  assay = DefaultAssay(srt)
  #}
  #data = GetData(srt, slot = 'data')

  data = GetAssayData(srt, layer = "data")#change to V5 seurat get data
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m){
    lapply(1:10^nrand, function(i){
      used = vector()
      unused = binned
      for (g in db[[m]]){
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}

#' Modules to cells
#'
#' @keywords internal
GeneToEnrichment = function(
    srt,
    type = 'GO',
    db = NULL,
    method = 'rand',
    genes = NULL,
    assay = NULL,
    do.rescale = FALSE,
    min.cells = 0,
    min.genes = 0,
    min.var = 0,
    min.var.rescaled = 0,
    auc_percentile = 0.05,
    db_rand = NULL,
    nrand = 4,
    nbin = 5,
    ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(db)){
    db = FindMSigDB(type)
  }



  #counts = as.matrix(GetData(srt, assay = assay, slot = 'counts'))

  counts = as.matrix(GetAssayData(srt, assay = assay, layer = "counts"))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]

  if (method == 'metagene'){

    data = as.matrix(GetAssayData(srt, assay = assay, layer = 'scale.data'))

    db = lapply(db, intersect, genes.expr)

    enrichment.profile = t(sapply(names(db), function(m){
      colMeans(data[db[[m]], ], na.rm = TRUE)
    }))

    enrichment.profile = enrichment.profile[sapply(names(db), function(x){
      v = var(enrichment.profile[x, ])
      l = length(db[[x]])
      return(l > min.genes
             && v > min.var
             && v*l^2 > min.var.rescaled)
    }), ]

    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }

    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }

  if (method == 'auc'){

    data = as.matrix(GetAssayData(srt, assay = assay, layer = 'data'))

    cells_rankings = AUCell_buildRankings(data)
    cells_AUC = AUCell_calcAUC(db, cells_rankings, aucMaxRank=nrow(cells_rankings)*auc_percentile)
    enrichment.profile = getAUC(cells_AUC)

    if (do.rescale){

      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }

    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }

  if (method == 'score'){

    temp = AddModuleScore(srt, features = db, assay = assay, name = names(db), nbin = nbin, ...)

    enrichment.profile = t(temp@meta.data[, names(db)])

    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }

    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }

  if (method == 'rand'){

    data = as.matrix(GetAssayData(srt, assay = assay, layer = 'scale.data'))

    db = lapply(db, intersect, genes)

    if (is.null(db_rand)){
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    } else {
      nrand = log10(length(db_rand[[1]]))
    }

    enrichment.profile = t(sapply(names(db), function(m){
      ra = sapply(db_rand[[m]], function(i){
        colMeans(data[i, ], na.rm = TRUE)
      })
      re = colMeans(data[db[[m]], ], na.rm = TRUE)
      p = rowMeans(ra >= re)
      p = -log10(p)
      return(p)
    }))
    enrichment.profile[is.infinite(enrichment.profile)] = nrand
    enrichment.profile = enrichment.profile/nrand

    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }

  return(srt)
}


#' Find neighbors
#'
#' @export

FindSTNeighbors = function(
    st,
    d_max,
    d_min = 0
){
  #coord = st@images$slice1@coordinates[,c('imagerow','imagecol')]
  coord = st@images$image@coordinates[,c('imagerow','imagecol')]

  distances = as.matrix(dist(coord))
  distances = distances/sort(distances[distances > 0], decreasing = FALSE)[10]
  #distances = round(distances, digits = 1)
  neighbors = apply(distances, 1, function(d){
    names(d)[d >= d_min & d <= d_max]
  })
  names(neighbors) = rownames(coord)
  return(neighbors)
}

#' ST random and scale
#'
#' @export
MakeSTRand = function(
    st
){
  data = GetAssayData(st, layer = 'SCT', layer = 'data')
  data = t(apply(data, 1, function(row){
    sample(row, size = length(row), replace = FALSE)
  }))
  colnames(data) = colnames(st)
  st@assays$SCT@data = data
  st = ScaleData(st, assay = 'SCT', do.center = TRUE, do.scale = FALSE)
}

#' Screen EV-mediate interactome
#'
#' @export
rawParse2 <- function(data, top_genes = 50, stats = "mean")
{
  res = NULL
  cell_group <- unique(data$cell_type)
  pb <- progress::progress_bar$new(total = length(cell_group))
  pb$tick(0)
  for (i in cell_group) {
    sub_data <- data[data$cell_type == i, ]
    counts <- t(subset(sub_data, select = -cell_type))
    counts <- apply(counts, 2, function(x) {
      storage.mode(x) <- "numeric"
      x
    })
    if (stats == "mean") {
      temp <- data.frame(rowMeans(counts), i, stringsAsFactors = FALSE)
    }
    else if (stats == "median") {
      temp <- data.frame(apply(counts, 1, FUN = median),
                         i, stringsAsFactors = FALSE)
    }
    else {
      print("error stats option")
    }
    temp <- temp[order(temp[, 1], decreasing = TRUE), ]
    temp <- temp[1:ceiling(nrow(temp) * top_genes/100), ]
    temp <- temp %>% tibble::rownames_to_column()
    res <- rbind(res, temp)
    pb$tick()
  }
  colnames(res) <- c("gene", "exprs", "cell_type")
  return(res)
}


#' Find LR recurrented in sender spot
#'
#' @export
FindLR2 <- function(data_1, data_2 = NULL, datatype, comm_type, database = NULL)
{
  if (is.null(database)) {
    #database <- iTALK:::database
    #database <- EV.spatalk.results@database.new
    database <- EV.spatalk.results@database#这边是S4对象自带的database，即EV_spatalkdb


  }
  database <- database[database$Classification == comm_type,]
  if (datatype == "mean count") {
    gene_list_1 <- data_1
    if (is.null(data_2)) {
      gene_list_2 <- gene_list_1
    }
    else {
      gene_list_2 <- data_2
    }
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in%
                          gene_list_1$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in%
                            gene_list_2$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_1 <- database[ind, c("Ligand.ApprovedSymbol",
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_1[,
                                                                                           c("gene", "exprs", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_from_mean_exprs = exprs, cell_from = cell_type) %>%
      left_join(gene_list_2[, c("gene", "exprs", "cell_type")],
                by = c(Receptor.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_to_mean_exprs = exprs, cell_to = cell_type)
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in%
                          gene_list_2$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in%
                            gene_list_1$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_2 <- database[ind, c("Ligand.ApprovedSymbol",
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_2[,
                                                                                           c("gene", "exprs", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_from_mean_exprs = exprs, cell_from = cell_type) %>%
      left_join(gene_list_1[, c("gene", "exprs", "cell_type")],
                by = c(Receptor.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_to_mean_exprs = exprs, cell_to = cell_type)
    FilterTable <- rbind(FilterTable_1, FilterTable_2)
  }
  else if (datatype == "DEG") {
    gene_list_1 <- data_1
    if (is.null(data_2)) {
      gene_list_2 <- gene_list_1
    }
    else {
      gene_list_2 <- data_2
    }
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in%
                          gene_list_1$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in%
                            gene_list_2$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_1 <- database[ind, c("Ligand.ApprovedSymbol",
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_1[,
                                                                                           c("gene", "logFC", "q.value", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_from_logFC = logFC, cell_from_q.value = q.value,
                    cell_from = cell_type) %>% left_join(gene_list_2[,
                                                                     c("gene", "logFC", "q.value", "cell_type")], by = c(Receptor.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_to_logFC = logFC, cell_to_q.value = q.value,
                    cell_to = cell_type)
    ligand_ind <- which(database$Ligand.ApprovedSymbol %in%
                          gene_list_2$gene)
    receptor_ind <- which(database$Receptor.ApprovedSymbol %in%
                            gene_list_1$gene)
    ind <- intersect(ligand_ind, receptor_ind)
    FilterTable_2 <- database[ind, c("Ligand.ApprovedSymbol",
                                     "Receptor.ApprovedSymbol")] %>% left_join(gene_list_2[,
                                                                                           c("gene", "logFC", "q.value", "cell_type")], by = c(Ligand.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_from_logFC = logFC, cell_from_q.value = q.value,
                    cell_from = cell_type) %>% left_join(gene_list_1[,
                                                                     c("gene", "logFC", "q.value", "cell_type")], by = c(Receptor.ApprovedSymbol = "gene")) %>%
      dplyr::rename(cell_to_logFC = logFC, cell_to_q.value = q.value,
                    cell_to = cell_type)
    FilterTable <- rbind(FilterTable_1, FilterTable_2)
  }
  else {
    stop("Error: invalid data type")
  }
  FilterTable <- FilterTable[!duplicated(FilterTable), ]
  res <- as.data.frame(FilterTable) %>% dplyr::rename(ligand = Ligand.ApprovedSymbol,
                                                      receptor = Receptor.ApprovedSymbol)
  if (datatype == "DEG") {
    res <- res[!(res$cell_from_logFC == 1e-04 & res$cell_to_logFC ==
                   1e-04), ]
  }
  res <- res %>% mutate(comm_type = comm_type)
  return(res)
}


#' Calculate LR mean multiply
#'
#' @export
get_lr_multiply <- function(lr, lr_list) {
  sapply(lr_list, function(df) {
    if (any(df$interaction == lr)) {
      df$multiply[df$interaction == lr]
    } else {
      0  # 如果LR不存在于某个spot中，设为0
    }
  })
}

#' Calculate LR pseudo_count
#'
#' @export
calculate_interaction_strength <- function(L, R, pseudo_count = 0.01) {
  # 保证log2转换后不会出现负数
  pseudo_count_adjusted <- max(pseudo_count, 1)
  interaction_strength <- log2(L + pseudo_count_adjusted) + log2(R + pseudo_count_adjusted)
  return(interaction_strength)
}

#这步function可以放在functions_public.R
#' Include the non-zero LR
#'
#' @export
increment_frequency <- function(df, frequencies) {
  # Find the unique interactions with non-zero 'multiply' values
  unique_interactions <- unique(df$interaction[df$multiply > 0])
  # Increment the frequency count for these interactions
  if(length(unique_interactions) > 0) {
    frequencies[unique_interactions] <- frequencies[unique_interactions] + 1
  }
  return(frequencies)
}

#' Storage the LR interaction score
#'
#' @export
process_spot <- function(spot_id, lr_list, lr_ids) {
  spot_vector <- numeric(length(lr_ids)) # 初始化一个向量来存储multiply值
  names(spot_vector) <- lr_ids
  for (lr_id in lr_ids) {
    lr_row <- lr_list[[spot_id]][lr_list[[spot_id]]$interaction == lr_id, ]
    if (nrow(lr_row) == 1) {
      spot_vector[lr_id] <- lr_row$multiply
    } else {
      spot_vector[lr_id] <- NA # 如果没有找到相应的lr_id，赋值为NA
    }
  }
  return(spot_vector)
}


#作图
#' Spatial indensity_plot
#'
#' @export
LR_spatial_indensity_plot <- function(EV_spatalk_object, s.LR.pair) {
  # 将特征字符串分割为配体和受体
  st.seurat.obj <- EV_spatalk_object@st.seurat.obj

  #只考虑sender就可以
  sender.cell.id <- EV_spatalk_object@Sender.spot.id
  st.seurat.obj <- subset(st.seurat.obj, cells=sender.cell.id)

  feature_string <- s.LR.pair
  features_split <- strsplit(feature_string, "_")[[1]]
  ligand <- features_split[1]
  receptor <- features_split[2]

  print(paste("Sender cell number:", dim(st.seurat.obj)[2], sep = ""))

  # 创建副标题
  subtitle <- paste("Ligand:", ligand, "| Receptor:", receptor)

  # 生成SpatialFeaturePlot
  p <- SpatialFeaturePlot(st.seurat.obj, features = feature_string) +
    ggtitle(paste("Sender:", s.cell.type[1], "| Receiver:", s.cell.type[2])) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Helvetica", size = 14),
      legend.position = "right" # 设置图例在右边
    ) +
    labs(subtitle = subtitle, color = "Intensity") + # 添加副标题和修改图例标题
    theme(
      plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 12), # 设置副标题样式
      legend.title = element_text(size = 8) # 设置图例标题文字大小
    ) +
    guides(color = guide_legend(title = "Intensity")) # 修改图例标题

  return(p)
}

#' Pie plot of the statistical results of EV-release related or distance related LRs
#'
#' @export
LR_pie.plot <- function(EV_spatalk_object=EV.spatalk.results){
  distance.correlation_results_df <- EV.spatalk.results@Distance.corr.LR.results
  EVrelease.correlation_results_df <- EV.spatalk.results@EVrelease.corr.LR.results
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


  pie_data <- data.frame(
    Category = c("Positive Correlation", "Negative Correlation"),
    Count = c(length(common_positive_LR), length(common_negative_LR))
  )
  pie_data$percent <- pie_data$Count/sum(pie_data$Count)
  pie_data$label <- sprintf("%s \n#%d(%.1f%%)", pie_data$Category, pie_data$Count, pie_data$percent * 100)

  # 绘制饼图，并调整标签使其居中显示
  ggplot(pie_data, aes(x = "", y = Count, fill = Category, label = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 4,
              family = "sans",
              lineheight = 0.9,   # 调整行间距
              hjust = 0.5) +      # 水平居中对齐
    scale_fill_manual(values = c("Positive Correlation" = "#F47378", "Negative Correlation" = "#94D8F6")) +
    theme(legend.position = "right") +
    guides(fill = guide_legend(title = "Correlation Type")) +
    labs(fill = "Correlation Type")
}

#venn plot
#' Venn plot of the LRs overlapping
#'
#' @export
LR_venn.plot <- function(EV_spatalk_object=EV.spatalk.results){
  distance.correlation_results_df <- EV.spatalk.results@Distance.corr.LR.results
  EVrelease.correlation_results_df <- EV.spatalk.results@EVrelease.corr.LR.results
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

  #
  venn.plot <- venn.diagram(
    x = list(
      pos.Distance = positive_distance$LR,
      pos.EVrelease = positive_EVrelease$LR,
      neg.Distance = negative_distance$LR,
      neg.EVrelease = negative_EVrelease$LR

    ),
    category.names = c("Spat_dist+", "EV_release+", "Spat_dist-", "EV_release-"),
    output = F, # 修改这里为FALSE
    filename=NULL,
    fontfamily = "Arial",
    fontface = "bold",
    fill = c("#A349A4", "#138535", "gray90", "gray90"),
  )

  # 显示韦恩图
  grid.draw(venn.plot)
}

#circlize
#' circos plot of candidate LRs
#'
#' @export
LR_circos.plot <- function(EV_spatalk_object=EV.spatalk.results, select.LR.id=select.LR.id, specific.LR.id="CD96_PVR"){
  set.seed(2)
LR.freq.res <- EV_spatalk_object@LR.in.spot.frequency_table
LR.stat.res <- EV_spatalk_object@EV_spatalk_stat_results
inter.LR.id <- intersect(LR.freq.res$LR, LR.stat.res$LR_pairs_ID)
LR.circ.data <- cbind.data.frame(Frequency=LR.freq.res[inter.LR.id,][,-2], LR.stat.res[inter.LR.id,][,-1])

#if(length(select.LR.id)>=0){
common_positive_LR <- EV_spatalk_object@inter.LR.results$common_positive_LR
#  common_negative_LR <- EV_spatalk_object@inter.LR.results$common_negative_LR
#  LR.circ.data.1 <- LR.circ.data[match(c(common_positive_LR, common_negative_LR), LR.circ.data$Frequency.LR),]
#}
LR.circ.data.1 <- LR.circ.data[select.LR.id,]
LR.circ.data.1$RRA_Score <- ifelse(LR.circ.data.1$RRA_Score<0.05, 1, 0)
LR.circ.data.1$dist_p.value <- ifelse(LR.circ.data.1$dist_p.value<0.05, 1, 0)
LR.circ.data.1$EV.release_p.value <- ifelse(LR.circ.data.1$EV.release_p.value<0.05, 1, 0)

LR.circ.data.1 <- cbind.data.frame(group="non.sig", LR.circ.data.1)

LR.circ.data.1$group[LR.circ.data.1$Frequency.LR %in% common_positive_LR] <- "sig"
#LR.circ.data.1$group[which(LR.circ.data.1$RRA_Score==0)] <- "non.sig"
LR.circ.data.1 <- LR.circ.data.1[order(LR.circ.data.1$Frequency.Frequency, decreasing = F),]
LR.circ.data.1 <- LR.circ.data.1[order(LR.circ.data.1$Frequency.LR, decreasing = F),]
LR.circ.data.1 <- LR.circ.data.1[order(LR.circ.data.1$group, decreasing = T),]


num_rows <- nrow(LR.circ.data.1)
# Create the sequences for 'start' and 'end' columns
start_seq <- seq(1, by=2, length.out=num_rows)
end_seq <- start_seq + 2
order.df <- cbind.data.frame(start_seq, end_seq)
# Add the sequences as columns to your data frame
LR.circ.data.2 <- cbind.data.frame(group=LR.circ.data.1[,1], order.df, LR.circ.data.1[,-c(1)])
LR.circ.data.2$dist.upper <- LR.circ.data.2$dist_estimate
LR.circ.data.2$dist.lower <- LR.circ.data.2$dist_estimate
LR.circ.data.2$EV.upper <- LR.circ.data.2$EV.release_estimate
LR.circ.data.2$EV.lower <- LR.circ.data.2$EV.release_estimate

LR.circ.data.2$dist.upper <- ifelse(LR.circ.data.2$dist.upper>0, LR.circ.data.2$dist.upper, 0)
LR.circ.data.2$dist.lower <- ifelse(LR.circ.data.2$dist.lower<0, LR.circ.data.2$dist.lower*1, 0)
LR.circ.data.2$EV.upper <- ifelse(LR.circ.data.2$EV.upper>0, LR.circ.data.2$EV.upper, 0)
LR.circ.data.2$EV.lower <- ifelse(LR.circ.data.2$EV.lower<0, LR.circ.data.2$EV.lower*1, 0)

circos.rawdata.1  <- LR.circ.data.2
#
rand_col = function(k) {
  return(rgb(runif(k), runif(k), runif(k)))
}
posTransform.fun = function(region) {
  return(region)
}
bed <- as.data.frame(circos.rawdata.1[,1:3])
circos.clear()
circos.par("cell.padding"=c(0,0,0,0))
circos.initializeWithIdeogram(bed, plotType = NA)
bed <- circos.rawdata.1[,1:3]
bed <- cbind.data.frame(bed, rownames(bed))
colnames(bed)[4] <- 'value1'
bed$color <- circos.rawdata.1$group
bed$color <- ifelse(bed$color == "sig", "black", "gray80")  # color1 和 color2 替换成你选择的颜色
circos.genomicLabels(bed, labels.column = 4, side="outside", niceFacing = TRUE, col = bed$color, cex = 0.4, connection_height = 0.01, font = 4)

#RAA
bed <- cbind.data.frame(circos.rawdata.1[,c(1,2,3,6)])
colnames(bed)[4] <- c('value1')
for (i in 4) {
  exp <- bed[,i]
  bed[,i] <- (exp-min(exp))/max(exp)
}
f = colorRamp2(breaks = c(summary(bed$value1)[1], summary(bed$value1)[6]), colors = c("gray80", "#FFB27D"))
circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = f(value[[1]]),
                     border = 'black', lwd = 0.1, posTransform = posTransform.default, ...)
}, bg.border = NA, track.height = 0.05)

#freq
bed_list = circos.rawdata.1[,c(1,2,3,5)]
circos.genomicTrackPlotRegion(bed_list, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "s", col = "#377EB8", cex = 0.8, lwd = 2, pt.col = "#FB8072", pch = 19, boder = NA ,...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  #circos.lines(cell.xlim, c(1, 1), col = "gray70", lty = 3)
}, track.height = 0.13, bg.border = 'black', bg.col='gray95')


#hr
bed_list = list(circos.rawdata.1[,c(1,2,3,16)],
                circos.rawdata.1[,c(1,2,3,15)])
col = c("#367EB8", "#FB8C62")
circos.genomicTrackPlotRegion(bed_list, ylim = c(min(circos.rawdata.1[,16]), max(circos.rawdata.1[,15])), panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicLines(region, value, type = "h", col = col[i], baseline = 0, lwd = 1.8, boder = NA ,...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  circos.lines(cell.xlim, c(0, 0), col = "gray70", lty = 3)
}, track.height = 0.15, bg.border = 'black', bg.col='gray95')


#highlight sector
s.im.type.id <- names(summary(as.factor(circos.rawdata.1$group)))
highlight.sector(s.im.type.id[1], track.index=c(1), col = "#00FF0000", border = "black", lwd=0.3)
highlight.sector(s.im.type.id[2], track.index=c(1), col = "#00FF0000", border = "black", lwd=0.3)


#LR之间的相关性
if(!is.na(specific.LR.id))
{
  s.lr.indensity.score <- EV_spatalk_object@interaction_df_raw.indensity[EV_spatalk_object@Sender.spot.id, circos.rawdata.1$Frequency.LR]

  cor_pmat_with_specific <- function(data, specific_col) {
    # Get the index of the specific column
    specific_col_index <- which(colnames(data) == specific_col)

    # Error handling if the specific column is not found
    if (length(specific_col_index) == 0) {
      stop("The specific column does not exist in the data frame.")
    }

    # Calculate the correlation and p-values against the specific column
    cor_pvalues <- apply(data, 2, function(x) {
      cor.test(x, data[, specific_col_index])$p.value
    })

    # Return the p-values
    return(cor_pvalues)
  }


  # Now you can use cor_pmat_with_specific() to get the p-values for correlations with the specific column
  p_values_with_specific <- cor_pmat_with_specific(s.lr.indensity.score, specific.LR.id)

  # Combine the correlations and p-values
  combined <- data.frame(correlation = cor_matrix[, specific.LR.id], p_value = p_values_with_specific)

  # Filter for significant correlations (p < 0.05)
  significant <- combined[combined$p_value < 0.05, ]

  # Sort by the absolute value of correlation in descending order
  significant <- significant[significant$correlation>0,]
  significant <- significant[order(-abs(significant$correlation)), ]

  # If there are less than 10, take as many as there are; otherwise, take the top 10
  top_correlations <- head(significant, min(11, nrow(significant)))

  # Print the top correlations
  print(top_correlations)

  colnames(bed) <- c("chr", "start", "end", "value")
  region2 <- bed[c(specific.LR.id),]
  region1 <- bed[c(setdiff(rownames(top_correlations), specific.LR.id)),]

  for (i in 1:nrow(region1)) {
    circos.link(sector.index1 = region1[i, 1], point1 = region1[i, 2],
                sector.index2 = region2[1], point2 = region2[2],
                col = "#FB8C62", lwd = 1.5, directional = -1)
  }
}

}

