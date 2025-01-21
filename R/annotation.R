#  Perform weighted average scoring on the data based on the frequency of gene occurrence in the meta program.
#' Title
#'
#' @param object A seurat object of single cell data
#' @param meta_module_genes A list named after SSpMosaic meta modules,where each item is a character array storing the gene set corresponding to the meta module.
#' @param assay Name of assay to pull the gene expressiong data from.
#' @param layer Name of layer to fetch data from,for Seurat v5 only.
#' @param slot Name of slot to fetch data from,for Seurat v4 only.
#' @param normalize_method A string representing the normalization method, which should be selected from ‘none’, ‘zscore’, or ‘min-max’.
#' @param anno A boolean variable that specifies whether to annotate the cells with the highest-scoring meta-module.
#' @import SeuratObject
#' @import Seurat
#' @return A Seurat object, in which the metadata has been updated with the scoring values from the SSpMosaic meta module, assigning labels to each spot.
#' @export
#'
#'
weighted_score <- function(object,meta_module_genes ,assay = 'RNA',layer = 'counts',slot = 'counts',normalize_method = 'none',anno = F)
{
  table.list <- list()
  f = 0
  for(m in 1:length(names(meta_module_genes)))
  {
    n <- names(meta_module_genes)[m]
    ta <- table(meta_module_genes[[n]])
    table.list[[n]] <- ta
  }
  for(n in names(table.list))
  {
    table.list[[n]][table.list[[n]] == 1] = 0.5
  }
  if (IsSeuratV5()) {
    new_data <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  }else {
    new_data <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  }
  score_matrix <- data.frame(matrix(NA, nrow = length(colnames(new_data)), ncol = length(table.list)))
  rownames(score_matrix) <- colnames(new_data)
  colnames(score_matrix) <- names(table.list)
  for(n in names(table.list))
  {
    cat(n)
    cat('\n')
    gene <- names(table.list[[n]])
    gene <- intersect(gene,rownames(new_data))
    gm <- new_data[gene,]
    gm <- as.data.frame(gm)
    w <-  table.list[[n]][gene]
    gm$weight <- w
    w <- gm$weight
    gm <- dplyr::select(gm,-weight)
    result_matrix <- sweep(gm, 2, w, `*`)
    wm <- colMeans( result_matrix)
    score_matrix[,n] <- wm
  }
  data <- score_matrix
  if(normalize_method == 'zscore')
  {
    zs = scale(data)
    res <- zs
  }
  if(normalize_method == 'min-max')
  {
    min_vals <- apply(data, 2, min)
    max_vals <- apply(data, 2, max)
    normalized_data <- t(apply(data, 1, function(x) (x - min_vals) / (max_vals - min_vals)))
    res <- normalized_data
  }
  if(normalize_method == 'none')
  {
    res <- data
  }
  if(anno == T)
  {
    colnames(res) <- gsub("meta_program_","",colnames(res))
    res = as.data.frame(res)
    res$annotation <- apply(res, 1, function(x) colnames(res)[which.max(x)])
    res$annotation <- as.character(res$annotation)
    object <- SeuratObject::AddMetaData(object,res,col.name = colnames(res))
  }
  res = as.data.frame(res)
  object <- SeuratObject::AddMetaData(object,res,col.name = colnames(res))
  return(object)
}

#  Annotate the data based on the weighted average results.
#' Title
#'
#' @param object A Seurat object containing unsupervised clustering results
#' @param meta_module_genes A list named after SSpMosaic meta modules,where each item is a character array storing the gene set corresponding to the meta module.
#' @param cluster_col A string corresponding to one colname of Seurat object metadata,indicating the cell cluster assignment.
#' @param sd_thres A non-negative double,indicating the threshold of the difference between two standard deviations.
#' @param mean_thres A non-negative double,indicating the threshold of the difference between two mean values.
#' @param annotation_name A string,indicating the column in meta.data to store the annotation information.
#' @return A list,'object' is a Seurat object, in which the metadata has been updated with the SSpMosaic annotation;'mean_dif' is a list,in which the differences between mean values are stored;'sd_dif' is a list,in which the differences between standard deviations values are stored
#' @export
#'
#'
SSpMosaic_sc_annotation <- function(object,meta_module_genes,cluster_col,sd_thres = 0,mean_thres = 0,annotation_name = 'annotation'){
  clus = unique(object@meta.data[,cluster_col])
  pros = names(meta_module_genes)
  dif_list <- list()
  dif_list2 <- list()
  meta = object@meta.data
  for(i in clus){
    cat(i)
    cat('\n')
    dif_vec <-  c()
    dif_vec2 <-  c()
    ctv = meta[,cluster_col]
    ctv=  as.character(ctv)
    meta1 = meta[(ctv == i),]
    meta2 = meta[(ctv != i),]
    for(j in pros)
    {
      print(j)
      val1 = mean(meta1[,j])
      val2 = mean(meta2[,j])
      val3 <- sd(meta[,j])
      val4 <- sd(meta2[,j])
      differnce1 <- val1 - val2
      differnce2 <- val3 - val4
      dif_vec[j] <- differnce1
      dif_vec2[j] <- differnce2
    }
    names(dif_vec) <- pros
    names(dif_vec2) <- pros
    dif_list[[i]] = dif_vec
    dif_list2[[i]] = dif_vec2
  }
  for(i in names(dif_list2))
  {
    for(j in names(dif_list2[[i]]))
    {
      if(dif_list[[i]][j] < 0 )
      {
        dif_list2[[i]][j] = 0
      }
    }
  }
  max_matrix2 <- data.frame(matrix(NA, nrow = 1, ncol = length(clus)))
  rownames(max_matrix2) <- "cluster_id"
  colnames(max_matrix2) <- clus
  for(j in 1:length(colnames(max_matrix2)))
  {
    a <-  which.max(dif_list2[[colnames(max_matrix2)[j]]])
    if(dif_list2[[colnames(max_matrix2)[j]]][names(a)] > 0)
    {
      max_matrix2[1,colnames(max_matrix2)[j]] <- names(a)
    }
  }
  result <- max_matrix2
  result[1,] <- gsub("meta_program_","",result[1,])
  rownames(result)[1] <- "celltype_assignment"
  seurat_object = object
  meta_a = seurat_object@meta.data
  ctv = meta_a[,cluster_col]
  ctv=  as.character(ctv)
  meta_a[,annotation_name] = "NA"
  for(i in clus){
    cells = rownames(meta_a[(ctv == i),])
    meta_a[cells,annotation_name] = result[1, as.character(i)]
  }
  seurat_object@meta.data = meta_a
  res_list = list(object = seurat_object,mean_dif = dif_list,sd_dif = dif_list2)
  return(res_list)
}
