#' @Title read_txt_to_list read a txt file of gene sets and convert it into a list.
#'
#' @param filename Path to the txt file.
#'
#' @return A list named after the gene set, where each item is a character vector storing the genes contained in the corresponding gene set.
#'
#'
#'
read_txt_to_list <- function(filename) {
  con <- file(filename, "r")
  subset_list <- list()
  while (length(line <- readLines(con, n = 1)) > 0) {
    split_line <- strsplit(line, ",")[[1]]
    cluster_name <- split_line[1]
    genes_str <- split_line[2:length(split_line)]
    genes <- strsplit(genes_str, ",")
    genes <- as.character(genes)
    subset_list[[cluster_name]] <- genes
  }
  close(con)
  return(subset_list)
}

#  Perform the non-negative least-squares method
#' @Title .pred_prop
#'
#' @param x A matrix,gene expression matrix,rows are genes,cols are cells/spots.
#' @param mod A matrix,each SSpMosaic meta modules' gene frequency matrix,rows are genes,cols are SSpMosaic meta modules.
#' @import nnls
#' @return A matrix,rows are genes,cols are SSpMosaic meta modules.
#'
#'
#'
.pred_prop <- function(x, mod) {
  W <- mod
  x <- x[rownames(W), ]
  y <- vapply(
    seq_len(ncol(x)),
    function(i) {

      if (i %% 100 == 0) {
        cat("Completed", i, "columns\n")
      }
      nnls::nnls(W, x[, i])$x
    },
    numeric(ncol(W)))
  rownames(y) <- colnames(mod)
  colnames(y) <- colnames(x)
  return(y)
}

#  Read a set of samples’ SSpMosaic gene modules
#' @Title get_all_modules reads the result of filter_module.
#'
#' @param sample_names A character vector corresponding to a group of samples' name.
#' @param read_dir The output directory path of fuction filter_module.
#'
#' @return A list named after SSpMosaic meta modules,where each item is a character array storing the gene set corresponding to the meta module.
#' @export
#'
#'
get_all_modules <- function(sample_names,read_dir)
{
  sample_names = gsub("[- _/]", ".", sample_names)
  dirs = paste(read_dir,sample_names,"_program_chosen.txt",sep = "")
  sample_list <- list()
  for( i in dirs)
  {
    tryCatch({
      sample_list[[i]] = read_txt_to_list(i)
    },error = function(e){
      print(paste("Error occurred:", e))})
  }
  pr_list <- sample_list[[1]]
  if(length(sample_list) > 1)
  {
    for( i in 2:length(sample_list))
    {
      pr_list <- c(pr_list,sample_list[[i]])
    }
  }
  return(pr_list)
}

# Check if the current Seurat package is version 5.
IsSeuratV5 <- function ()
{
  startsWith(as.character(packageVersion("Seurat")), "5")
}

#  Perform deconvolution on spatial transcriptomics data using NNLS based on SSpMosaic gene modules.
#' Title
#'
#' @param seurat_object A seurat object of spatial transcriptomics data.
#' @param meta_moudule_genes A list named after SSpMosaic meta modules,where each item is a character array storing the gene set corresponding to the meta module.
#' @param normalize_method A string representing the normalization method, which should be selected from ‘none’, ‘zscore’, or ‘min-max’.
#' @param assay Name of assay to pull the gene expressiong data from.
#' @param layer Name of layer to fetch data from,for Seurat v5 only.
#' @param slot Name of slot to fetch data from,for Seurat v4 only.
#' @import SeuratObject
#' @import nnls
#' @import Seurat
#' @return A Seurat object, in which the metadata has been updated with the scoring values from the SSpMosaic meta module, assigning labels to each spot.
#' @export
#'
#'
SSpMosaic_decon <- function(seurat_object,meta_moudule_genes,normalize_method = 'zscore',assay = 'Spatial',layer = 'counts',slot = 'counts')
{
  table.list <- list()
  f = 0
  for(m in 1:length(names(meta_moudule_genes)))
  {
    n <- names(meta_moudule_genes)[m]
    ta <- table(meta_moudule_genes[[n]])
    table.list[[n]] <- ta
  }
  for(n in names(table.list))
  {
    table.list[[n]][table.list[[n]] == 1] = 0.5
  }
  if (IsSeuratV5()) {
    new_data <- SeuratObject::LayerData(seurat_object, assay = assay, layer = layer)
  }
  else {
    new_data <- Seurat::GetAssayData(seurat_object, assay = assay, slot = slot)}
  score_matrix <- data.frame(matrix(NA, nrow = length(colnames(new_data)), ncol = length(table.list)))
  rownames(score_matrix) <- colnames(new_data)
  colnames(score_matrix) <- names(table.list)
  pr_list = meta_moudule_genes
  gene_left <- c()
  for(m in pr_list)
  {
    gene_left <- c(gene_left,unique(m))
  }
  int_gene = intersect(rownames(new_data),gene_left)
  B = data.frame(matrix(0, nrow = length(pr_list), ncol = length(int_gene)))
  colnames(B) <- int_gene
  rownames(B) <- names(pr_list)
  A = new_data[rownames(new_data) %in% int_gene,]
  A = as.matrix(A)
  for(a in 1:length(table.list))
  {
    for(b in 1:length(table.list[[a]]))
    {
      B[names(table.list)[a],names(table.list[[a]][b])] <- table.list[[a]][b]
    }
  }
  B <- as.matrix(B)
  B <- t(B)
  B <- as.data.frame(B)
  B <- B[rownames(B)%in%rownames(A),]
  B <- as.matrix(B)
  res = .pred_prop(A,B)
  res = t(res)
  data <- res
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
  colnames(res) <- gsub("meta_program_","",colnames(res))
  res = as.data.frame(res)
  res$decon <- apply(res, 1, function(x) colnames(res)[which.max(x)])
  res$decon <- as.character(res$decon)
  seurat_object <- SeuratObject::AddMetaData(seurat_object,res,col.name = colnames(res))
  return(seurat_object)
}
