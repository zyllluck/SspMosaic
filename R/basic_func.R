#  Calculate PCA dimensionality reduction results for spatial transcriptomics data.
#' @Title st_pca_procedure perform the SCTransform and PCA reduction step for a Seurat object
#'
#' @param seurat_object A Seurat object
#' @param assay Name of assay to pull the count data from
#'
#' @return A Seurat object with the PCA calculation stored in the reductions slot and a new assay (named SCT) with counts being (corrected) counts, data being log1p(counts), scale.data being pearson residuals.
#' @import Seurat
#'
#'
st_pca_procedure <- function(seurat_object,assay = 'Spatial')
{
  sct_obj<-Seurat::SCTransform(seurat_object,assay = assay)
  sct_obj <- Seurat::RunPCA(sct_obj, assay = "SCT", verbose = FALSE)
  return(sct_obj)
}

#  Calculate PCA&UMAP dimensionality reduction results for single-cell transcriptomics data.
#' @Title pca_procedure perform standard procedure for a Seurat object
#'
#' @param merge_seurat A Seurat object
#'
#' @return A Seurat object with PCA and UMAP resuction
#' @import Seurat
#'
#'
pca_procedure <- function(merge_seurat)
{
  merge_seurat <- Seurat::NormalizeData(merge_seurat, verbose = FALSE)
  merge_seurat <- Seurat::FindVariableFeatures(merge_seurat, selection.method = "vst", verbose = FALSE)
  merge_seurat <- Seurat::ScaleData(merge_seurat, verbose = FALSE)
  merge_seurat <- Seurat::RunPCA(merge_seurat, npcs = 30, verbose = FALSE)
  merge_seurat <- Seurat::RunUMAP(merge_seurat, reduction = "pca", dims = 1:30)
  return(merge_seurat)
}

#  Count the occurrences of a specific string in a vector.
#' @Title Calculate the number of times a certain string appears in a character vector
#'
#' @param arr A character vector
#' @param target A certain string
#'
#' @return The number of times a certain string appears in a character vector
#'
#'
#'
count_str<- function(arr, target) {
  count <- base::sum(base::grepl(target, arr))
  return(count)
}


#  Calculate PCA&UMAP dimensionality reduction results for single-cell transcriptomics data.
#' @Title standard_procedure perform standard procedure for a Seurat object
#'
#' @param merge_seurat A Seurat object
#'
#' @return A Seurat object with PCA and UMAP resuction
#' @import Seurat
#'
#'
standard_procedure <- function(merge_seurat)
{
  merge_seurat <- NormalizeData(merge_seurat, verbose = FALSE)
  merge_seurat <- FindVariableFeatures(merge_seurat, selection.method = "vst", verbose = FALSE)
  merge_seurat <- ScaleData(merge_seurat, verbose = FALSE)
  merge_seurat <- RunPCA(merge_seurat, npcs = 30, verbose = FALSE)
  merge_seurat <- RunUMAP(merge_seurat, reduction = "pca", dims = 1:30)
  return(merge_seurat)
}


#  Read gene modules from a txt file.
#' @Title gene_list read a txt file of gene sets and convert it into a list.
#'
#' @param filename Path to the txt file.
#'
#' @return A list named after the gene set, where each item is a character vector storing the genes contained in the corresponding gene set.
#' @export
#'
#'
gene_list<- function(filename) {
  con <- file(filename, "r")
  subset_list <- list()
  while (length(line <- readLines(con, n = 1)) > 0) {
    split_line <- strsplit(line, ",")[[1]]
    cluster_name <- split_line[1]
    genes_str <- split_line[2:length(split_line)]
    genes <- strsplit(genes_str, ",")
    genes = as.character(genes)
    subset_list[[cluster_name]] <- genes
  }
  close(con)
  return(subset_list)
}


#' @Title Read a set of samples’ SSpMosaic gene modules
#'
#' @param sample_names A vector indicating a set of samples' names
#' @param read_dir Path to the SSpMosaic gene modules
#'
#' @return A list named after the sample names, where each item is a list that stores all the gene modules corresponding to the clusters for that specific sample.
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
      sample_list[[i]] = gene_list(i)
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

#' @Title Check if the current Seurat package is version 5.
#' @import utils
#' @return TRUE when the current Seurat package is version 5,FALSE when the current Seurat package is not version 5.
#'
#'
#'
IsSeuratV5 <- function ()
{
  startsWith(as.character(utils::packageVersion("Seurat")), "5")
}

#  Calculate the cosine distance between two vectors.
#' Title
#'
#' @param a A numeric vector.
#' @param b A numeric vector of same length.
#'
#' @return The cosine_distance of a,b.
#'
#'
#'
cosine_distance <- function(a, b) {
  sim = sum(a * b) / (sqrt(sum(a ^ 2)) * sqrt(sum(b ^ 2)))
  return(1-sim)
}
