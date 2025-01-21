#  Read all SSpMosaic module candidates for a specific sample.
#' @Title get_module reads the result of generate_module.
#'
#' @param sample_name A string indicating the sample name ,corresponding to the input of generate_module.
#' @param read_dir The output directory path of fuction generate_module.
#' @param desc A boolean variable that determines whether to sort the genes within a module in descending order.
#' @import SeuratObject
#' @import Seurat
#' @import dplyr
#' @import hdWGCNA
#' @return A list named after gene programs,where each item is a character array storing the gene set corresponding to the gene program.
#' @export
#'
#'
get_module <- function(sample_name,read_dir,desc = FALSE)
{
  sample_name = gsub("[- _/]", ".", sample_name)
  index = paste(sample_name,'_',sep = '')
  dir_path <- read_dir
  file_names <- list.files(dir_path, full.names = FALSE)
  filtered_names1 <- file_names[grep(paste(index,".*_hdWGCNA_object\\.rds",sep= ""), file_names)]
  filtered_names1 <- paste(dir_path,filtered_names1,sep = '')
  merg_ <- list()
  for (nn in filtered_names1){
    seurat_obj <- readRDS(nn)
    celseq_ <- hdWGCNA::GetHubGenes(seurat_obj, n_hubs = 50)
    if(desc == TRUE)
    {
      celseq_ <- dplyr::arrange(celseq_ ,module, dplyr::desc(kME))
    }else{
      celseq_ <- dplyr::arrange(celseq_ ,module, kME)
    }
    unique_modules <- unique(celseq_$module)
    unique_modules <- unique_modules[unique_modules != 'grey']
    celseq_a <- list()
    for (module in unique_modules) {
      module_rows <- celseq_[celseq_$module == module, "gene_name"]
      celseq_a[[module]] <- c(module_rows)
    }
    merg_ <- c(merg_,celseq_a)
  }
  names(merg_) = gsub("[- _/]", ".",names(merg_))
  return(merg_)
}

#  Score all candidate modules in the sample data
#' Title
#'
#' @param module A list,the result of get_module,the input sample_name should be consistent.
#' @param sample_name A string indicating the sample name ,corresponding to the input of generate_module.
#' @param read_dir The output directory path of fuction generate_module.
#' @param out_dir The output directory of module score result.
#' @param cluster_col A string corresponding to one colname of Seurat object metadata,indicating the cell cluster assignment.
#' @param nbin Number of bins of aggregate expression levels for all analyzed features used in Seurat::AddModuleScore
#' @import utils
#' @import stringr
#' @import Seurat
#' @import SeuratObject
#' @return This fuction has no return.
#' @export
#'
#'
score_module <- function(module ,sample_name,read_dir,out_dir = NULL,cluster_col,nbin = 2)
{
  sample_name = gsub("[- _/]", ".", sample_name)
  index = paste(sample_name,'_',sep = '')
  dir_path <- read_dir
  if(is.null(out_dir) == T)
  {
    out_dir = read_dir
  }
  file_names <- list.files(dir_path, full.names = FALSE)
  filtered_names1 <- file_names[grep(paste(index,".*_hdWGCNA_object\\.rds",sep= ""), file_names)]
  filtered_names1 <- paste(dir_path,filtered_names1,sep = '')
  merg_ <- module
  if(is.null(out_dir) == T)
  {
    out_dir = read_dir
  }
  result_path <- out_dir
  for (nn in filtered_names1){
    seurat_obj <- readRDS(nn)
    v = stringr::str_extract(nn, "(?<=/)[^/]+(?=_hdWGCNA)")
    cat(v)
    cat('\n')
    meta_df = data.frame(matrix(NA, nrow = length(colnames(seurat_obj)), ncol = length(names(merg_))))
    colnames(meta_df) <- names(merg_)
    rownames(meta_df) <- colnames(seurat_obj)
    for (m in 1:length(names(merg_)))
    {
      cat(m)
      cat('\n')
      module_name <-names(merg_)[m]
      genes_in_module <- c(merg_[[module_name]])
      intersect_gene <- intersect(genes_in_module,rownames(seurat_obj))
      if(length(intersect_gene) == 1)
      {
        meta_df[, m] <- SeuratObject::FetchData(seurat_obj,vars = intersect_gene)
      }
      if(length(intersect_gene) >= 2)
      {
        intersect_list <- list()
        intersect_list[[1]] = intersect_gene
        seurat_obj <- Seurat::AddModuleScore(seurat_obj,features = intersect_list,name = "score",nbin = nbin)
        meta_df[, m] <- seurat_obj@meta.data[,"score1"]
      }
    }
    colnames(meta_df) <- gsub('\\+','.',colnames(meta_df))
    seurat_obj <- AddMetaData(seurat_obj,meta_df,col.name = colnames(meta_df))
    cols_to_keep <- c(colnames(meta_df),cluster_col)
    output_df <- seurat_obj@meta.data[,cols_to_keep,drop = F]
    utils::write.table(output_df, file = paste(result_path,stringr::str_extract(nn, "(?<=/)[^/]+(?=_hdWGCNA)"),'_program.txt',sep=''), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}

#' @Title extract_cell_type extracts celltype/cluster from the name of a certain gene program.
#'
#' @param prefix The prefix corresponding to a certain sample name
#' @param colName A string to extract celltype/cluster from,usually the name of a certain gene program
#'
#' @return A string corresponding to a certain celltype/cluster
#'
#'
#'
extract_cell_type <- function(prefix, colName) {
  if (substr(prefix, nchar(prefix), nchar(prefix)) == ".") {
    prefix <- substr(prefix, 1, nchar(prefix) - 1)
  }
  escapedPrefix <- gsub("\\.", "\\\\.", prefix)
  pattern <- paste0("^", escapedPrefix, "\\.(.+)\\.[0-9]+$")
  matches <- regexec(pattern, colName)
  result <- regmatches(colName, matches)
  if(length(result[[1]]) > 1) {
    cellType <- result[[1]][2]
  } else {
    cellType <- NA
  }
  return(cellType)
}

#' @Title write_list_to_txt writes a list to a txt file.
#'
#' @param subset_list A list named after gene programs,where each item is a character array storing the gene set corresponding to the gene program.
#' @param filename The name of the txt file.
#'
#' @return This fuction has no return.
#' @export
#'
#'
write_list_to_txt <- function(subset_list, filename) {
  con <- file(filename, "w")
  for (cluster_name in names(subset_list)) {
    genes <- as.character(subset_list[[cluster_name]])
    genes_str <- paste(genes, collapse = ",")
    writeLines(paste(cluster_name, ",", genes_str,sep = ""), con)
  }
  close(con)
}

#  Filter SSpMosaic candidate modules based on scoring results and output the filtered modules.
#' Title
#'
#' @param object A seurat object
#' @param module A list,the result of get_module
#' @param sample_name A string indicating the sample name ,corresponding to the input of generate_module&score_module.
#' @param read_dir The output directory path of fuction score_module.
#' @param out_dir The output directory of filtered module.
#' @param cluster_col A string corresponding to one colname of Seurat object metadata,indicating the cell cluster assignment.
#' @param sd_thres A non-negative double,indicating the threshold of the difference between two standard deviations.
#' @param mean_thres A non-negative double,indicating the threshold of the difference between two mean values.
#' @param pct_thres A non-negative double,indicating the threshold of the proportion of cells with program score greater than zero.
#' @param save_fig A boolean variable that specifies whether to save the visualization of the program score.
#' @param number_per_fig A positive integer that specifies the number of programs saving in one figure.
#' @param fig_dir A string indicating the directory path of the figures.
#' @param save_fig_spatial A boolean variable that specifies whether to save the spatial visualization of the program score.FALSE by default,only saving the UMAP visualization of the program score,TRUE,saving the spatial dimplot of the program score.
#' @param spot_size A positive double,indicating the spot size of the spatial dimplot.
#' @param merge_module A boolean variable that specifies whether to merge modules generated from the same clusters.
#' @import utils
#' @import Seurat
#' @import SeuratObject
#' @return A list,'chosen' is a dataframe,storing the information of the programs kept and their corresponding clusters;'obj' is a Seurat object, in which the metadata has been updated with the scoring values from the SSpMosaic programs;'sd'/'fc'/'pct' are lists named after clusters,where each item is a vector named after gene programs where each item is the program's metric value.
#' @export
#'
#'
filter_module <- function(object,module ,sample_name,read_dir,out_dir = NULL,cluster_col,sd_thres = 0.01,mean_thres = NA,pct_thres = NA,save_fig = FALSE,number_per_fig = 9,fig_dir = NULL,save_fig_spatial = FALSE,spot_size,merge_module = TRUE)
{
  sample_name = gsub("[- _/]", ".", sample_name)
  index = paste(sample_name,'_',sep = '')
  prefix = gsub('_','.',index )
  dir_path <- read_dir
  if(is.null(out_dir) == T)
  {
    out_dir = read_dir
  }
  file_names <- list.files(dir_path, full.names = FALSE)
  filtered_names <- file_names[grep(paste(index,".*_program\\.txt",sep= ""), file_names)]
  filtered_names <- paste(dir_path,filtered_names,sep = '')
  cat('Reading module score')
  cat('\n')
  mg = data.frame()
  feature.name <- names(module)
  meta.name <- cluster_col
  for(i in 1:length(filtered_names)){
    nn = filtered_names[i]
    meta_program_data <- utils::read.table(file = nn,sep = "\t")
    feature.count <- meta_program_data[,feature.name,drop = F]
    mg = rbind(mg,feature.count)
    new_meta <- meta_program_data[,meta.name,drop = F]
  }
  smart_seq2 <- object
  smart_seq2 <- smart_seq2[,colnames(smart_seq2) %in% rownames(mg)]
  smart_seq2 <- SeuratObject::AddMetaData(smart_seq2,mg,col.name = colnames(mg))
  smart_seq2@meta.data[,cluster_col] = gsub("[- _/]", ".", smart_seq2@meta.data[,cluster_col])
  act = unique(smart_seq2@meta.data[,cluster_col])
  max_matrix2 <- data.frame(matrix(NA, nrow = 1, ncol = length(colnames(mg))))
  rownames(max_matrix2) <- "cluster_id"
  colnames(max_matrix2) <- colnames(mg)
  for(pr in colnames(mg)){
    for(ct in act){
      if(extract_cell_type(prefix, pr) == ct)
      {
        max_matrix2[1,pr] <- ct
        break
      }
    }
  }
  max_matrix2 <- t(max_matrix2)
  max_matrix2 <- data.frame(max_matrix2)
  pro_list <- list()
  for(j in 1:length(act)){
    pro_list[[j]] <- subset(max_matrix2,cluster_id == act[j])
  }
  cat('Filtering module')
  cat('\n')
  dif_list <- list()
  dif_list2 <- list()
  dif_list3 <- list()
  for(j in 1:length(pro_list)){
    pro <- rownames(pro_list[[j]])
    clu <- pro_list[[j]][1,1]
    dif_vec <-  array()
    dif_vec2 <-  array()
    dif_vec3 <- array()
    meta = smart_seq2@meta.data
    ctv = meta[,cluster_col]
    ctv=  as.character(ctv)
    meta1 = meta[(ctv == clu),]
    meta2 = meta[(ctv != clu),]
    for(i in 1:length(pro)){
      n <- pro[i]
      print(n)
      val1 <- mean(meta1[,n])
      val3 <- sd(meta[,n])
      val2 <- mean(meta2[,n])
      val4 <- sd(meta2[,n])
      mr_array = meta1[,n]
      positive_count <- sum(mr_array > 0)/length(mr_array)
      differnce <- val1-val2
      differnce2 <- val3 - val4
      dif_vec[i] <- differnce
      dif_vec2[i] <- differnce2
      dif_vec3[i] <- positive_count
    }
    names(dif_vec) <- pro
    names(dif_vec2) <- pro
    names(dif_vec3) <- pro
    dif_list[[j]] = dif_vec
    dif_list2[[j]] = dif_vec2
    dif_list3[[j]] = dif_vec3
  }
  max_matrix <- max_matrix2
  if(is.na(sd_thres) == F)
  {
    for(j in 1:length(dif_list2)){
      for(l in names(dif_list2[[j]])){
        if(dif_list2[[j]][[l]] < sd_thres){
          max_matrix = subset(max_matrix,row.names(max_matrix) != l)
        }
      }
    }
  }
  if(is.na(mean_thres) == F)
  {
    for(j in 1:length(dif_list)){
      for(l in names(dif_list[[j]])){
        if(dif_list[[j]][[l]] < mean_thres){
          max_matrix = subset(max_matrix,row.names(max_matrix) != l)
        }
      }
    }
  }
  if(is.na(pct_thres) == F)
  {
    for(j in 1:length(dif_list3)){
      for(l in names(dif_list3[[j]])){
        if(dif_list3[[j]][[l]] < pct_thres){
          max_matrix = subset(max_matrix,row.names(max_matrix) != l)
        }
      }
    }
  }
  if(nrow(max_matrix) >0 )
  {
    if(merge_module == T)
    {
      meta_moudule <- list()
      for(ct in unique(max_matrix$cluster_id))
      {
        i = 1
        for(pr in rownames(max_matrix))
        {
          if(max_matrix[pr,1] == ct)
          {
            meta_moudule[[ct]][i] = pr
            i = i + 1
          }
        }
      }
      meta_moudule_genes <- list()
      for(meta in names(meta_moudule))
      {
        for(program in meta_moudule[[meta]])
        {
          meta_moudule_genes[[meta]] = c(meta_moudule_genes[[meta]],module[[program]])
        }
      }
      for(i in 1:length(names(meta_moudule_genes)))
      {
        names(meta_moudule_genes)[[i]] <- paste(index,names(meta_moudule_genes)[[i]],sep = "")
      }
      new_merg <- meta_moudule_genes
      cat('Saving module')
      cat('\n')
      write_list_to_txt(new_merg, paste(out_dir,index,"program_chosen.txt",sep = ''))
    }else{
      unmerge_module = list()
      for(pr in rownames(max_matrix))
      {
        unmerge_module[[pr]] = module[[pr]]
      }
      cat('Saving module without merging')
      cat('\n')
      write_list_to_txt(unmerge_module, paste(out_dir,index,"program_chosen.txt",sep = ''))
    }
    if(save_fig == T)
    {
      cat('Saving fig')
      cat('\n')
      if(is.null(fig_dir) == T)
      {
        fig_dir = out_dir
      }
      fig_path = fig_dir
      nr = ceiling(nrow(max_matrix)/number_per_fig)
      for(i in 1:nr)
      {
        p3 = Seurat::FeaturePlot(smart_seq2,features = rownames(max_matrix)[(number_per_fig*i-number_per_fig + 1):(number_per_fig*i)])
        ggsave(paste(fig_path,index,'_',(number_per_fig*i-number_per_fig + 1),"_",(number_per_fig*i),'program_distribution.png',sep=''),p3,width = 30,height = 40,dpi = 300)
        ggsave(paste(fig_path,index,'_',(number_per_fig*i-number_per_fig + 1),"_",(number_per_fig*i),'program_distribution.pdf',sep=''),p3,width = 30,height = 40,dpi = 300)
      }
    }
    if(save_fig_spatial == T)
    {
      cat('Saving spatail fig')
      cat('\n')
      if(is.null(fig_dir) == T)
      {
        fig_dir = out_dir
      }
      fig_path = fig_dir
      nr = ceiling(nrow(max_matrix)/number_per_fig)
      for(i in 1:nr)
      {
        p3 = Seurat::SpatialFeaturePlot(smart_seq2,features = rownames(max_matrix)[(number_per_fig*i-number_per_fig + 1):(number_per_fig*i)],pt.size.factor = spot_size,image.alpha = 0)
        ggsave(paste(fig_path,index,'_',(number_per_fig*i-number_per_fig + 1),"_",(number_per_fig*i),'program_spatial_distribution.png',sep=''),p3,width = 30,height = 40,dpi = 300)
        ggsave(paste(fig_path,index,'_',(number_per_fig*i-number_per_fig + 1),"_",(number_per_fig*i),'program_spatial_distribution.pdf',sep=''),p3,width = 30,height = 40,dpi = 300)
      }
    }
  }
  else{cat('No module were kept')
    cat('\n')}
  output = list(sd = dif_list2,fc = dif_list,pct = dif_list3,chosen = max_matrix,obj = smart_seq2)
  return(output)
}
