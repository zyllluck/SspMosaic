# `new_ReassignModules` is a modification of `hdWGCNA::ReassignModules`, which fixes the bug that occurs when there is only one gene module.
#' Title
#'
#' @param seurat_obj A Seurat object
#' @param harmonized logical indicating whether or not to use harmonized MEs
#' @param features character vector containing features for manual module reassignment
#' @param new_modules character vector containing modules to reassign the genes
#' @param ignore logical indicating whether or not to ignore error message about reassigning non-grey features
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @import dplyr
#' @import hdWGCNA
#' @return A Seurat object with the an updated modules table for the selected wgcna experiment
#'
#'
#'
new_ReassignModules <- function (seurat_obj, harmonized = TRUE, features = NULL, new_modules = NULL,
                                 ignore = FALSE, wgcna_name = NULL)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- hdWGCNA::GetModules(seurat_obj, wgcna_name)
  MEs <- hdWGCNA::GetMEs(seurat_obj, harmonized, wgcna_name)
  mod_colors <- dplyr::select(modules, c(module, color)) %>%
    distinct()
  mod_cp <- mod_colors$color
  names(mod_cp) <- as.character(mod_colors$module)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  genes_use <- hdWGCNA::GetWGCNAGenes(seurat_obj, wgcna_name)
  if (!is.null(features) & !is.null(new_modules)) {
    if (!all(features %in% genes_use)) {
      stop("Some features are not found in GetWGCNAGenes(seurat_obj).")
    }
    if (!all(new_modules %in% levels(modules$module))) {
      stop("Some module names in new_modules are invalid. Features must be reassigned to existing modules found in GetModules(seurat_obj)")
    }
    orig_mods <- subset(modules, gene_name %in% features) %>%
      .$module %>% as.character()
    if (!all(orig_mods == "grey") & !ignore) {
      stop("Attempting to reassign non-grey genes to new modules. Proceed with caution. If you wish to reassign these genes, re-run this function and set ignore=TRUE")
    }
    reassign_df <- data.frame(gene_name = as.character(features),
                              module = as.character(new_modules))
    reassign_df$module <- factor(as.character(reassign_df$module),
                                 levels = levels(modules$module))
    reassign_df$color <- as.character(mod_cp[as.character(reassign_df$module)])
    modules[reassign_df$gene_name, "module"] <- reassign_df$module
    modules[reassign_df$gene_name, "color"] <- reassign_df$color
  }else {
    neg_df <- do.call(rbind, lapply(mods, function(cur_mod) {
      cur <- subset(modules, module == cur_mod)
      cur <- cur[, c("gene_name", "module", paste0("kME_",
                                                   cur_mod))]
      names(cur)[3] <- "kME"
      cur %>% subset(kME < 0)
    }))
    if (nrow(neg_df) == 0) {
      return(seurat_obj)
    }
    rownames(neg_df) <- 1:nrow(neg_df)
    kMEs <- modules[, 4:ncol(modules),drop = F]
    kMEs <- kMEs[, colnames(kMEs) != "kME_grey",drop = F]
    reassigned <- sapply(neg_df$gene_name, function(cur_gene) {
      cur_kMEs <- kMEs[cur_gene, ]
      max_kME <- max(cur_kMEs)
      if (max_kME < 0) {
        return("kME_grey")
      }
      colnames(kMEs)[which(cur_kMEs == max_kME)]
    })
    reassigned <- do.call(rbind, strsplit(reassigned, "kME_"))[,
                                                               2]
    reassigned <- factor(as.character(reassigned), levels = levels(modules$module))
    reassigned_colors <- as.character(mod_cp[as.character(reassigned)])
    modules[neg_df$gene_name, "module"] <- reassigned
    modules[neg_df$gene_name, "color"] <- reassigned_colors
  }
  seurat_obj <- hdWGCNA::SetModules(seurat_obj, modules, wgcna_name)
  seurat_obj
}

# `new_ModuleConnectivity` is a modification of `hdWGCNA::ModuleConnectivity`, where the internally called `hdWGCNA::ReassignModules` has been replaced with `new_ReassignModules`.
#' Title
#'
#' @param seurat_obj A Seurat object
#' @param group.by column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
#' @param group_name name of the group(s) in group.by to use for kME calculation
#' @param corFnc character string specifying the function to be used to calculate co-expression similarity.
#' @param corOptions character string specifying additional arguments to be passed to the function given by corFnc
#' @param harmonized logical determining whether to use harmonized MEs for kME calculation
#' @param assay Assay in seurat_obj containing expression information
#' @param slot Slot in seurat_obj, default to normalized 'data' slot
#' @param sparse logical indicating whether or not to run the correlation using a sparse matrix
#' @param reassign_modules logical indicating whether or not to reassign genes to different co-expression modules if their kME value in the assigned module is negative
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @import qlcMatrix
#' @import hdWGCNA
#' @import utils
#' @import WGCNA
#' @import dplyr
#' @import SeuratObject
#' @return A Seurat object with the kMEs computed for the selected wgcna experiment
#'
#'
#'
new_ModuleConnectivity <- function (seurat_obj, group.by = NULL, group_name = NULL, corFnc = "bicor",
                                    corOptions = "use='p'", harmonized = TRUE, assay = NULL,
                                    slot = "data", sparse = TRUE, reassign_modules = TRUE, wgcna_name = NULL,
                                    ...)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- hdWGCNA::GetModules(seurat_obj, wgcna_name)
  MEs <- hdWGCNA::GetMEs(seurat_obj, harmonized, wgcna_name)
  genes_use <- as.character(modules$gene_name)
  params <- hdWGCNA::GetWGCNAParams(seurat_obj, wgcna_name)
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }
  if (!is.null(group.by)) {
    cells.use <- seurat_obj@meta.data %>% subset(get(group.by) %in%
                                                   group_name) %>% rownames
    MEs <- MEs[cells.use, ]
  }
  else {
    cells.use <- colnames(seurat_obj)
  }
  exp_mat <- SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot)[genes_use,
                                                                  cells.use]
  if (sparse) {
    if (!("qlcMatrix" %in% utils::installed.packages())) {
      stop("Need to install package qlcMatrix if sparse=TRUE. Either set sparse=FALSE or install qlcMatrix package.")
    }
    kMEs <- qlcMatrix::corSparse(X = Matrix::t(exp_mat),
                                 Y = as.matrix(MEs))
    rownames(kMEs) <- genes_use
    kMEs <- as.data.frame(kMEs)
  }else {
    datExpr <- t(as.matrix(exp_mat))
    kMEs <- WGCNA::signedKME(datExpr, MEs, outputColumnName = "kME",
                             corFnc = corFnc, corOptions = corOptions, ...)
  }
  modules <- modules[, 1:3]
  mods <- levels(modules$module)
  colnames(kMEs) <- colnames(MEs)
  kMEs <- kMEs[, mods]
  colnames(kMEs) <- paste0("kME_", colnames(kMEs))
  kMEs <- cbind(modules, kMEs)
  seurat_obj <- hdWGCNA::SetModules(seurat_obj, kMEs, wgcna_name)
  if (reassign_modules) {
    seurat_obj <- new_ReassignModules(seurat_obj, harmonized = harmonized,
                                      wgcna_name = wgcna_name)
  }
  seurat_obj
}

#  Generate a set of SSpMosaic gene module candidates for each cluster.
#' @Title generate_module
#'
#' @param object A Seurat object with PCA reduction and cell cluster information
#' @param cluster_col A string corresponding to one colname of Seurat object metadata,indicating the cell cluster assignment
#' @param meta_cell A named integer vector,named by cell clusters,indicating the number of meta cell in hdWGCNA for each cluster,or NULL to use the default value.
#' @param soft_power A named integer vector,named by cell clusters,indicating the soft_power in hdWGCNA for each cluster,or NULL to use the default value.
#' @param cluster_chosen A chatacter vector ,specify the cell clusters to generate SSpMosaic modules,or NULL to choose all clusters
#' @param min_cell_number A positive integer,indicating the min cell number for a cell cluster
#' @param log2FC_thres A positive double,indicating log2FC threshold when finding cell cluster markers
#' @param sample_name A string,indicating sample name
#' @param out_dir A string,indicating the output path
#' @param normalize_metacell A named boolean vector that specifies whether to normalize metacell in hdWGCNA for each cluster,named by cell clusters,or NULL to use the default value.
#' @param min_metacell A positive integer,indicating the min meta_cell number
#' @param min.pct A non-negative double,only test genes that are detected in a minimum fraction of min.pct cells in the cell cluster.
#' @param max_share A named integer vector,named by cell clusters,indicating the number of the maximum number of shared cells for meta cells in hdWGCNA for each cluster,or NULL to use the default value
#' @param assay A named vector,named by cell clusters,indicating the assay used for hdWGCNA.
#' @param slot A named vector,named by cell clusters,indicating slot to extract data.
#' @param layer A named vector,named by cell clusters,indicating layer to extract data,only used for seurat v5.
#' @param gene_use A named list,indicating the gene candidates used for SSpMosaic,or NULL to use marker genes.
#' @param verbose Whether to print detailed information during the generation of the module.
#' @import Seurat
#' @import patchwork
#' @import ggplot2
#' @import grDevices
#' @import hdWGCNA
#' @return This function save modules to out_dir,without returning.
#' @export
#'
#'
generate_module <- function(object,cluster_col,meta_cell = NULL,soft_power = NULL,cluster_chosen = NULL,min_cell_number = 25,log2FC_thres = 0.25,sample_name,out_dir,normalize_metacell = NULL,min_metacell = 6,min.pct = NULL,max_share = NULL,assay = NULL,slot = NULL,layer = NULL,gene_use = NULL,verbose = TRUE)
{
  batch.list <- list()
  batch.list[[1]] = object
  sample_name = gsub("[- _/]", ".", sample_name)
  index = paste(sample_name,'_',sep = '')
  result_path = out_dir
  setwd(result_path)
  batch.list[[1]]@meta.data[,cluster_col] = gsub("[- _/]", ".", batch.list[[1]]@meta.data[,cluster_col])
  Idents(batch.list[[1]]) <- batch.list[[1]]@meta.data[,cluster_col]
  wrong_id <- ""
  cat('Filtering clusters')
  cat('\n')
  cell_counts = batch.list[[1]]@meta.data[,cluster_col] %>% table()
  keep_clusters <- names(cell_counts[cell_counts >= min_cell_number])
  batch.list[[1]] <- subset(batch.list[[1]], idents = keep_clusters)
  drop_clusters <- names(cell_counts[cell_counts < min_cell_number])
  cat(paste('Cluster',drop_clusters,'removed',sep = ' '))
  cat('\n')
  if(is.null(min.pct))
  {
    if(IsSeuratV5())
    {
      min.pct = 0.01
    }else{
      min.pct = 0.1
    }
  }
  for(i in 1:1)
  {
    cat('Spliting clusters')
    cat('\n')
    cluster.list <- SplitObject(batch.list[[i]], split.by = cluster_col)
    cluster <- names(cluster.list)
    if(is.null(cluster_chosen))
    {
      cluster_chosen = cluster
    }
    else{
      cluster_chosen = gsub("[- _/]", ".", cluster_chosen)
    }

    if(is.null(meta_cell) == T)
    {
      meta_cell = cell_counts[keep_clusters]
      meta_cell = ceiling(meta_cell/30)
      meta_cell = meta_cell[cluster_chosen]
    }
    else{
      names(meta_cell) = gsub("[- _/]", ".", names(meta_cell))
      meta_cell = meta_cell[!names(meta_cell) %in% drop_clusters]
      if(setequal(names(meta_cell) ,cluster_chosen) == F)
      {
        stop("meta_cell names must be consistent with cluster_chosen!", call. = FALSE)
      }
      if(is.numeric(meta_cell) == F )
      {
        stop("meta_cell must be numeric!", call. = FALSE)
      }
      for(mcn in names(meta_cell))
      {
        if(is.na(meta_cell[mcn]) == T)
        {
          meta_cell[mcn] = ceiling(cell_counts[mcn]/30)
        }
      }
    }
    if(is.null(normalize_metacell) == T)
    {
      normalize_metacell = rep(NA,length(cluster_chosen))
      names(normalize_metacell) = cluster_chosen
    }
    else{
      names(normalize_metacell) = gsub("[- _/]", ".", names(normalize_metacell))
      normalize_metacell = normalize_metacell[!names(normalize_metacell) %in% drop_clusters]
      if(setequal(names(normalize_metacell) ,cluster_chosen) == F)
      {
        stop("normalize_metacell names must be consistent with cluster_chosen!", call. = FALSE)
      }
      if(all(is.na(normalize_metacell) | normalize_metacell %in% c(TRUE, FALSE)) == F )
      {
        stop("normalize_metacell must be bool!", call. = FALSE)
      }
    }
    if(is.null(soft_power) == T)
    {
      soft_power = rep(NA,length(cluster_chosen))
      names(soft_power) = cluster_chosen
    }
    else{
      names(soft_power) = gsub("[- _/]", ".", names(soft_power))
      if(setequal(names(soft_power) ,cluster_chosen) == F)
      {
        stop("soft_power names must be consistent with cluster_chosen!", call. = FALSE)
      }
      if(is.numeric(soft_power) == F )
      {
        stop("soft_power must be numeric!", call. = FALSE)
      }
    }
    if(is.null(max_share) == T)
    {
      max_share = rep(NA,length(cluster_chosen))
      names(max_share) = cluster_chosen
    }
    else{
      names(max_share) = gsub("[- _/]", ".", names(max_share))
      if(setequal(names(max_share) ,cluster_chosen) == F)
      {
        stop("max_share names must be consistent with cluster_chosen!", call. = FALSE)
      }
      if(is.numeric(max_share) == F )
      {
        stop("max_share must be numeric!", call. = FALSE)
      }
    }
    if(is.null(slot) == T)
    {
      slot = rep('data',length(cluster_chosen))
      names(slot) = cluster_chosen
    }
    else{
      names(slot) = gsub("[- _/]", ".", names(slot))
      if(setequal(names(slot) ,cluster_chosen) == F)
      {
        stop("slot names must be consistent with cluster_chosen!", call. = FALSE)
      }
    }
    if(is.null(layer) == T)
    {
      layer = rep('counts',length(cluster_chosen))
      names(layer) = cluster_chosen
    }
    else{
      names(layer) = gsub("[- _/]", ".", names(layer))
      if(setequal(names(layer) ,cluster_chosen) == F)
      {
        stop("layer names must be consistent with cluster_chosen!", call. = FALSE)
      }
    }
    if(is.null(assay) == T)
    {
      assay = rep(NA,length(cluster_chosen))
      names(assay) = cluster_chosen
    }else{
      names(assay) = gsub("[- _/]", ".", names(assay))
    }
    if(is.null(gene_use) == T)
    {
      gene_use = list()
    }else{
      names(gene_use) = gsub("[- _/]", ".", names(gene_use))
    }
    for(j in cluster_chosen[cluster_chosen %in% keep_clusters])
    {
      tryCatch({
        cat(j)
        cat('\n')
        if(is.null(gene_use[[j]])){
          cat('Finding markers')
          cat('\n')
        markers <- Seurat::FindMarkers(batch.list[[i]],ident.1 = j,only.pos = T, logfc.threshold = log2FC_thres,min.pct = min.pct)
        markers <- subset(markers, avg_log2FC > 0)
        markers <- rownames(markers)
        }else{
          markers <- gene_use[[j]]
        }
        seurat_obj <- cluster.list[[j]]
        seurat_obj <- hdWGCNA::SetupForWGCNA(
          seurat_obj,
          features = markers,
          wgcna_name = "tutorial"
        )
        metacell_number <-  meta_cell[j]#ceiling(cell_number/30)
        if(metacell_number < min_metacell)
          metacell_number = min_metacell
        if(is.na(max_share[j]))
        {
          max_sharec <- ceiling(metacell_number/2)
        }else{
          max_sharec <- max_share[j]
        }
        cat('Constructing metacell')
        cat('\n')
        seurat_obj <- hdWGCNA::MetacellsByGroups(
          min_cells = min_cell_number,
          seurat_obj = seurat_obj,
          group.by = cluster_col,
          reduction = "pca",
          k = metacell_number,
          max_shared =  max_sharec,
          ident.group = cluster_col
        )
        cell_number = cell_counts[j]
        if(is.na(normalize_metacell[j]) == T){
          if(cell_number > 1000)
          {
            cat('Normalizing metacell')
            cat('\n')
            seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
          }
        }else{
          if(normalize_metacell[j] == T)
          {
            cat('Normalizing metacell')
            cat('\n')
            seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
          }
        }
        if(verbose == TRUE){
          seurat_obj@misc[["tutorial"]][["wgcna_metacell_obj"]]

          head(seurat_obj@misc$tutorial$wgcna_metacell_obj, 2)
        }
        if(is.na(assay[j])){
        seurat_obj <- hdWGCNA::SetDatExpr(
          seurat_obj,
          assay = DefaultAssay(seurat_obj),
          slot = slot[j],layer = layer[j]
        )}else{
          seurat_obj <- hdWGCNA::SetDatExpr(
            seurat_obj,
            assay = assay[j],
            slot = slot[j],layer = layer[j]
          )
        }
        if(verbose == TRUE){
        seurat_obj <- hdWGCNA::TestSoftPowers(
          seurat_obj,powers = c(seq(1,10,by = 1),seq(12,50,by = 2)),
          networkType = "signed"
        )
        plot_list <- hdWGCNA::PlotSoftPowers(seurat_obj)
        }else{
          tmp_file <- tempfile()
          sink(tmp_file)
          seurat_obj <- hdWGCNA::TestSoftPowers(
            seurat_obj,powers = c(seq(1,10,by = 1),seq(12,50,by = 2)),
            networkType = "signed"
          )
          plot_list <- hdWGCNA::PlotSoftPowers(seurat_obj)
          sink()
          unlink(tmp_file)
      }
        p <- patchwork::wrap_plots(plot_list, ncol=2)
        ggplot2::ggsave(paste('./',index,j,'_wrap_plots.pdf',sep=''),p)
        cat('Constructing network')
        cat('\n')
        if(is.na(soft_power[j]) == T){
          seurat_obj <- hdWGCNA::ConstructNetwork(
            seurat_obj,
            setDatExpr = FALSE,
            tom_name =paste(index,j,'___',sep = ''),
            minModuleSize=50,
            mergeCutHeight=0.1,
            deepSplit=4,
            maxModuleSize=300
            ,overwrite = TRUE
          )
        }
        else{
          seurat_obj <- hdWGCNA::ConstructNetwork(
            seurat_obj,
            soft_power=soft_power[j],
            setDatExpr = FALSE,
            tom_name =paste(index,j,'___',sep = ''),
            minModuleSize=50,
            mergeCutHeight=0.1,
            deepSplit=4,
            maxModuleSize=300
            ,overwrite = TRUE
          )
        }
        grDevices::pdf(paste('./',index,j,'_PlotDendrogram.pdf',sep='') )
        hdWGCNA::PlotDendrogram(seurat_obj, main=paste(index,j,'hdWGCNA Dendrogram'))
        grDevices::dev.off()
        TOM <- GetTOM(seurat_obj)
        if(verbose == TRUE){
          seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj),do.irlba = FALSE)
        seurat_obj <- hdWGCNA::ModuleEigengenes(
          seurat_obj,npcs = 3
        )
        seurat_obj <- new_ModuleConnectivity(
          seurat_obj,
          group.by = NULL, group_name = NULL)
        }else{
          tmp_file <- tempfile()
          sink(tmp_file)
          suppressMessages({
          seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj),do.irlba = FALSE)
          seurat_obj <- hdWGCNA::ModuleEigengenes(
            seurat_obj,npcs = 3
          )
          seurat_obj <- new_ModuleConnectivity(
            seurat_obj,
            group.by = NULL, group_name = NULL)
          })
          sink()
          unlink(tmp_file)
        }
        seurat_obj <- hdWGCNA::ResetModuleNames(
          seurat_obj,
          new_name = paste(index,j,'_',sep='')
        )
        select <- dplyr::select
        p <- hdWGCNA::PlotKMEs(seurat_obj, ncol = 5, n_hubs = 10)
        ggsave(paste('./',index,j,'_plot_kmes.pdf',sep=''),p)
        modules <- hdWGCNA::GetModules(seurat_obj)
        write.table(modules, file = paste('./',index,j,'_hdWGCNA_modules.txt',sep = ''), sep = "\t", quote = FALSE, row.names = TRUE)
        hub_df <- hdWGCNA::GetHubGenes(seurat_obj = seurat_obj, n_hubs = 10)
        module_names <- levels(modules$module)
        filtered_module_names <- module_names[module_names != "grey"]
        module_data <- data.frame(matrix(NA, nrow = length(filtered_module_names), ncol = ncol(seurat_obj)))
        rownames(module_data) <- filtered_module_names
        colnames(module_data) <- colnames(seurat_obj)
        for (m in 1:length(filtered_module_names)) {
          module_name <-filtered_module_names[m]
          genes_in_module <- rownames(modules[modules$module == module_name, ])
          mgdf = as.data.frame(t(FetchData(seurat_obj,vars = genes_in_module) ))
          module_data[m, ] <- colMeans(mgdf, na.rm = TRUE)
        }
        module_data <- t(module_data)
        write.table(module_data, file = paste('./',index,j,'_hdWGCNA_module_data.txt',sep = ''), sep = "\t", quote = FALSE, row.names = TRUE)
        seurat_obj <- AddMetaData(seurat_obj, module_data, col.name = colnames(module_data))
        cat('Saving object')
        cat('\n')
        cat('\n')
        saveRDS(seurat_obj, file = paste('./',index,j,'_hdWGCNA_object.rds',sep = ''))
      },error = function(e){
        cat(paste("Error occurred:", e))
        wrong_id <- paste(wrong_id,j,sep = "+")})
    }
  }
}
