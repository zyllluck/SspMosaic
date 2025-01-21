#' @Title getSparseMatrixData is used to handle sparse matrices to speed up the integration process.
#'
#' @param sparse_matrix A dgCMatrix object
#'
#' @return A list containing three items: values, which includes all non-zero values of the sparse matrix; rows, which contains the corresponding row names; and cols, which contains the corresponding column names.
#'
#'
#'
getSparseMatrixData <- function(sparse_matrix) {
  if (!inherits(sparse_matrix, "dgCMatrix")) {
    stop("The provided matrix is not a 'dgCMatrix' object.")
  }
  values <- sparse_matrix@x
  rows_indices <- sparse_matrix@i + 1
  col_pointers <- sparse_matrix@p
  cols_indices <- rep(1:ncol(sparse_matrix), times = diff(col_pointers))
  rownames_matrix <- rownames(sparse_matrix)
  colnames_matrix <- colnames(sparse_matrix)
  if (is.null(rownames_matrix)) {
    rownames_matrix <- as.character(seq_len(nrow(sparse_matrix)))
  }
  if (is.null(colnames_matrix)) {
    colnames_matrix <- as.character(seq_len(ncol(sparse_matrix)))
  }
  rows <- rownames_matrix[rows_indices]
  cols <- colnames_matrix[cols_indices]
  result <- list(
    values = values,
    rows = rows,
    cols = cols
  )
  return(result)
}

#' @Title find_group is used to obtain the modality/species information for each batch.
#'
#' @param search_string A batch name
#' @param named_list A list named after all the modalities or species contained in the dataset, where each item is a character array storing the names of all batches corresponding to that modality or species.
#'
#' @return The modality/species this batch belongs to.
#'
#'
#'
find_group <- function(search_string, named_list) {
  for (name in names(named_list)) {
    if (search_string %in% named_list[[name]]) {
      return(name)
    }
  }
  return(NULL)
}

#' @Title locate_meta_neighbors finds the nearest neighbors for each SSpMosaic meta based on PCA reduction
#'
#' @param pca A dataframe , rows are cells and cols are PCs.
#' @param batch_list A vector indicating the batch to which each cell belongs.
#' @param celltype_list A vector indicating the cluster to which each cell belongs.
#' @param n_meta_neighbors A positive integer indicating the number of meta neighbors.
#' @param UNASSIGN The temporary label name for cells that have not been assigned a tag.
#'
#' @return A list named after the names of the clusters(SSpMosaic meta),where each item is a character array storing the names of its meta neighbors
#' @import dplyr
#' @import stats
#'
locate_meta_neighbors <- function(pca, batch_list, celltype_list, n_meta_neighbors, UNASSIGN = 'UNASSIGN') {
  # 1. Get all unique cell types
  celltypes <- unique(celltype_list)
  # Initialize array to hold cell type neighbors
  cns <- vector("list", length(celltypes))
  names(cns) <- celltypes
  # Handle case when n_meta_neighbors is 1 or there is only one cell type
  if (n_meta_neighbors == 1 || length(celltypes) == 1) {
    cns <- stats::setNames(lapply(celltypes, function(celltype) {
      return(c(celltype))
    }), celltypes)
  }
  else {
    # 2. Compute average PCA per cell type per batch
    c2b <- lapply(celltypes, function(celltype) {
      return(unique(batch_list[celltype_list == celltype]))
    })
    names(c2b) <- celltypes
    b2d <- list()
    for (batch in unique(batch_list)) {
      batch_flag <- batch_list == batch
      sub_celltypes <- unique(celltype_list[batch_flag])
      pcs <- do.call(rbind, lapply(sub_celltypes, function(sub_celltype) {
        pca_subset <- pca[batch_flag & celltype_list == sub_celltype, , drop = FALSE]
        colMeans(pca_subset)
      }))
      dists <- as.matrix(dist(pcs))
      b2d[[batch]] <- dists
      rownames(b2d[[batch]]) <- sub_celltypes
      colnames(b2d[[batch]]) <- sub_celltypes
    }
    # 3. For each batch, sort by distance and choose top n neighbors
    for (celltype in celltypes) {
      celltype = as.character(celltype)
      neighbors <- unique(unlist(lapply(c2b[[celltype]], function(b) {
        sort(b2d[[b]][celltype,])[1:n_meta_neighbors]
      })) %>% names())
      neighbors <- neighbors[neighbors != '']
      cns[[celltype]] <- neighbors
    }
  }
  # 4. Assign neighbors to special UNASSIGN cell type, if present
  if (UNASSIGN %in% celltypes) {
    cns <- lapply(names(cns), function(cns_key) {
      if (cns_key == UNASSIGN) {
        return(celltypes)
      } else {
        return(union(cns[[cns_key]], UNASSIGN))
      }
    })
    names(cns) <- names(cns)
  }
  return(cns)
}

#' @Title query_annoy_tree query an annoy tree,return a group of cell's KNN result.
#'
#' @param query a dataframe stores a group of cell's PCA reduction,rows are cells,cols are PCs
#' @param ckd an annoy tree
#' @param k the number of cell neighbors
#'
#' @return A list with two data frame.index is a dataframe,rows are cells,cols are the neighbors'index.dist is a dataframe,rows are cells,cols are the neighbors'distance
#' @import RcppAnnoy
#'
#'
query_annoy_tree <- function(query, ckd, k = 3) {
  index <- matrix(nrow = nrow(x = query), ncol = k)
  dist <- matrix(nrow = nrow(x = query), ncol = k)
  for (i in seq(nrow(x = query))) {
    holder <- ckd$getNNsByVectorList(
      (query[i, ] %>% as.numeric()),
      n = k,
      search_k = -1,
      include_distances = TRUE
    )
    index[i, ] <- holder$item
    dist[i, ] <- holder$distance
  }
  return(list(index = index + 1L, dist = dist))
}

#  Run the KNN algorithm.
#' @Title
#'
#' @param pca A dataframe stores the pca reduction of cells
#' @param batch_vec A vector stores the bacth information of cells
#' @param seed A integer , Random seed
#' @param n_pcs A vector of positive integer indicating the PCs used for calculating distance
#' @param neighbors_within_modality A positive integer,the number of cell neighbors from the same modality/species
#' @param neighbors_across_modality A positive integer,the number of cell neighbors from the different modality/species
#' @param annoy_n_trees A positive integer,the number of annoy trees
#' @param batch_group_by_modality A list named by modality/species names, containing the batch names belonging to this modality/species.
#' @import methods
#' @return A sparseMatrix,rows are cells,cols are cells,for each row,the cols of the corresponding cell's nearest neighbors stores the distance,other cols are zero
#'
#'
#'
generate_knn_dis_mat <- function(pca,batch_vec,seed = 42,n_pcs = 1:30,neighbors_within_modality = 3,neighbors_across_modality = 3,annoy_n_trees = 10L,batch_group_by_modality = NULL){
  metric = "euclidean"
  batch_list <- as.character(x = batch_vec)
  batch_counts <- table(batch_list)
  batches <- names(x = batch_counts)
  if(is.null(batch_group_by_modality) == TRUE)
  {
    neighbors_within_batch <- neighbors_within_modality
  }
  if(is.null(batch_group_by_modality) == FALSE)
  {
    neighbors_within_batch <- max(c(neighbors_within_modality,neighbors_across_modality))
  }
  differ <- abs(neighbors_within_modality - neighbors_across_modality)
  k <- neighbors_within_batch * length(x = batches)
  min_batch <- min(batch_counts)
  if (min_batch < neighbors_within_batch) {
    b <- batches[which.min(x = batch_counts)]
    stop(
      "Not all batches have at least 'neighbors_within_batch' cells in them: ",
      "\n neighbors_within_batch: ", neighbors_within_batch,
      "\n cells of ", b, ": ", min_batch
    )
  }
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  knn_dist <- matrix(0, nrow = nrow(x = pca), ncol = k)
  knn_index <- matrix(0L, nrow = nrow(x = pca), ncol = k)
  message('Finding nearest neighbor for each cell in each batch')
  for (to_ind in seq_along(along.with = batches)) {
    message(paste('Calclulating batch ',to_ind,sep = ''))
    batch_to <- batches[to_ind]
    if(is.null(batch_group_by_modality) == FALSE)
    {
      to_modality <- find_group(batch_to,batch_group_by_modality)
    }
    mask_to <- batch_list == batch_to
    ind_to <- seq(length(x = batch_list))[mask_to]
    data <- pca[mask_to, n_pcs]
    ckd <- switch(
      EXPR = metric,
      "euclidean" = methods::new(RcppAnnoy::AnnoyEuclidean, ncol(x = data)),
      "angular" = methods::new(RcppAnnoy::AnnoyAngular, ncol(x = data)),
      "manhattan" = methods::new(RcppAnnoy::AnnoyManhattan, ncol(x = data)),
      "hamming" = methods::new(RcppAnnoy::AnnoyHamming, ncol(x = data)),
    )
    ckd$setSeed(seed)
    for (i in seq(nrow(x = data))) {
      ckd$addItem(i - 1, (data[i, ] %>% as.numeric()))
    }
    ckd$build(annoy_n_trees)
    for (from_ind in seq_along(along.with = batches)) {
      batch_from <- batches[from_ind]
      if(is.null(batch_group_by_modality) == FALSE)
      {
        from_modality <- find_group(batch_from,batch_group_by_modality)
      }
      mask_from <- batch_list == batch_from
      ind_from <- seq(length(x = batch_list))[mask_from]
      ckdout <- query_annoy_tree(
        query = pca[mask_from, n_pcs],
        ckd = ckd,
        k = neighbors_within_batch
      )
      col_range <- seq(
        (to_ind - 1) * neighbors_within_batch + 1,
        to_ind * neighbors_within_batch
      )
      try1 <- ckdout$index
      try2 <- ckdout$dist
      for(l in 1:nrow(try1))
      {
        try1[l,] = ind_to[try1[l,]]
      }
      knn_index[ind_from,col_range]  = try1
      knn_dist[ind_from,col_range]  = try2
      if(is.null(batch_group_by_modality) == FALSE)
      {
        if(from_modality == to_modality)
        {
          if(differ > 0 )
          {
            zero_range <- seq((max(col_range) - differ + 1),max(col_range))
            knn_dist[ind_from,zero_range]  = 0
          }
        }
      }
    }
    ckd$unload()
  }
  for (i in seq(nrow(x = knn_index))) {
    knn_index[i, ] <- knn_index[i, ][order(knn_dist[i, ])]
    knn_dist[i, ] <- knn_dist[i, ][order(knn_dist[i, ])]
  }
  message('Finding nearest neighbor finished,generating distance matrix')
  n_obs = nrow(x = knn_index)
  n_neighbors = ncol(x = knn_index)
  if (!is.null(x = seed)) {
    set.seed(seed)
  }
  if(n_obs > 150000)
  {
    row_index <- c()
    col_index <- c()
    m_val <- c()
    for(i in 1:n_obs)
    {
      row_index <- c(row_index,as.numeric(knn_index[i, ]))
      col_index <- c(col_index,rep(i,length(knn_index[i, ])))
      m_val <- c(m_val,as.numeric(knn_dist[i, ]))
    }
    coln <- rownames(pca)
    rown <- rownames(pca)
    dis_mat <- sparseMatrix(row_index, col_index, x=m_val, dims=c(length(rown), length(coln)),dimnames=list(rown, coln))
  }
  if(n_obs <= 150000)
  {
    dis_mat <- matrix(0, nrow = n_obs, ncol = n_obs)
    for(i in 1:n_obs)
    {
      dis_mat[knn_index[i, ],i] <- knn_dist[i, ]
    }
    rownames(dis_mat) <- rownames(pca)
    colnames(dis_mat) <- rownames(pca)
    dis_mat <- methods::as(dis_mat,'sparseMatrix')
  }
  return(dis_mat)
}

#  Perform data integration using SSpMosaic meta labels.
#' Title
#'
#' @param object A Seurat object
#' @param col_batch one of the colnames of the object's metadata,storing the batch information of the cells
#' @param col_SSpMosaic_cluster one of the colnames of the object's metadata,storing the batch information of the cells
#' @param reduction_use a string indicating the name of the reductio result used for calculating distance
#' @param reduction_dimention A vector of positive integer indicating the PCs used for calculating distance
#' @param seed a integer,random seed
#' @param neighbors_within_modality A positive integer,the number of cell neighbors from the same modality/species
#' @param neighbors_across_modality A positive integer,the number of cell neighbors from the different modality/species
#' @param graph_name a string indicating the name used to save the SSpMosaic integration result
#' @param n_meta_neighbors A positive integer indicating the number of meta neighbors.
#' @param batch_group_by_modality A list named by modality/species names, containing the batch names belonging to this modality/species.
#' @param calculate_umap A boolean variable that specifies whether to compute the UMAP dimensionality reduction results after integration.
#' @import dplyr
#' @import methods
#' @import uwot
#' @import SeuratObject
#' @import Seurat
#' @return A Seurat object after integration
#' @export
#'
#'
SSpMosaic_meta_integration.Seurat <- function(object,col_batch,col_SSpMosaic_cluster = 'SSpMosaic_cluster',reduction_use = 'pca',reduction_dimention = 1:30,seed = 42,neighbors_within_modality = 3,neighbors_across_modality = 3,graph_name = "SSpMosaic",n_meta_neighbors,batch_group_by_modality = NULL,calculate_umap = T)
{
  meta = object@meta.data
  if(is.null(batch_group_by_modality) == TRUE)
  {
    neighbors_within_batch <- neighbors_within_modality
  }
  if(is.null(batch_group_by_modality) == FALSE)
  {
    neighbors_within_batch <- max(c(neighbors_within_modality,neighbors_across_modality))
  }
  column_index <- which(names(meta) == col_SSpMosaic_cluster)
  meta[, column_index] <- as.character(meta[, column_index])
  embeddings <- SeuratObject::Embeddings(object = object, reduction = reduction_use)
  batch_vec <- SeuratObject::FetchData(object = object, vars = col_batch)
  batch_num  <-  unique(batch_vec[,1]) %>% length()
  tn <- neighbors_within_batch * batch_num * 10
  SSpMosaic_cluster_vec <- as.character(FetchData(object = object, vars = col_SSpMosaic_cluster)[, 1])
  SSpMosaic_cluster <- unique(SSpMosaic_cluster_vec)
  SSpMosaic_cluster.list <- list()
  for(i in SSpMosaic_cluster)
  {
    subset_df <- meta[meta[, column_index] == i, ]
    SSpMosaic_cluster.list[[i]] <- rownames(subset_df)
  }
  meta_neighbors = locate_meta_neighbors(pca = embeddings[,reduction_dimention],batch_list = batch_vec[,1],celltype_list = SSpMosaic_cluster_vec,n_meta_neighbors = n_meta_neighbors)
  meta_cluster <- list()
  for(i in names(meta_neighbors))
  {
    cts <- meta_neighbors[[i]]
    cell_names <- c()
    for(j in cts)
    {
      cell_names <- c(cell_names,SSpMosaic_cluster.list[[j]])
    }
    meta_cluster[[i]] <- cell_names
  }
  each_cluster_res <- list()
  for(i in names(meta_cluster))
  {
    message(paste('Working on SSpMosaic cluster',i,sep = ''))
    chosen_col <- SSpMosaic_cluster.list[[i]]
    embeddings_use <- embeddings[meta_cluster[[i]],]
    batch_vec_use <- batch_vec[meta_cluster[[i]],1]
    res_for_each_cluster = generate_knn_dis_mat(pca = embeddings_use,batch_vec = batch_vec_use,seed = seed,neighbors_within_modality = neighbors_within_modality,neighbors_across_modality = neighbors_across_modality,batch_group_by_modality = batch_group_by_modality,n_pcs = reduction_dimention)
    sparse_res <- res_for_each_cluster[,chosen_col]
    each_cluster_res[[i]] <- sparse_res %>% getSparseMatrixData()
    message('\n')
  }
  message('Merging the results from the SSpMosaic clusters')
  row <- c()
  col <- c()
  m_value <- c()
  for(i in 1:length(each_cluster_res))
  {
    message(i)
    row <- c(row,each_cluster_res[[i]][['rows']])
    col <- c(col,each_cluster_res[[i]][['cols']])
    m_value <- c(m_value,each_cluster_res[[i]][['values']])
  }
  coln <- rownames(embeddings)
  rown <- rownames(embeddings)
  col_index <- match(col,coln)
  row_index <- match(row,rown)
  dis_mat <- sparseMatrix(row_index, col_index, x=m_value, dims=c(length(rown), length(coln)),dimnames=list(rown, coln))
  message('Converting distance matrix to connectivity matrix')
  X <- matrix(0, nrow = nrow(embeddings), ncol = 2)
  set_op_mix_ratio = 1
  local_connectivity = 1
  cnts <- uwot::similarity_graph(
    X = X,
    nn_method = dis_mat,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity
  )
  rownames(cnts) <- rownames(embeddings)
  colnames(cnts) <- rownames(embeddings)
  assay <- NULL
  assay <- assay %||% DefaultAssay(object = object)
  if(calculate_umap == T)
  {
    message('Computing umap reduction')
    min_dist = 0.3
    spread = 1
    UMAP_key = "UMAP_"
    UMAP_name = 'umap'
    umap <- uwot::umap( X = embeddings,
                  n_threads = 1,
                  n_neighbors = (3*4),
                  nn_method = dis_mat,
                  set_op_mix_ratio = set_op_mix_ratio,
                  local_connectivity = local_connectivity,
                  min_dist = min_dist,
                  spread = spread)
    colnames(x = umap) <- paste0(UMAP_key, c(1, 2))
    rownames(x = umap) <- rownames(x = cnts)
    object[[UMAP_name]] <- SeuratObject::CreateDimReducObject(
      embeddings = umap,
      assay = assay,
      key = UMAP_key
    )
  }
  message('Adding the knn graph to seurat object')
  graph <- as.Graph(x = cnts)
  methods::slot(object = graph, name = "assay.used") <- assay
  object[[graph_name]] <- graph
  message('All done')
  return(object)
}
