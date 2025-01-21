#  Normalize the adjacency matrix of the input biological network
#' Title
#'
#' @param graph A igraph object,biological network,such as PPI network.
#' @param conserve_heat A boolean variable that specifies whether to maintain heat conservation during network propagation.
#' @param weighted A boolean variable that specifies whether the biological network is weighted.
#' @import Matrix
#' @import igraph
#' @return A matrix,normalized adjacency matrix
#' @export
#'
#'
get_normalized_adjacency_matrix <- function(graph, conserve_heat=TRUE, weighted=FALSE) {
  stopifnot("graph cannot have nodes with degree=zero" = all(igraph::degree(graph) != 0))
  if (conserve_heat) {
    graph_weighted <- igraph::make_empty_graph(directed = TRUE)
  } else {
    graph_weighted <- igraph::make_empty_graph(directed = FALSE)
  }
  vertices <- igraph::V(graph)
  graph_weighted <- igraph::add_vertices(graph_weighted, nv = length(vertices))
  edge_weights <- c()
  node_to_degree <- igraph::degree(graph)
  if (weighted && !igraph::is_weighted(graph)) {
    warning("Input graph is not weighted. All edge weights will be set to 1.")
  }
  if (conserve_heat == TRUE && weighted == FALSE){
    mat = igraph::as_adjacency_matrix(graph, attr="weight")
    mat@x = ifelse(mat@x != 0, 1, 0)
    row_sums <- Matrix::rowSums(mat)
    w_prime = mat / row_sums
  }else{
    for (e in igraph::E(graph)) {
      v1 <- igraph::ends(graph, e)[1]
      v2 <- igraph::ends(graph, e)[2]
      deg1 <- node_to_degree[v1]
      deg2 <- node_to_degree[v2]
      weight <- if (weighted) igraph::E(graph)[e]$weight else 1
      if (conserve_heat) {
        edge_weights <- rbind(edge_weights, c(v1, v2, weight / deg1))
        edge_weights <- rbind(edge_weights, c(v2, v1, weight / deg2))
      } else {
        edge_weights <- rbind(edge_weights, c(v1, v2, weight / sqrt(deg1 * deg2)))
      }
    }
    for (i in 1:nrow(edge_weights)) {
      graph_weighted <- igraph::add_edges(graph_weighted, c(edge_weights[i, 1], edge_weights[i, 2]), attr=list(weight=as.numeric(edge_weights[i, 3])))
    }
    w_prime <- igraph::as_adjacency_matrix(graph_weighted, attr="weight")
  }
  return(w_prime)
}

#  Calculate the contributions of each individual gene to the final heat of each other gene  after propagation
#' Title
#'
#' @param nam_or_graph A igraph object or a normalized_adjacency_matrix
#' @param alpha A double number,between 0  and 1
#' @param conserve_heat A boolean variable that specifies whether to maintain heat conservation during network propagation.
#' @param weighted A boolean variable that specifies whether the biological network is weighted.
#' @import Matrix
#' @return A matrix,individual_heats_matrix
#' @export
#'
#'
get_individual_heats_matrix <- function(nam_or_graph, alpha=0.5, conserve_heat=TRUE, weighted=FALSE) {
  stopifnot("Alpha must be between 0 and 1" = 1 >= alpha && alpha >= 0)
  if (is.igraph(nam_or_graph)) {
    nam <- get_normalized_adjacency_matrix(nam_or_graph, conserve_heat=conserve_heat, weighted=weighted)
  } else {
    nam <- nam_or_graph
  }
  nam <- t(nam)
  identity_matrix <- Matrix::Diagonal(n=dim(nam)[1])
  que = identity_matrix - alpha * nam
  que = as.matrix(que)
  d_name <- Matrix::solve(que) * (1 - alpha)
  return(d_name)
}

#  Calculate the corresponding vector of each gene set after network propagation.
#' Title
#'
#' @param individual_heats_matrix A matrix,the result of function get_individual_heats_matrix
#' @param nodes A character vector represents all the genes included in the biological network, and the order should be consistent with the individual_heats_matrix.
#' @param seed_genes A character vector represents the genes included in the gene program
#' @import stats
#' @return A numeric vector,the result of network_propagation
#' @export
#'
#'
network_propagation <- function(individual_heats_matrix, nodes, seed_genes) {
  seed_genes <- intersect(nodes, seed_genes)
  F <- rep(0, length(nodes))
  indices <- match(seed_genes, nodes)
  for (index in indices) {
    F <- F + individual_heats_matrix[, index]
  }
  F <- F / length(seed_genes)
  heat_series <- stats::setNames(F, nodes)
  return(heat_series)
}

#' @Title Input a biological network and a gene set corresponding to gene programs, and output a distance matrix representing the cosine distances between gene programs after network propagation.
#'
#' @param gene_programs A dataframe where the row names are the names of gene programs, with each row corresponding to a gene set of a gene module. Gene programs with fewer genes are filled with NA to fill in the gaps.
#' @param biological_network A symmetric adjacency matrix representing a biological network, with row and column names as gene symbols.
#' @import igraph
#' @return A distance matrix where the row and column names are the names of gene programs, and the values represent the cosine distances between the gene programs.
#' @export
#'
#'
gene_program_network_propagation_dist <- function(gene_programs,biological_network){
  gene_list <- rownames(biological_network)
  G <- igraph::graph_from_adjacency_matrix(biological_network, mode = "undirected", weighted = TRUE)
  adj <- get_normalized_adjacency_matrix(G)
  heat <- get_individual_heats_matrix(adj)
  vec_matrix <- matrix(0, ncol = length(gene_list), nrow = nrow(gene_programs))
  for (i in 1:nrow(gene_programs)) {
    message("Processing row:", i, "\n")
    g1 <- unlist(gene_programs[i, ])
    g1 <- g1[!is.na(g1)]
    vec <- network_propagation(heat,gene_list,g1)
    vec_matrix[i, ] <- vec
  }
  rownames(vec_matrix) = rownames(gene_programs)
  distance_mat <- apply(vec_matrix, 1, function(row) {
    row = as.numeric(row)
    row_distance <- apply(vec_matrix, 1, function(other_row) {
      other_row = as.numeric(other_row)
      cosine_distance(row, other_row)
    })
    return(row_distance)
  })
  return(distance_mat)
}
