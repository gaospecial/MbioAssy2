#' Construct Co-occurrence Network
#'
#' @param comm Abundance matrix (rows are samples, columns are OTUs).
#' @param cor_cutoff Correlation coefficient cutoff (optional if use_rmt is TRUE).
#' @param p_cutoff P-value cutoff.
#' @param use_rmt Logical, whether to use RMT to find cutoff.
#' @param method Correlation method, default "spearman".
#' @param rmt_interval Interval for RMT threshold search.
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom igraph graph_from_adjacency_matrix simplify V degree write_graph cluster_walktrap modularity membership transitivity mean_distance edge_density diameter ecount vcount betweenness closeness
#' @importFrom RMThreshold rm.get.threshold
#' @export
construct_network <- function(comm, cor_cutoff = 0.8, p_cutoff = 0.05, use_rmt = FALSE, method = "spearman", rmt_interval = c(0.6, 0.99)) {
  
  matrix <- comm # Input is expected to be samples x OTUs. Original script reads table then as.matrix.
  # rcorr expects columns as variables. If comm is samples x OTUs, we need to transpose it because rcorr computes correlations of columns.
  # Original script: rcorr(t(matrix)) where matrix was read from OTU.txt (rows samples, cols OTUs)
  # Wait, original script: table = read.table(...). table is samples x OTUs. 
  # Then rcorr(t(matrix)). t(matrix) is OTUs x samples. 
  # rcorr(x) "x must be a matrix... Computes ... matrix of correlations ... for each pair of columns of x".
  # So if we want correlation between OTUs, we input samples x OTUs (columns are OTUs).
  # The original script does `table <- t(read.table(...))` in some places but `read.table` in others.
  # In 5.1: `table = read.table('example_input/OTU.txt'...)`. OTU.txt has samples as rows? No, usually OTU table has OTUs as rows.
  # Let's check environment details for file content of 5.1_Co-occurrence_network.R
  # In 5.1: `table = read.table('example_input/OTU.txt',sep = '	',header = T, row.names = 1)`.
  # Then `rcorr(t(matrix))`.
  # If OTU.txt has OTUs as rows (standard), then `table` is OTUs x Samples. `t(table)` is Samples x OTUs.
  # `rcorr` on Samples x OTUs calculates correlation between OTUs (columns).
  # So we need to ensure input to rcorr is Samples x OTUs.
  
  # Assuming input `comm` is Samples x OTUs (standard for vegan etc).
  # If comm is Samples x OTUs, rcorr(comm) gives correlations between OTUs.
  # Original script 5.1: table = read.table... (OTU.txt). Then rcorr(t(table)).
  # If OTU.txt is standard (OTU x Sample), then table is OTU x Sample. t(table) is Sample x OTU.
  
  # Let's assume standard input for this function is Samples x OTUs.
  matrix_for_cor <- comm
  
  # Correlation analysis
  matrix.dist <- Hmisc::rcorr(matrix_for_cor, type = method)
  matrix.cor <- matrix.dist$r
  matrix.cor.p <- matrix.dist$P
  
  # RMT thresholding
  final_cor_cutoff <- cor_cutoff
  if (use_rmt) {
    # Using RMT theory to find the optimized correlation threshold
    # Note: rm.get.threshold might be interactive or print plots. We set interactive=F.
    res <- RMThreshold::rm.get.threshold(
      matrix.cor, 
      interactive = FALSE, 
      dist.method = 'LL',
      plot.comp = FALSE, 
      save.fit = FALSE, 
      plot.spacing = FALSE,
      interval = rmt_interval
    )
    
    if (max(res$p.ks > 0.05)) {
      pks.id <- which(res$p.ks >= 0.05)[1]
    } else {
      pks.id <- which.max(res$p.ks)
    }
    
    RMT.cutoff <- res$tested.thresholds[pks.id]
    RMT.cutoff.round <- round(RMT.cutoff, 2)
    if (isTRUE(RMT.cutoff.round < RMT.cutoff)) {
      RMT.cutoff.round <- RMT.cutoff.round + 0.01
    }
    final_cor_cutoff <- RMT.cutoff.round
  }
  
  # FDR correction
  matrix.cor.p <- stats::p.adjust(matrix.cor.p, method = "BH")
  
  # Apply cutoffs
  matrix.cor1 <- matrix.cor
  matrix.cor1.p <- matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= final_cor_cutoff)] <- 0
  matrix.cor1[which(matrix.cor1.p > p_cutoff)] <- 0
  
  # Remove isolated nodes (rows/cols with sum = 0, but diag is 1, so rowSum != 1)
  # Diagonals in correlation matrix are 1.
  # If a node has no edges, its row sum in adjacency matrix (with 0 on diag) would be 0.
  # But matrix.cor1 has 1 on diagonal.
  # Original script: matrix.cor1[which(rowSums(matrix.cor1)!=1),]
  # This implies removing rows where the only non-zero value is the diagonal.
  keep_idx <- which(rowSums(abs(matrix.cor1)) > 1)
  matrix.cor1 <- matrix.cor1[keep_idx, keep_idx]
  
  # Generate graph
  g1 <- igraph::graph_from_adjacency_matrix(matrix.cor1, weight = TRUE, mode = "undirected", diag = FALSE)
  g1 <- igraph::simplify(g1)
  igraph::V(g1)$label <- igraph::V(g1)$name
  igraph::V(g1)$degree <- igraph::degree(g1)
  
  return(list(graph = g1, cor_matrix = matrix.cor1, cutoff = final_cor_cutoff))
}

#' Calculate Network Topology
#'
#' @param g igraph object.
#' @param output_dir Directory to save topology files (optional).
#' @param prefix Prefix for output files.
#' @importFrom igraph cluster_walktrap modularity membership transitivity mean_distance edge_density diameter degree ecount vcount betweenness closeness
#' @importFrom utils write.csv
#' @export
calculate_topology <- function(g, output_dir = NULL, prefix = "Network") {
  
  # Global topology
  c <- igraph::cluster_walktrap(g)
  md <- igraph::modularity(g, igraph::membership(c), weights = NULL)
  cc <- igraph::transitivity(g, vids = NULL, weights = NULL)
  spl <- igraph::mean_distance(g, directed = FALSE, unconnected = TRUE)
  gd <- igraph::edge_density(g, loops = FALSE)
  nd <- igraph::diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)
  node.degree <- igraph::degree(g, v = igraph::V(g), mode = "all")
  ad <- mean(node.degree)
  e <- igraph::ecount(g)
  v <- igraph::vcount(g)
  
  global.topology <- data.frame(e, v, cc, spl, md, gd, nd, ad)
  
  # Node topology
  betweenness.centrality <- igraph::betweenness(g, v = igraph::V(g), directed = FALSE, weights = NA, normalized = FALSE)
  closeness.centrality <- igraph::closeness(g, vids = igraph::V(g), weights = NA, normalized = FALSE)
  node.transitivity <- igraph::transitivity(g, type = "local", vids = NULL, weights = NA)
  
  node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write.csv(global.topology, file = file.path(output_dir, paste0(prefix, ".global.topology.csv")), row.names = FALSE)
    write.csv(node.topology, file = file.path(output_dir, paste0(prefix, ".node.topology.csv")), row.names = TRUE)
  }
  
  return(list(global = global.topology, node = node.topology))
}
