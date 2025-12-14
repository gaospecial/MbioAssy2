#' Generate Random Networks
#'
#' @param n Number of nodes.
#' @param e Number of edges.
#' @param reps Number of random networks to generate, default 1000.
#' @param output_file Output file path for summary (optional).
#' @importFrom igraph erdos.renyi.game cluster_walktrap modularity transitivity mean_distance edge_density diameter degree ecount vcount
#' @importFrom utils write.table
#' @export
generate_random_networks <- function(n, e, reps = 1000, output_file = NULL) {
  
  results <- list()
  
  for (i in 1:reps) {
    g <- igraph::erdos.renyi.game(n, e, 'gnm')
    
    # Global topological features
    c <- igraph::cluster_walktrap(g)
    md <- igraph::modularity(g, igraph::membership(c), weights = NULL)
    cc <- igraph::transitivity(g, vids = NULL, weights = NULL)
    spl <- igraph::mean_distance(g, directed = FALSE, unconnected = TRUE)
    gd <- igraph::edge_density(g, loops = FALSE)
    nd <- igraph::diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
    
    ND <- igraph::degree(g, v = igraph::V(g), mode = "all")
    ad <- mean(ND)
    
    global.topol <- data.frame(n, e, cc, spl, md, gd, nd, ad)
    results[[i]] <- global.topol
    
    if (!is.null(output_file)) {
      write.table(global.topol, file = output_file, 
                  append = TRUE, sep = "	", row.names = FALSE, 
                  col.names = (i == 1))
    }
  }
  
  final_res <- do.call(rbind, results)
  return(final_res)
}
