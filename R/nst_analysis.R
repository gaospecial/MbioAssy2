#' Run Normalized Stochasticity Ratio (NST) Analysis
#'
#' @param comm Community abundance matrix (rows are samples, columns are OTUs).
#' @param group Group information (data frame with row names matching samples).
#' @param dist_method Distance method, default "jaccard".
#' @param abundance_weighted Logical, default TRUE.
#' @param rand Number of randomizations, default 20.
#' @param output_file Output file path for the summary.
#' @importFrom NST tNST
#' @importFrom utils write.table
#' @export
run_nst <- function(comm, group, dist_method = "jaccard", abundance_weighted = TRUE, rand = 20, output_file = "NST.output.txt") {
  
  # Preprocessing
  comm <- as.matrix(comm)
  comm <- comm[which(rowSums(comm) > 0), ]
  comm <- comm[, which(colSums(comm) > 0)]
  
  # NST calculation
  nst_res <- NST::tNST(
    comm = comm,
    group = group,
    dist.method = dist_method,
    abundance.weighted = abundance_weighted,
    rand = rand,
    null.model = "PF"
  )
  
  nst_sum <- nst_res$index.grp
  
  if (!is.null(output_file)) {
    write.table(nst_sum, output_file, sep = "	", quote = FALSE, row.names = FALSE)
  }
  
  return(nst_sum)
}
