#' Run C-score-var Analysis
#'
#' @param comm Abundance matrix (rows are samples, columns are OTUs).
#' @param nReps Number of replications, default 500.
#' @param save_plot Logical, whether to save the histogram plot.
#' @param output_file Output file path for summary.
#' @param plot_file Output file path for plot (EMF format).
#' @importFrom EcoSimR cooc_null_model
#' @importFrom devEMF emf
#' @importFrom graphics plot
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @export
run_cscore_var <- function(comm, nReps = 500, save_plot = FALSE, output_file = "c-score-var.summary.txt", plot_file = "c-score-var.hist.emf") {
  
  # Set seed for reproducibility as in original script
  set.seed(56)
  
  # Preprocessing
  table01 <- t(comm)
  table01[table01 > 0] <- 1
  # Filter out empty rows
  table01 <- table01[which(rowSums(table01) > 0), ]
  
  # C-score-var calculation
  csvarModel <- EcoSimR::cooc_null_model(
    table01,
    algo = "sim9",
    metric = "c_score_var",
    nReps = nReps,
    saveSeed = FALSE,
    burn_in = 500,
    algoOpts = list(),
    metricOpts = list(),
    suppressProg = FALSE
  )
  
  # Output results
  if (!is.null(output_file)) {
    write.table('C-score-var summary', output_file, append = TRUE)
    sink(output_file, append = TRUE)
    print(summary(csvarModel))
    sink(NULL)
  }
  
  if (save_plot && !is.null(plot_file)) {
    devEMF::emf(file = plot_file, width = 7, height = 7,
                bg = "transparent", fg = "black", pointsize = 12,
                family = "Helvetica", custom.lty = FALSE)
    plot(csvarModel, type = "hist")
    dev.off()
  }
  
  return(csvarModel)
}
