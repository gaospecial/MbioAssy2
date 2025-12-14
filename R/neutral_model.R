#' Fit Neutral Community Model
#'
#' @param spp Abundance matrix (rows are samples, columns are OTUs).
#' @param pool Metacommunity abundance (optional).
#' @param taxon Taxonomy information (optional).
#' @importFrom minpack.lm nlsLM
#' @importFrom Hmisc binconf
#' @importFrom stats4 mle
#' @importFrom stats coef confint dnorm pbeta ppois AIC BIC
#' @importFrom ggplot2 ggplot geom_point scale_fill_manual scale_fill_discrete geom_line xlab ylab ggtitle annotate theme_bw theme element_blank element_line aes
#' @importFrom rlang .data
#' @export
fit_sncm <- function(spp, pool = NULL, taxon = NULL) {
  options(warn = -1)
  
  # Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  # Calculate the average relative abundance of each taxa across communities
  if (is.null(pool)) {
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  }
  
  # Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1 * (spp > 0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # Combine
  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]), ]
  C <- as.data.frame(C)
  # Removes rows with any zero (absent in either source pool or local communities)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  # Calculate the limit of detection
  d <- 1 / N
  
  # Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- tryCatch({
    minpack.lm::nlsLM(freq ~ stats::pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.001))
  }, error = function(e) {
    message("Error in nlsLM: ", e$message)
    return(NULL)
  })
  
  if (is.null(m.fit)) return(NULL)
  
  m.ci <- confint(m.fit, 'm', level = 0.95)
  
  # Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- stats::pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq - freq.pred)^2) / (length(freq) - 1))
  
  pred.ci <- Hmisc::binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
  
  # Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma) {
    R <- freq - stats::ppois(d, N * p, lower.tail = FALSE)
    R <- stats::dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- stats4::mle(pois.LL, start = list(mu = 0, sigma = 0.1), nobs = length(p))
  
  aic.pois <- stats::AIC(pois.mle, k = 2)
  bic.pois <- stats::BIC(pois.mle)
  
  # Goodness of fit for Poisson model
  pois.pred <- stats::ppois(d, N * p, lower.tail = FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2)) / (sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2) / (length(freq) - 1))
  
  # Results
  fitstats <- data.frame(
    m = as.numeric(coef(m.fit)),
    m.ci = as.numeric(coef(m.fit) - m.ci[1]),
    poisLL = as.numeric(pois.mle@details$value),
    Rsqr = as.numeric(Rsqr),
    Rsqr.pois = as.numeric(Rsqr.pois),
    RMSE = as.numeric(RMSE),
    RMSE.pois = as.numeric(RMSE.pois),
    AIC.pois = as.numeric(aic.pois),
    BIC.pois = as.numeric(bic.pois),
    N = as.numeric(N),
    Samples = as.numeric(nrow(spp)),
    Richness = as.numeric(length(p)),
    Detect = as.numeric(d)
  )
  
  A <- cbind(p, freq, freq.pred, pred.ci[, 2:3])
  A <- as.data.frame(A)
  colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr')
  if (is.null(taxon)) {
    B <- A[order(A[, 1]), ]
  } else {
    B <- merge(A, taxon, by = 0, all = TRUE)
    row.names(B) <- B[, 1]
    B <- B[, -1]
    B <- B[order(B[, 1]), ]
  }
  B <- B[!is.na(B$freq), ]
  # fit_class for graphing
  B$fit_class <- "As predicted"
  B[which(B$freq < B$pred.lwr), "fit_class"] <- "Below prediction"
  B[which(B$freq > B$pred.upr), "fit_class"] <- "Above prediction"
  B[which(is.na(B$freq)), "fit_class"] <- "NA"
  
  # combine fit stats and predictions into list
  i <- list(fitstats = fitstats, predictions = B)
  return(i)
}

#' Plot Neutral Community Model Fit
#'
#' @param spp_out Output object from fit_sncm.
#' @param fill Variable to fill points with, default "fit_class".
#' @param title Plot title.
#' @import ggplot2
#' @importFrom rlang .data
#' @export
plot_sncm_fit <- function(spp_out, fill = NULL, title = NULL) {
  
  # Extract taxonomy levels if available
  # Assuming predictions has columns: p, freq, freq.pred, pred.lwr, pred.upr, fit_class, and optionally taxonomy
  # fit_sncm result 'predictions' is a data frame.
  
  if (is.null(fill)) {
    fill <- "fit_class"
  }
  
  # Bind global variables to avoid R CMD check notes
  p <- freq <- freq.pred <- pred.lwr <- pred.upr <- NULL
  
  # Only use tax_levels if there are extra columns
  cols <- colnames(spp_out$predictions)
  # Standard columns are p, freq, freq.pred, pred.lwr, pred.upr, fit_class (6 columns)
  # Anything else is likely taxonomy
  
  r2_val <- paste("r^2 ==", round(spp_out$fitstats$Rsqr, 4))
  m_val <- paste("m ==", round(spp_out$fitstats$m, 4))
  
  # Create summary for legend labels
  df <- data.frame(t(table(spp_out$predictions$fit_class)))
  # Check if table has results
  if (ncol(df) >= 3) {
      df <- df[, c(2, 3)]
      colnames(df) <- c("Prediction", "AVS Abundance")
  } else {
      # Handle case where not all classes are present or structure is different
      df <- data.frame(Prediction = names(table(spp_out$predictions$fit_class)), 
                       "AVS Abundance" = as.vector(table(spp_out$predictions$fit_class)))
  }

  
  plt <- ggplot(data = spp_out$predictions)
  
  # Helper to safe-get count
  get_count <- function(pred_class) {
      val <- df[df$Prediction == pred_class, 2]
      if (length(val) == 0) return(0)
      return(val)
  }

  if (fill == "fit_class") {
    plt <- plt + geom_point(aes(x = log(p), y = freq, fill = .data[[fill]]), shape = 21, color = "black", size = 2, alpha = 0.75)
    
    # Calculate percentages
    richness <- spp_out$fitstats$Richness
    
    labels_vec <- c(
        paste0("Above prediction (", round((get_count("Above prediction") / richness) * 100, 1), "%)"),
        paste0("As predicted (", round((get_count("As predicted") / richness) * 100, 1), "%)"),
        paste0("Below Prediction (", round((get_count("Below prediction") / richness) * 100, 1), "%)"),
        paste0("NA (", get_count("NA"), ")")
    )
    names(labels_vec) <- c("Above prediction", "As predicted", "Below prediction", "NA")
    
    plt <- plt + scale_fill_manual(
      name = "Prediction",
      values = c("Above prediction" = "seagreen", "As predicted" = "black", "Below prediction" = "tan1", "NA" = "white"),
      breaks = c("Above prediction", "As predicted", "Below prediction", "NA"),
      labels = labels_vec
    )
    
  } else {
    # Check if fill is in columns
    if (fill %in% colnames(spp_out$predictions)) {
        plt <- plt + geom_point(aes(x = log(p), y = freq, fill = .data[[fill]]), shape = 21, color = "black", size = 2, alpha = 0.75)
        plt <- plt + scale_fill_discrete(name = fill)
    } else {
        warning(paste0("fill variable: ", fill, " is not a valid column in predictions"))
        # Fallback
        plt <- plt + geom_point(aes(x = log(p), y = freq), shape = 21, color = "black", fill="grey", size = 2, alpha = 0.75)
    }
  }
  
  plt <- plt + geom_line(aes(x = log(p), y = freq.pred), color = "dodgerblue4", lwd = 1.5)
  plt <- plt + geom_line(aes(x = log(p), y = pred.lwr), color = "dodgerblue4", linetype = "dashed", lwd = 1.5)
  plt <- plt + geom_line(aes(x = log(p), y = pred.upr), color = "dodgerblue4", linetype = "dashed", lwd = 1.5)
  plt <- plt + xlab("log(Mean Relative Abundance)")
  plt <- plt + ylab("Frequency")
  if (!is.null(title)) plt <- plt + ggtitle(title)
  plt <- plt + annotate("text", x = -5, y = 0.65, size = 5, label = r2_val, parse = TRUE)
  plt <- plt + annotate("text", x = -5, y = 0.5, size = 5, label = m_val, parse = TRUE)
  plt <- plt + theme_bw()
  plt <- plt + theme(panel.grid = element_blank(), element_line(size = 1, colour = "black"))
  return(plt)
}
