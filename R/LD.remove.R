#' Select the most important marker in each cluster
#'
#' @param index index of markers
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param value improtance
#' @param LD.threshold cluster threshold
#' @param verbose logical, TRUE for print the output
#'
#' @return index
#' @export
#'
#' @examples
#' #\donttest{
#' #index <- LD.remove(index = index, geno = geno, value = value,
#' #LD.threshold = LD.threshold, verbose = verbose)
#' #}
LD.remove <- function(index = NULL, geno = NULL, value = NULL,
                      LD.threshold = 0.7, verbose = TRUE) {
  if (length(index) < 2) {
    return(index)
  }else {
    ldr <- function(index, value) {
      X1 <- geno[index, ]
      X1<-sweep(X1[, ], 1, rowMeans(X1))
      rownames(X1) <- index
      sigma <- tcrossprod(X1) / (dim(geno)[2] - 1)
      sigma.distance <- stats::as.dist(1 - abs(stats::cov2cor(sigma)))
      fit <- stats::hclust(sigma.distance, method = "single")
      clusters <- stats::cutree(fit, h = 1 - LD.threshold)
      names(value) <- rownames(X1)
      top.selected <- sapply(1:max(clusters), function(i) {
	cluster_elements <- clusters == i
	top_within <- which.min(value[cluster_elements])
	if (length(top_within) == 0) top_within <- 1
	return(which(cluster_elements)[top_within])
      })
      index<-as.numeric(names(top.selected))
      rm(sigma, sigma.distance, fit, clusters, top.selected)
      return(index)
    }
    if (length(index) <= 10000) {
      if (verbose) {
        cat("LD removing...\n")
        pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      index <- ldr(index = index, value = value)
      if (verbose) {
        pbapply::setTimerProgressBar(pb, 1)
      }
    }else {
      nbreak <- ceiling(length(index) / 10000)
      cutind <- cut(1:length(index), breaks = nbreak, labels = FALSE)
      if (verbose) {
        cat("LD removing...\n")
        pb <- pbapply::timerProgressBar(max = nbreak + 1, width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      indi <- lapply(1:nbreak, function(i) {
        idi <- index[which(cutind == i, arr.ind = TRUE)]
        vi <- value[which(cutind == i, arr.ind = TRUE)]
	if (verbose) {
          pbapply::setTimerProgressBar(pb, i)
        }
        return(ldr(index = idi, value = vi))
      })
      indi <- unlist(indi)
      vali <- value[which(index %in% indi)]
      index <- ldr(index = indi, value = vali)
      if (verbose) {
        pbapply::setTimerProgressBar(pb, nbreak + 1)
      }
    }	
    return(index)
  }
}

