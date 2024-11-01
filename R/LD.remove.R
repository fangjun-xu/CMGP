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
                      LD.threshold = 0.7, ncpus = 1, verbose = TRUE) {
  if (length(index) != length(value)) {
    stop("length of index and value not equal")
  }
  if (sum(is.na(index), is.na(value)) != 0) {
    stop("NAs is not allowed")
  }
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
      index <- as.numeric(names(top.selected))
      rm(sigma, sigma.distance, fit, clusters, top.selected)
      return(index)
    }
    if (length(index) <= 6000) {
      if (verbose) {
        cat("LD removing...\n")
        pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      index_out <- ldr(index = index, value = value)
      if (verbose) {
        pbapply::setTimerProgressBar(pb, 1)
      }
      return(index_out)
    }else {
      mcs <- ifelse(Sys.info()[["sysname"]] == "Linux", ncpus, 1)
      nbreak <- ceiling(length(index) / 3000)
      cutind <- cut(seq_len(length(index)), breaks = nbreak, labels = FALSE)
      com <- utils::combn(1 : nbreak, 2)
      if (ncol(com) > mcs && mcs > 1) {
        pbseq <- seq(ncol(com), 1, -ceiling(ncol(com) / ncpus))[1:ncpus]
      }else {
        pbseq <- seq_len(ncol(com))
      }
      if (verbose) {
        cat("LD removing...\n")
        maxpb <- length(pbseq) + 1
        pb <- pbapply::timerProgressBar(max = maxpb, width = 30,
                                        char = "-", style = 3)
        on.exit(close(pb))
      }
      indmc <- parallel::mclapply(seq_len(ncol(com)), function(i) {
        a <- com[1, i]
        b <- com[2, i]
        ida <- index[which(cutind == a, arr.ind = TRUE)]
        va <- value[which(cutind == a, arr.ind = TRUE)]
        idb <- index[which(cutind == b, arr.ind = TRUE)]
        vb <- value[which(cutind == b, arr.ind = TRUE)]
        idi <- c(ida, idb)
        vi <- c(va, vb)
        outi <- ldr(index = idi, value = vi)
        sign <- rep(NA, length(idi))
        sign[which(idi %in% outi)] <- 1
        if (verbose && i %in% pbseq) {
          pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
        }
        return(rbind(idi, sign))
      }, mc.cores = mcs)
      indmc <- do.call(cbind, indmc)
      if (sum(is.na(indmc[2, ])) == 0) {
        index_out <- unique(indmc[1, ])
      }else {
        idNA <- indmc[1, which(is.na(indmc[2, ]))]
        id1 <- indmc[1, which(!is.na(indmc[2, ]))]
        index_out <- setdiff(id1, idNA)
        index_out <- unique(index_out)
      }
      if (verbose) {
        pbapply::setTimerProgressBar(pb, maxpb)
      }
      return(index_out)
    }
  }
}

