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
    if (length(index) <= 5000) {
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
      nbreak <- ceiling(length(index) / 2500)
      cutind <- cut(1:length(index), breaks = nbreak, labels = FALSE)
      com <- 
    }

    if (verbose) {
      cat("LD removing...\n")
      t1 <- as.numeric(Sys.time())
      if (length(index) <= 5000) {
        pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
    }
    
    indi <- index
    vali <- value
    cutn <- 1
    do.repeat <- ifelse(length(indi) > 5000, TRUE, FALSE)
    while(do.repeat) {
      nbreak <- ceiling(length(indi) / 5000)
      cutind <- cut(1:length(indi), breaks = nbreak, labels = FALSE)
      if (nbreak > ncpus && ncpus > 1) {
        pbseq <- seq(nbreak, 1, -ceiling(nbreak / ncpus))[1:ncpus]
      }else {
        pbseq <- 1:nbreak
      }
      if (verbose) {
        cat("Repeat:", cutn, ",", "Num.:", length(indi), "\n")
        maxpb <- length(pbseq) + 1
        pb <- pbapply::timerProgressBar(max = maxpb, width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      
      mc <- ifelse(Sys.info()[["sysname"]] == "Linux", ncpus, 1)
      indmc <- parallel::mclapply(1:nbreak, function(i) {
        idi <- indi[which(cutind == i, arr.ind = TRUE)]
        vi <- vali[which(cutind == i, arr.ind = TRUE)]
	      if (verbose && i %in% pbseq) {
          pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
        }
        return(ldr(index = idi, value = vi))
      }, mc.cores = mc)
      indmc <- unlist(indmc)
      Nin <- length(indi)
      Nout <- length(indmc)
      
      indi <- index[which(index %in% indmc)]
      vali <- value[which(index %in% indmc)]
      cutn <- cutn + 1
      if (verbose) {
        pbapply::setTimerProgressBar(pb, maxpb)
        cat("\n")
      }
      do.repeat <- ifelse((Nout > 5000 && Nout != Nin), TRUE, FALSE)
    }

    index_out <- ldr(index = indi, value = vali)

    if (verbose) {
     if (length(index) <= 5000) {
        pbapply::setTimerProgressBar(pb, 1)
      }else {
        t2 <- as.numeric(Sys.time())
        tud <- time.trans(round(t2 - t1))
        cat(paste("  Finished, elapsed=", tud, sep = ""))
      }
    }
    return(index_out)
  }
}


te=LD.remove(index = 1:70000, geno = geno, value = 1:70000,
                      LD.threshold = 0.7, ncpus = 10, verbose = TRUE)