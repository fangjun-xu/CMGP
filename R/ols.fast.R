#' Parallel OLS excluding the additional covariates
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param ncpus integer value, default as 1, how many threads used
#' @param verbose logical, TRUE for print the output
#'
#' @return matrix Beta, SE, Pvalue
#' @export
#'
#' @examples
#' #\donttest{
#' #revl <- ols.fast(y = y, CV = CV, geno = geno,
#' #ncpus = ncpus, verbose = verbose)
#' #}
ols.fast <- function(y = NULL, CV = NULL, geno = NULL,
                     ncpus = 1, verbose = TRUE) {
  #---------------------------------------------------------#
  #GLM model
  if (!is.null(CV) && length(y) != dim(CV)[1]) {
    stop("Number of individuals don't match between phenotype and covariates!")
  }
  if (!bigmemory::is.big.matrix(geno) && !is.matrix(geno)) {
    stop("Wrong format of geno inputed! matrix/big.matrix is allowed.")
  }
  if (sum(is.na(geno[, ])) != 0) stop("'NA' isn't allowed in geno!")
  if (length(y) != dim(geno)[2]) {
    stop("Number of individuals don't match between y and geno!", "\n",
         "In addition: geno should be (marker size) * (individual size)!")
  }
  if (dim(geno)[1] > ncpus && ncpus > 1) {
    pbseq <- seq(dim(geno)[1], 1, -ceiling(dim(geno)[1] / ncpus))[1:ncpus]
  }else {
    pbseq <- 1:dim(geno)[1]
  }
  if (verbose) {
    cat("Calculating marginal effects...\n")
    maxpb <- length(pbseq) + 1
    pb <- pbapply::timerProgressBar(max = maxpb, width = 30,
                                    char = "-", style = 3)
    on.exit(close(pb))
  }
  ret <- which(!is.na(y))
  y <- y[ret]
  CV <- CV[ret, ]
  if (is.matrix(geno)) {
    geno <- geno[, ret]
  }else {
    geno <- bigmemory::deepcopy(geno, col = ret)
  }
  y <- matrix(y, ncol = 1)
  w <- CovMatrix(CV = CV, n = length(y))
  DF <- length(y) - ncol(w) - 1
  wwi <- try(solve(crossprod(w)), silent = TRUE)
  if (inherits( wwi, "try-error")) {
    warning("The cor-matrix of CV is singular!")
    wwi <- MASS::ginv(crossprod(w))
  }
  My <- y - tcrossprod(tcrossprod(w, wwi), w) %*% y
  yy <- sum(My ^ 2)
  

  OL <- function(i) {
    x <- matrix(geno[i, ], ncol = 1)
	xx <- sum(x ^ 2)
	xw <- crossprod(x, w)
	xy <- crossprod(x, My)
	xhxh <- xx - tcrossprod(tcrossprod(xw, wwi), xw)
	
	b <- xy / xhxh
	RSS <- yy - b * xy	
	se <- sqrt(RSS / (xhxh * DF))
    
    pvalue <- 2 * stats::pt(abs(b / se), DF, lower.tail = FALSE)
    pvalue[which(is.na(pvalue))] <- 1
    revl <- c(b, se, pvalue)
    if (verbose && i %in% pbseq) {
      pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
    }
    return(revl)
  }
  if (dim(geno)[1] == 1 || ncpus == 1) {
    revl <- lapply(1:dim(geno)[1], FUN = OL)
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      revl <- parallel::mclapply(1:dim(geno)[1], FUN = OL,
                                 mc.cores = ncpus,
                                 mc.preschedule = TRUE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      suppressMessages(snowfall::sfExport("geno", "My", "w", "yy", "DF",
                                          "wwi","pbseq", "pb", "verbose"))
      suppressMessages(snowfall::sfLibrary(pbapply))
      suppressMessages(snowfall::sfLibrary(stats))
      revl <- snowfall::sfLapply(1:dim(geno)[1], fun = OL)
      suppressMessages(snowfall::sfStop())
    }
  }
  revl <- do.call(rbind, revl)
  colnames(revl) <- c("Effect", "SE", "Pvalue")
  if (verbose) {
    Sys.sleep(0.1)
    pbapply::setTimerProgressBar(pb, maxpb)
  }
  return(revl)
}
