#' Calculate the PVE
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param map geno information, SNP, Chr, Pos
#' @param random list, non-additive genetic random effect
#' @param alpha genetype scaling value
#' @param ldscore LD score
#' @param bin region size for correlation calculation
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param ncpus integer value, default as 1, how many threads used
#' @param verbose logical, TRUE for print the output
#'
#' @return matrix, Beta, SE, Pvalue
#' @export
#'
#' @examples
#' #\donttest{
#' #revl <- h2.m(y = y, CV = CV, geno = geno, map = map,
#' #random = random, alpha = alpha, ldscore = ldscore,
#' #bin = bin, EMsteps = EMsteps, EM_alpha = EM_alpha,
#' #EMsteps_fail = EMsteps_fail, max_iter = max_iter,
#' #eps = eps, ncpus = ncpus, verbose = verbose)
#' #}
h2.m <- function(y = NULL, CV = NULL, geno = NULL, map = NULL,
                 random = NULL, alpha = 0, ldscore = NULL, bin = 1e6,
                 EMsteps = 0, EMsteps_fail = 10, EM_alpha = 1, max_iter = 50,
                 eps = 1e-02, ncpus = 1, verbose = TRUE) {

  n <- dim(geno)[2]
  m <- dim(geno)[1]
  if (dim(map)[1] != m) {
    stop("Number of markers don't match between map and geno!")
  }
  if (sum(is.na(map)) != 0) stop("Missing value in map!")

  glm <- ols.fast(y = y, CV = CV, geno = geno,
                  ncpus = ncpus, verbose = verbose)
  z_2 <- (glm[, 1] / glm[, 2]) ^ 2

  if (verbose) {
    cat("Calculating MAF...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  maf <- rowMeans(geno[, ]) / 2
  maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]
  if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
  }
  if (sum(is.na(z_2)) != 0) {
    warning(sum(is.na(z_2)), " NAs in OLS result.")
  }
  if (length(which(maf == 0)) != 0) {
    warning(length(which(maf == 0)), " SNPs are homozygote.")
  }
  z_2[which(is.na(z_2))] <- 0
  maf[which(maf == 0.5)] <- 0.01 / n

  if (is.null(ldscore)) {
    if (is.integer(y)) {
      if (verbose) {
        cat("LD scores not recommended for discrete traits", "\n")
      }
      ldscore <- rep(0, length(maf))
    }else {
      if (is.null(map)) {
        stop("geno map is required for LD score caculating!")
      }
      ldscore <- LDSC(geno = geno, map = map, bin = bin,
                      ncpus = ncpus, verbose = verbose)
    }
  }else {
    if (is.integer(y)) {
      if (verbose) {
        cat("LD scores not recommended for discrete traits", "\n")
      }
      ldscore <- rep(0, length(maf))
    }else {
      ldscore <- as.vector(ldscore)
      if (length(ldscore) != m) {
        stop("Wrong LD score provided!")
      }else {
        if (verbose) {
          cat("LD score has been provided\n")
        }
      }
    }
  }
  vy <- stats::var(y, na.rm = TRUE)
  b2_adj <- (1 - (n / (n + z_2))) / (1 + ldscore)
  if (length(alpha) == 1) {
    if (verbose) {
      cat("PVE calculated with alpha =", alpha)
    }
    pve <- b2_adj * ((2 * maf * (1 - maf)) ^ (1 + alpha))
    pve <- pve / vy
    pve[which(is.na(pve))] <- min(pve, na.rm = TRUE)
  }else {
    mina <- min(alpha, na.rm = TRUE)
    maxa <- max(alpha, na.rm = TRUE)
    if (verbose) {
      t1 <- as.numeric(Sys.time())
      cat("Optimizing alpha at interval(", mina, ",", maxa, ")...\n")
    }
    op_a <- function(al) {
      wei <- (b2_adj * ((2 * maf * (1 - maf)) ^ (1 + al))) / vy
      wei[which(is.na(wei))] <- min(wei, na.rm = TRUE)
      mix <- EBV.trans(y = y, CV = CV, geno = list(geno), weight = wei,
                       random = random, EMsteps = EMsteps, EM_alpha = EM_alpha,
                       EMsteps_fail = EMsteps_fail, max_iter = max_iter,
                       eps = eps, get.KR = FALSE, get.P = FALSE,
                       doTrans = FALSE, verbose = FALSE)
      ll <- mix[[1]]$logL
      if (verbose) {
        cat("  [", "alpha:", al, ";", "LogL:", ll, "]\n")
      }
      rm(mix, wei)
      return(ll)
    }
    opa <- stats::optimize(op_a, c(mina, maxa), maximum = TRUE)$maximum
    if (verbose) {
      t2 <- as.numeric(Sys.time())
      cat(paste("  | 100% elapsed=", time.trans(round(t2 - t1)),sep = ""), "\n")
      cat("PVE calculated with alpha =", opa)
    }
    pve <- b2_adj * ((2 * maf * (1 - maf)) ^ (1 + opa))
    pve <- pve / vy
    pve[which(is.na(pve))] <- min(pve, na.rm = TRUE)
  }
  revl <- cbind(map, glm, maf, ldscore, pve)
  colnames(revl) <- c(colnames(map), colnames(glm), "MAF", "LD_score", "PVE")
  return(revl)
}
