#' Marker Assisted BLUP
#'
#' @param y phenotype, n*1 vector, NAs refer to validation
#' @param CV covariates, n*q matrix, numeric format for quantitative, character or factor for discrete
#' @param geno genotype, m * n matrix/bigmatrix
#' @param map marker information, SNP, Chr, Pos
#' @param random list, non-additive genetic effects
#' @param alpha a numeric value for genotype scaling
#' @param ldscore vector, the sum of adjust squared correlation between a marker and markers within the region(in bin size)
#' @param bin region size for correlation calculation
#' @param ncpus integer value for number of threads used
#' @param part.method c("OLS", "ebvTrans", "h2"), the first is default
#' @param part.threshold significant level for statistic test
#' @param LD.threshold threshold for clustering when LD removing is working
#' @param doTrans logical, FALSE is default, if TRUE, gwas by linear transform will work
#' @param get.KR logical, should the Kinship be output
#' @param get.P see ?gaston::lmm.aireml
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param verbose logical, should the running progress be printed
#'
#' @return list
#' @export
#'
#' @examples
#' #\donttest{
#' #revl <- CMGP.MA(y = y, CV = CV, geno = geno, map = map, random = random,
#' #alpha = alpha, ldscore = ldscore, bin = bin, ncpus = ncpus,
#' #part.method = part.method, part.threshold = part.threshold,
#' #LD.threshold = LD.threshold, doTrans = doTrans,
#' #get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
#' #EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
#' #max_iter = max_iter, eps = eps, verbose = verbose)
#' #}
CMGP.MA <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
                    alpha = c(-3, 1), ldscore = NULL, bin = 1e6, ncpus = 1,
                    part.method = "ebvTrans", part.threshold = 0.01,
                    LD.threshold = 0.7, doTrans = FALSE, get.KR = FALSE,
                    get.P = FALSE, EMsteps = 0L, EMsteps_fail = 10L,
                    EM_alpha = 1, max_iter = 50L, eps = 1e-02, verbose = TRUE) {

  #---------------------------------------------------------#
  #MAGBLUP (Marker-assisted) for Genomic Prediction
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "MAGBLUP Running",
              paste(rep("-", 22), collapse = ""), sep = ""), "\n")
  }

  sp <- SNP.prioritized(y = y, CV = CV, geno = geno, map = map,
                        random = random, alpha = alpha, ldscore = ldscore,
                        bin = bin, method = part.method, EMsteps = EMsteps,
                        EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
                        max_iter = max_iter, eps = eps,
                        ncpus = ncpus, verbose = verbose)
  sp <- as.data.frame(sp)
  if (part.method == "OLS" || part.method == "ebvTrans") {
    if (verbose) {
      cat("Significant level:", part.threshold / nrow(sp), "\n")
    }
    index <- which(sp$Pvalue <= (part.threshold / nrow(map)))
    if (length(index) == 0) {
      if (verbose) {
        cat("No signals under Bonferroni correction\nSelected the top markers\n")
      }
      index <- which(sp$Pvalue <= stats::quantile(sp$Pvalue,
                                                  part.threshold,
                                                  na.rm = TRUE))
    }
    value <- sp$Pvalue[index]
  }else {
    index <- which(sp$PVE >= stats::quantile(sp$PVE,
                                             1 - part.threshold,
                                             na.rm = TRUE))
    value <- sp$PVE[index]
  }
  if (verbose) {
    cat(length(index), "signals were preliminarily screened out\n")
  }
  if (length(index) >= 2) {
    index <- LD.remove(index = index, geno = geno, value = value,
                       LD.threshold = LD.threshold, verbose = verbose)
    if (length(index) > 300) {
      if (verbose) {
        cat("The remaining signals were plethora, select the top 300\n")
      }
      if (part.method == "OLS" || part.method == "ebvTrans") {

        index300 <- which(sp$Pvalue[index] <= stats::quantile(sp$Pvalue[index],
                                                              300 / length(index),
                                                              na.rm = TRUE))
        index <- index[index300]
      }else {
        index300 <- which(sp$PVE[index] >= stats::quantile(sp$PVE[index],
                                                           (length(index) - 300) / length(index),
                                                           na.rm = TRUE))
        index <- index[index300]
      }
    }

    if (verbose) {
      cat(length(index), "signals were selected as fixed effect\n")
    }
  }
  if (length(index) == 1) {
    ma.cv <- cbind(geno[index, ], CV)
  }else {
    ma.cv <- cbind(t(geno[index, ]), CV)
  }

  ma.geno <- list(geno[-index, ])

  ma.fit <- EBV.trans(y = y, CV = ma.cv, geno = ma.geno, random = random,
                      ncpus = ncpus, max_iter = max_iter, eps = eps,
                      EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                      EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                      doTrans = doTrans, verbose = verbose)
  if (doTrans) {
    gwas <- ma.fit$GWAS
    beta <- ma.fit$GP$beta
    se <- sqrt(diag(ma.fit$GP$varbeta))
    DF <- sum(!is.na(y)) - ncol(ma.cv) - 1
    pv <- 2*stats::pt(abs(beta / se), DF, lower.tail = FALSE)
    magwas <- cbind(beta, se, pv)
    magwas <- magwas[2:(length(index) + 1), ]
    colnames(magwas) <- c("Effect", "SE", "Pvalue")
    row.names(magwas) <- rownames(geno)[index]
    gwas <- rbind(gwas, magwas)
    ma.fit$GWAS <- gwas[rownames(geno), ]
  }
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "MAGBLUP  FINISH",
              paste(rep("-", 22), collapse = ""), sep = ""), "\n")
  }
  revl <- list(GP = ma.fit$GP, GWAS = ma.fit$GWAS, MA.index = index)
  return(revl)
}
