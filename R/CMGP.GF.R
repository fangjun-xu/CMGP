#' Genomin factor BLUP
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param map geno information, SNP, Chr, Pos
#' @param random list, non-additive genetic random effect
#' @param alpha genetype scaling value
#' @param ldscore LD score
#' @param bin region size for correlation calculation
#' @param ncpus integer value, default as 1, how many threads used
#' @param part.method prioritization method
#' @param part.threshold significant level
#' @param LD.threshold cluster threshold
#' @param doTrans logical, TRUE for print the output
#' @param get.KR logical, should the Kinship be returned
#' @param get.P see ?gaston::lmm.aireml
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param verbose logical, TRUE for print the output
#'
#' @return list
#' @export
#'
#' @examples
#' #\donttest{
#' #revl <- CMGP.GF(y = y, CV = CV, geno = geno, map = map, random = random,
#' #alpha = alpha, ldscore = ldscore, bin = bin, ncpus = ncpus,
#' #part.method = part.method, part.threshold = part.threshold,
#' #LD.threshold = LD.threshold, doTrans = doTrans,
#' #get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
#' #EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
#' #max_iter = max_iter, eps = eps, verbose = verbose)
#' #}
CMGP.GF <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
                    alpha = c(-3, 1), ldscore = NULL, bin = 1e6, ncpus = 1,
                    part.method = "ebvTrans", part.threshold = 0.01,
                    LD.threshold = 0.7, doTrans = FALSE, get.KR = FALSE,
                    get.P = FALSE, EMsteps = 0L, EMsteps_fail = 10L,
                    EM_alpha = 1, max_iter = 50L, eps = 1e-02, verbose = TRUE) {

  #---------------------------------------------------------#
  #GFBLUP (Genomic Factor) for Genomic Prediction
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "GFBLUP  Running",
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
    index <- which(sp$PVE <= stats::quantile(sp$PVE,
                                             part.threshold,
                                             na.rm = TRUE))
    value <- sp$PVE[index]
  }
  if (verbose) {
    cat(length(index), "signals were preliminarily screened out\n")
  }
  if (length(index) >= 2) {
    index <- LD.remove(index = index, geno = geno, value = value,
                       LD.threshold = LD.threshold, verbose = verbose)
    if (verbose) {
      cat(length(index), "signals retained as genomic factor\n")
    }
  }
  gfg <- matrix(geno[index, ], nrow = length(index))
  rownames(gfg) <- rownames(geno)[index]
  colnames(gfg) <- colnames(geno)
  gf.geno <- list(gfg, geno[-index, ])
  rm(gfg)
  gf.fit <- EBV.trans(y = y, CV = CV, geno = gf.geno, random = random,
                      ncpus = ncpus, max_iter = max_iter, eps = eps,
                      EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                      EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                      doTrans = doTrans, verbose = verbose)
  if (doTrans) {
    gf.fit$GWAS <- gf.fit$GWAS[rownames(geno), ]
  }
  if (verbose) {
    cat(paste(paste(rep("-", 23), collapse = ""),
              "GFBLUP FINISH",
              paste(rep("-", 23), collapse = ""), sep = ""), "\n")
  }
  revl <- list(GP = gf.fit$GP, GWAS = gf.fit$GWAS, GF.index = index)
  return(revl)
}
