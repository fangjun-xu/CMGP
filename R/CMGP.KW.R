#' Kinship Weighted BLUP
#'
#' @param y phenotype, n*1 vector, NAs refer to validation
#' @param CV covariates, n*q matrix, numeric format for quantitative, character or factor for discrete
#' @param geno genotype, m * n matrix/bigmatrix
#' @param map marker information, SNP, Chr, Pos
#' @param random list, non-additive genetic effects
#' @param weight vector, weight of markers, using marker names to name it
#' @param alpha a numeric value for genotype scaling
#' @param ldscore vector, the sum of adjust squared correlation between a marker and markers within the region(in bin size)
#' @param bin region size for correlation calculation
#' @param ncpus integer value for number of threads used
#' @param W.method c("OLS", "ebvTrans", "h2"), the first is default
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
#' #revl <- CMGP.KW(y = y, CV = CV, geno = geno, map = map, random = random,
#' #alpha = alpha, ldscore = ldscore, bin = bin,
#' #weight = weight, W.method = part.method,
#' #doTrans = doTrans, get.KR = get.KR, get.P = get.P,
#' #EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
#' #EMsteps = EMsteps, max_iter = max_iter, eps = eps,
#' #ncpus = ncpus, verbose = verbose)
#' #}
CMGP.KW <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
                    weight = NULL, alpha = c(-3, 1), ldscore = NULL, bin = 1e6,
                    ncpus = 1, W.method = "h2", doTrans = FALSE,
                    get.KR = FALSE, get.P = FALSE, EMsteps = 0L,
                    EMsteps_fail = 10L, EM_alpha = 1,
                    max_iter = 50L, eps = 1e-02, verbose = TRUE) {

  #---------------------------------------------------------#
  #KWBLUP (Kinship-weighted) for Genomic Prediction
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "KWGBLUP Running",
              paste(rep("-", 22), collapse = ""), sep = ""), "\n")
  }
  if (!is.null(weight)) {
    if (length(weight) != dim(geno)[1]) {
      stop("Provided weight does not match the geno")
    }
    weight <- as.vector(weight)
  }else {
    sp <- SNP.prioritized(y = y, CV = CV, geno = geno, map = map,
                          random = random, alpha = alpha, ldscore = ldscore,
                          bin = bin, method = W.method, EMsteps = EMsteps,
                          EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
                          max_iter = max_iter, eps = eps,
                          ncpus = ncpus, verbose = verbose)
    sp <- as.data.frame(sp)
    if (W.method == "OLS" || W.method == "ebvTrans") {
      weight <- as.vector(sp$Pvalue)
      weight <- -log10(weight)
    }else {
      weight <- as.vector(sp$PVE)
    }
  }
  weight[which(is.na(weight))] <- min(weight, na.rm = TRUE)
  names(weight) <- rownames(geno)

  kw.fit <- EBV.trans(y = y, CV = CV, geno = list(geno), random = random,
                      weight = weight, ncpus = ncpus, max_iter = max_iter,
                      EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                      EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                      eps = eps, doTrans = doTrans, verbose = verbose)
  if (doTrans) {
    kw.fit$GWAS <- kw.fit$GWAS[rownames(geno), ]
  }
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "KWGBLUP  FINISH",
              paste(rep("-", 22), collapse = ""), sep = ""), "\n")
  }
  revl <- list(GP = kw.fit$GP, GWAS = kw.fit$GWAS, KW.weight = weight)
  return(revl)
}
