#' Prioritization of marker using OLS, ebvTrans or PVE
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param map geno information, SNP, Chr, Pos
#' @param random list, non-additive genetic random effect
#' @param alpha genetype scaling value
#' @param ldscore LD score
#' @param bin region size for correlation calculation
#' @param method c("OLS", "ebvTrans", "h2"), "OLS" was default
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param ncpus integer value, default as 1, how many threads used
#' @param verbose logical, TRUE for print the output
#'
#' @return matrix
#' @export
#'
#' @examples
#' #\donttest{
#' #sp <- SNP.prioritized(y = y, CV = CV, geno = geno, map = map,
#' #random = random, alpha = alpha, ldscore = ldscore,
#' #bin = bin, method = "h2", EMsteps = EMsteps,
#' #EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
#' #max_iter = max_iter, eps = eps,
#' #ncpus = ncpus, verbose = verbose)
#' #}
SNP.prioritized <- function(y = NULL, CV = NULL, geno = NULL, map = NULL,
                            random = NULL, alpha = 0, ldscore = NULL, bin = 1e6,
                            method = c("OLS", "ebvTrans", "h2"), EMsteps = 0L,
                            EMsteps_fail = 10L, EM_alpha = 1, max_iter = 50L,
                            eps = 1e-02, ncpus = 1, verbose = TRUE) {

  #---------------------------------------------------------#
  #Prioritization of SNP

  method <- method[1]
  if (verbose) {
    cat(paste(">>>", paste(rep("-", 3), collapse = ""),
              "Prioritization  Start",
              paste(rep("-", 29), collapse = ""), "<<<", sep = ""), "\n")
  }

  if (method == "OLS") {
    if (verbose) {
      cat("Prioritized method: GLM\n")
    }
    revl <- ols.fast(y = y, CV = CV, geno = geno,
                     ncpus = ncpus, verbose = verbose)
  }

  if (method == "ebvTrans") {
    if (verbose) {
      cat("Prioritized method: EBV transform\n")
    }
    revl <- EBV.trans(y = y, CV = CV, geno = list(geno), random = random,
                      ncpus = 1, max_iter = max_iter, eps = eps,
                      EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                      EMsteps = EMsteps, get.KR = FALSE, get.P = FALSE,
                      doTrans = TRUE, verbose = verbose)[[2]]
  }

  if (method == "h2") {
    if (verbose) {
      cat("Prioritized method: Heritability\n")
    }
    revl <- h2.m(y = y, CV = CV, geno = geno, map = map,
                 random = random, alpha = alpha, ldscore = ldscore,
                 bin = bin, EMsteps = EMsteps, EM_alpha = EM_alpha,
                 EMsteps_fail = EMsteps_fail, max_iter = max_iter,
                 eps = eps, ncpus = ncpus, verbose = verbose)
  }

  if (verbose) {
    cat(paste(">>>", paste(rep("-", 3), collapse = ""),
              "Prioritization  End",
              paste(rep("-", 31), collapse = ""), "<<<", sep = ""), "\n")
  }
  return(revl)
}
