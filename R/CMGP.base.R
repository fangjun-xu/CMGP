#' GBLUP model for genomic prediction
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param random list, non-additive genetic random effect
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param get.KR logical, should the Kinship be returned
#' @param get.P see ?gaston::lmm.aireml
#' @param doTrans logical, should marker effects be calculated
#' @param verbose logical, TRUE for print the output
#'
#' @return list
#' @export
#'
#' @examples
#' #\donttest{
#' #revl <- CMGP.base(y = y, CV = CV, geno = geno, random = random,
#' #EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
#' #EM_alpha = EM_alpha, max_iter = max_iter,
#' #eps = eps, doTrans = doTrans, get.KR = get.KR,
#' #get.P = get.P, verbose = verbose)
#' #}
CMGP.base <- function(y = NULL, CV = NULL, geno = NULL, random = NULL,
                      EMsteps = 0L, EMsteps_fail = 10L, EM_alpha = 1,
                      max_iter = 50L, eps = 1e-02, doTrans = FALSE,
                      get.KR = FALSE, get.P = FALSE, verbose = TRUE) {

  #---------------------------------------------------------#
  #Basic Model(GBLUP) for Genomic Prediction
  if (verbose) {
    cat(paste(paste(rep("-", 23), collapse = ""),
              "GBLUP Running",
              paste(rep("-", 23), collapse = ""), sep = ""), "\n")
  }

  base.fit <- EBV.trans(y = y, CV = CV, geno = list(geno), random = random,
                        ncpus = 1, max_iter = max_iter, eps = eps,
                        EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                        EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                        doTrans = doTrans, verbose = verbose)

  if (verbose) {
    cat(paste(paste(rep("-", 23), collapse = ""),
              "GBLUP  FINISH",
              paste(rep("-", 23), collapse = ""), sep = ""), "\n")
  }
  return(base.fit)
}
