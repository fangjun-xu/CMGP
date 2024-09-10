#' Marker Adaptive Prioritization and Stratification
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
#' @param interval_S interval for search the optimal number of stratification, only used when model was "MAPS"
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
#' #revl <- CMGP.MAPS(y = y, CV = CV, geno = geno, map = map, random = random,
#' #alpha = alpha, ldscore = ldscore, bin = bin,
#' #ncpus = ncpus, interval_S = interval_S, doTrans = doTrans,
#' #get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
#' #EMsteps_fail = EMsteps_fail, EMsteps = EMsteps,
#' #max_iter = max_iter, eps = eps, verbose = TRUE)
#' #}
CMGP.MAPS <- function(y = NULL, CV = NULL, geno = NULL, map = NULL,
                      random = NULL, alpha = c(-3, 1), ldscore = NULL, bin = 1e6,
                      ncpus = 1, interval_S = c(1,15), doTrans = FALSE,
                      get.KR = FALSE, get.P = FALSE, EMsteps = 0L,
                      EMsteps_fail = 10L, EM_alpha = 1,
                      max_iter = 50L, eps = 1e-02, verbose = TRUE) {

  #---------------------------------------------------------#
  #MAPS (Marker adaptive prioritization and stratification) for Genomic Prediction
  if (verbose) {
    cat(paste(paste(rep("-", 23), collapse = ""),
              "MAPS  Running",
              paste(rep("-", 23), collapse = ""), sep = ""), "\n")
  }
  options(bigmemory.allow.dimnames = TRUE)
  mp <- MAPS.partition(y = y, CV = CV, geno = geno, map = map, random = random,
                       alpha = alpha, ldscore = ldscore, bin = bin,
                       interval_S = interval_S, EMsteps = EMsteps, eps = eps,
                       EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
                       max_iter = max_iter, ncpus = ncpus, verbose = verbose)
  clusm <- as.data.frame(mp$Cluster.map)
  part <- as.numeric(clusm$Cluster)
  maps.geno <- lapply(1:max(part), function(i){
    index <- as.numeric(which(part == i))
    gi <- bigmemory::big.matrix(length(index), dim(geno)[2], type = "double",
                                init = NULL, shared = FALSE)
    gi[, ] <- geno[index, ]
    rownames(gi) <- rownames(geno)[index]
    colnames(gi) <- colnames(geno)
    return(gi)
  })
  weight <- clusm$PVE
  names(weight) <- rownames(geno)
  rm(clusm, part)
  maps.fit <- EBV.trans(y = y, CV = CV, geno = maps.geno, random = random,
                        weight = weight, ncpus = ncpus, max_iter = max_iter,
                        EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                        EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                        eps = eps, doTrans = doTrans, verbose = verbose)
  if (doTrans) {
    maps.fit$GWAS <- maps.fit$GWAS[rownames(geno), ]
  }
  if (verbose) {
    cat(paste(paste(rep("-", 24), collapse = ""),
              "MAPS FINISH",
              paste(rep("-", 24), collapse = ""), sep = ""), "\n")
  }
  revl <- list(GP = maps.fit$GP, GWAS = maps.fit$GWAS, MAPS.cluster = mp)
  return(revl)
}
