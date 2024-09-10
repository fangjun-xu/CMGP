#' A Comprehensive Mixed-linear-model-based Tool for Genomic Prediction
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
#' @param AMB.region region size for genome dividing, only used when model was "AMB"
#' @param LRT.threshold significant level for likelihood ratio test, only used when model was "AMB
#' @param model c("Base", "GF", "MA", "KW", "AMB", "MAPS"), the first is default
#' @param part.method c("OLS", "ebvTrans", "h2"), the first is default
#' @param interval_S interval for search the optimal number of stratification, only used when model was "MAPS"
#' @param part.threshold significant level for statistic test, only used when model was "GF" or "MA
#' @param LD.threshold threshold for clustering when LD removing is working
#' @param doTrans logical, FALSE is default, if TRUE, gwas by linear transform will work
#' @param get.KR logical, should the Kinship be output
#' @param get.P see ?gaston::lmm.aireml
#' @param EMsteps see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param ncpus integer value for number of threads used
#' @param verbose logical, should the running progress be printed
#'
#' @return list
#' @export
#'
#' @examples
#' #\donttest{
#' #CMGP <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
#' #weight = NULL, alpha = c(-3, 1), ldscore = NULL, bin = 1e6,
#' #AMB.region = 1e6, LRT.threshold = 0.01,
#' #model = c("Base", "GF", "MA", "KW", "AMB", "MAPS"),
#' #part.method = c("OLS", "ebvTrans", "h2"), interval_S = c(1,15),
#' #part.threshold = 0.01, LD.threshold = 0.7, doTrans = FALSE,
#' #get.KR = FALSE, get.P = FALSE, EMsteps = 0L, eps = 1e-02,
#' #EMsteps_fail = 10L, EM_alpha = 1, max_iter = 50L,
#' #ncpus = 1, verbose = TRUE)
#' #}
CMGP <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
                 weight = NULL, alpha = c(-3, 1), ldscore = NULL, bin = 1e6,
                 AMB.region = 1e6, LRT.threshold = 0.01,
                 model = c("Base", "GF", "MA", "KW", "AMB", "MAPS"),
                 part.method = c("OLS", "ebvTrans", "h2"), interval_S = c(1,15),
                 part.threshold = 0.01, LD.threshold = 0.7, doTrans = FALSE,
                 get.KR = FALSE, get.P = FALSE, EMsteps = 0L, eps = 1e-02,
                 EMsteps_fail = 10L, EM_alpha = 1, max_iter = 50L,
                 ncpus = 1, verbose = TRUE) {
  #---------------------------------------------------------#
  #Comprehensive Mixed-linear-model-based Tool for Genomic Prediction

  CMGP.version(verbose = verbose)
  model <- model[1]
  part.method <- part.method[1]
  if (bigmemory::is.big.matrix(geno)) {
    options(bigmemory.allow.dimnames = TRUE)
  }
  rownames(geno) <- as.character(map[, 1])
  map[, 1] <- rownames(geno)
  rownames(map) <- rownames(geno)
  if (verbose) {
    cat(paste(paste(rep("=", 14), collapse = ""),
              "CMGP Start(",
              paste(Sys.Date(), rev(strsplit(date(), " ")[[1]])[2]),
              ")",
              paste(rep("=", 14), collapse = ""), sep = ""),
        "\n")
    t1 <- as.numeric(Sys.time())
    n <- dim(geno)[2]
    m <- dim(geno)[1]
    cat("Genotype data has", n, "individuals,", m, "markers\n")
    cat("Phenotype data has", length(which(!is.na(y))), "reference,",
        length(which(is.na(y))), "validation\n")
    if (!is.null(CV)) {
      cat("Number of covariates:", ncol(CV), "\n")
    }
    if (!is.null(random)) {
      cat("Number of additional random effects:", length(random), "\n")
    }
    cat("Number of threads used:", ncpus, "\n")
    cat("Analyzed model:", model, "\n")
  }

  if (model == "Base") {
    revl <- CMGP.base(y = y, CV = CV, geno = geno, random = random,
                      EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
                      EM_alpha = EM_alpha, max_iter = max_iter,
                      eps = eps, doTrans = doTrans, get.KR = get.KR,
                      get.P = get.P, verbose = verbose)
  }
  if (model == "GF") {
    revl <- CMGP.GF(y = y, CV = CV, geno = geno, map = map, random = random,
                    alpha = alpha, ldscore = ldscore, bin = bin, ncpus = ncpus,
                    part.method = part.method, part.threshold = part.threshold,
                    LD.threshold = LD.threshold, doTrans = doTrans,
                    get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
                    EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
                    max_iter = max_iter, eps = eps, verbose = verbose)
  }
  if (model == "MA") {
    revl <- CMGP.MA(y = y, CV = CV, geno = geno, map = map, random = random,
                    alpha = alpha, ldscore = ldscore, bin = bin, ncpus = ncpus,
                    part.method = part.method, part.threshold = part.threshold,
                    LD.threshold = LD.threshold, doTrans = doTrans,
                    get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
                    EMsteps = EMsteps, EMsteps_fail = EMsteps_fail,
                    max_iter = max_iter, eps = eps, verbose = verbose)
  }
  if (model == "KW") {
    revl <- CMGP.KW(y = y, CV = CV, geno = geno, map = map, random = random,
                    alpha = alpha, ldscore = ldscore, bin = bin,
                    weight = weight, W.method = part.method,
                    doTrans = doTrans, get.KR = get.KR, get.P = get.P,
                    EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                    EMsteps = EMsteps, max_iter = max_iter, eps = eps,
                    ncpus = ncpus, verbose = verbose)
  }
  if (model == "AMB") {
    revl <- CMGP.AMB(y = y, CV = CV, geno = geno, map = map, random = random,
                     bin = AMB.region, LRT.threshold = LRT.threshold,
                     doTrans = doTrans, get.KR = get.KR, get.P = get.P,
                     EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                     EMsteps = EMsteps, max_iter = max_iter, eps = eps,
                     ncpus = ncpus, verbose = verbose)
  }
  if (model == "MAPS") {
    revl <- CMGP.MAPS(y = y, CV = CV, geno = geno, map = map, random = random,
                      alpha = alpha, ldscore = ldscore, bin = bin,
                      ncpus = ncpus, interval_S = interval_S, doTrans = doTrans,
                      get.KR = get.KR, get.P = get.P, EM_alpha = EM_alpha,
                      EMsteps_fail = EMsteps_fail, EMsteps = EMsteps,
                      max_iter = max_iter, eps = eps, verbose = TRUE)
  }

  if (verbose) {
    t2 <- as.numeric(Sys.time())
    cat("Time used:", time.trans(round(t2 - t1)), "\n")
    cat(paste(paste(rep("=", 10), collapse = ""),
              "CMGP  ACCOMPLISHED(",
              paste(Sys.Date(), rev(strsplit(date(), " ")[[1]])[2]),
              ")",
              paste(rep("=", 10), collapse = ""), sep = ""),
        "\n")
  }
  return(revl)
}
