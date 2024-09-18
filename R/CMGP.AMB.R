 #' Adaptive Multi-BLUP
#'
#' @param y phenotype, n*1 vector, NAs refer to validation
#' @param CV covariates, n*q matrix, numeric format for quantitative, character or factor for discrete
#' @param geno genotype, m * n matrix/bigmatrix
#' @param map marker information, SNP, Chr, Pos
#' @param random list, non-additive genetic effects
#' @param bin region size for genome dividing
#' @param LRT.threshold significant level for likelihood ratio test
#' @param doTrans logical, FALSE is default, if TRUE, gwas by linear transform will work
#' @param get.KR logical, should the Kinship be output
#' @param get.P see ?gaston::lmm.aireml
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#' @param ncpus integer value for number of threads used
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
CMGP.AMB <- function(y = NULL, CV = NULL, geno = NULL, map = NULL, random = NULL,
                     bin = 1e6, LRT.threshold = 0.01, doTrans = FALSE,
                     get.KR = FALSE, get.P = FALSE, EMsteps = 0L,
                     EMsteps_fail = 10L, EM_alpha = 1, max_iter = 50L,
                     eps = 1e-02, ncpus = 1, verbose = TRUE) {
  #---------------------------------------------------------#
  #AMB (Adaptive Multi-BLUP) for Genomic Prediction
  if (bigmemory::is.big.matrix(geno)) {
    options(bigmemory.allow.dimnames = TRUE)
  }
  rownames(geno) <- as.character(map[, 1])
  map[, 1] <- rownames(geno)
  rownames(map) <- rownames(geno)
  if (verbose) {
    cat(paste(paste(rep("-", 22), collapse = ""),
              "AMBLUP  Running",
              paste(rep("-", 22), collapse = ""), sep = ""), "\n")
  }

  genolist <- genome_cut(geno = geno, map = map, bin = bin,
                         ncpus = ncpus, verbose = verbose)
  if (verbose) {
    cat("Sorting the map...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  amb.map <- lapply(1:length(genolist), function(i) {
    ind <- rownames(genolist[[i]])
    mapi <- map[ind, ]
    mapi$Region <- rep(i, length(ind))
    return(mapi)
  })
  amb.map <- do.call(rbind, amb.map)
  if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
  }
  if (length(genolist) > ncpus && ncpus > 1) {
    pbseq <- seq(length(genolist), 1,
                 -ceiling(length(genolist) / ncpus))[1:ncpus]
  }else {
    pbseq <- 1:length(genolist)
  }
  if (verbose) {
    cat("Likelihood ratio test for genome regions...\n")
    cat("Significant level:", LRT.threshold / length(genolist), "\n")
    maxpb <- length(pbseq) + 1
    pb <- pbapply::timerProgressBar(max = maxpb,
                                    width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  L0 <- EBV.trans(y = y, CV = CV, geno = list(geno), random = random,
                  ncpus = 1, max_iter = max_iter, eps = eps,
                  EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                  EMsteps = EMsteps, get.KR = FALSE, get.P = FALSE,
                  doTrans = FALSE, verbose = FALSE)$GP$logL

  logl_1 <- function(i) {
    fg <- genolist[[i]]
    bg <- genolist[-i]
    bg <- do.call(rbind, bg)
    ll <- EBV.trans(y = y, CV = CV, geno = list(fg, bg), random = random,
                    ncpus = 1, max_iter = max_iter, eps = eps,
                    EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                    EMsteps = EMsteps, get.KR = FALSE, get.P = FALSE,
                    doTrans = FALSE, verbose = FALSE)$GP$logL
    if (verbose && i %in% pbseq) {
      pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
    }
    return(ll)
  }
  if (ncpus == 1 || length(genolist) == 1) {
    L1 <- lapply(1:length(genolist), logl_1)
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      L1 <- parallel::mclapply(1:length(genolist), logl_1,
                               mc.cores = ncpus,
                               mc.preschedule = length(genolist) > 2 * ncpus)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      suppressMessages(snowfall::sfExport("genolist", "y", "CV",
                                          "random", "max_iter", "EMsteps",
                                          "EM_alpha", "EMsteps_fail", "eps",
                                          "pbseq", "pb", "verbose"))
      suppressMessages(snowfall::sfLibrary(CMGP))
      suppressMessages(snowfall::sfLibrary(pbapply))
      L1 <- snowfall::sfLapply(1:length(genolist), fun = logl_1)
      suppressMessages(snowfall::sfStop())
    }
  }
  L1 <- unlist(L1)
  chis <- 2 * (L1 - L0)
  Pvalue <- stats::pchisq(chis, 1, lower.tail = FALSE)
  Pvalue[which(is.na(Pvalue))] <- 1
  chis.test <- cbind(L1, chis, Pvalue)
  ambp <- list(Region = amb.map, Statistic = chis.test)
  index <- which(Pvalue <= (LRT.threshold / length(genolist)))
  if (verbose) {
    pbapply::setTimerProgressBar(pb, maxpb)
    cat("\n")
    cat(length(index), "regions are significant\n")
  }
  if (length(index) == 0) {
    index <- which(Pvalue <= stats::quantile(Pvalue,
                                             LRT.threshold,
                                             na.rm = TRUE))
    if (verbose) {
      cat("Selected the top", length(index), "regions\n")
    }
  }

  if (verbose) {
    cat("Merging adjacent regions...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  index <- sort(index)
  if (length(index) == length(geno)) {
    geno.bg <- NULL
  }else {
    geno.bg <- do.call(rbind, genolist[-index])
  }
  if (length(index) == 1) {
    geno.fg <- genolist[index]
  }else {
    ind.mer <- rep(1, length(index))
    for (i in 2:length(index)) {
      if (index[i] == (index[i - 1] + 1)) {
        ind.mer[i] <- ind.mer[i - 1]
      }else {
        ind.mer[i] <- ind.mer[i - 1] + 1
      }
    }
    geno.fg <- lapply(unique(ind.mer), function(i) {
      exi <- index[which(ind.mer == i)]
      gi <- do.call(rbind, genolist[exi])
      gi <- as.matrix(gi)
      return(gi)
    })
  }
  if (is.null(geno.bg)) {
    amb.geno <- geno.fg
  }else {
    amb.geno <- c(geno.fg, list(geno.bg))
  }
  if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
    cat(length(amb.geno), "regions obtained after mergence\n")
  }

  amb.fit <- EBV.trans(y = y, CV = CV, geno = amb.geno, random = random,
                       ncpus = ncpus, max_iter = max_iter,
                       EM_alpha = EM_alpha, EMsteps_fail = EMsteps_fail,
                       EMsteps = EMsteps, get.KR = get.KR, get.P = get.P,
                       eps = eps, doTrans = doTrans, verbose = verbose)
  if (doTrans) {
    amb.fit$GWAS <- amb.fit$GWAS[rownames(geno), ]
  }

  if (verbose) {
    cat(paste(paste(rep("-", 23), collapse = ""),
              "AMBLUP FINISH",
              paste(rep("-", 23), collapse = ""), sep = ""), "\n")
  }
  revl <- list(GP = amb.fit$GP, GWAS = amb.fit$GWAS, AMB.cluster = ambp)
  return(revl)
}
