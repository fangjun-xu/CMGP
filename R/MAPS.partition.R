#' Sratification of markers based on PVE
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param map geno information, SNP, Chr, Pos
#' @param random list, non-additive genetic random effect
#' @param alpha genetype scaling value
#' @param ldscore LD score
#' @param bin region size for correlation calculation
#' @param interval_S interval of number of stratification when optimizing
#' @param ncpus integer value, default as 1, how many threads used
#' @param verbose logical, TRUE for print the output
#' @param EMsteps see ?gaston::lmm.aireml
#' @param EMsteps_fail see ?gaston::lmm.aireml
#' @param EM_alpha see ?gaston::lmm.aireml
#' @param max_iter see ?gaston::lmm.aireml
#' @param eps see ?gaston::lmm.aireml
#'
#' @return list
#' @export
#'
#' @examples
#' #\donttest{
#' #mp <- MAPS.partition(y = y, CV = CV, geno = geno, map = map, random = random,
#' #alpha = alpha, ldscore = ldscore, bin = bin,
#' #interval_S = interval_S, EMsteps = EMsteps, eps = eps,
#' #EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
#' #max_iter = max_iter, ncpus = ncpus, verbose = verbose)
#' #}
MAPS.partition <- function(y = NULL, CV = NULL, geno = NULL, map = NULL,
                           random = NULL, alpha = 0, ldscore = NULL, bin = 1e6,
                           interval_S = c(1,15), ncpus = 1, verbose = TRUE,
                           EMsteps = 0L, EMsteps_fail = 10L, EM_alpha = 1,
                           max_iter = 50L, eps = 1e-02) {

  #---------------------------------------------------------#
  #Stratification of SNP
  if (bigmemory::is.big.matrix(geno)) {
    options(bigmemory.allow.dimnames = TRUE)
  }
  rownames(geno) <- as.character(map[, 1])
  map[, 1] <- rownames(geno)
  rownames(map) <- rownames(geno)
  sp <- SNP.prioritized(y = y, CV = CV, geno = geno, map = map,
                        random = random, alpha = alpha, ldscore = ldscore,
                        bin = bin, method = "h2", EMsteps = EMsteps,
                        EMsteps_fail = EMsteps_fail, EM_alpha = EM_alpha,
                        max_iter = max_iter, eps = eps,
                        ncpus = ncpus, verbose = verbose)

  if (verbose) {
    cat(paste(">>>", paste(rep("-", 3), collapse = ""),
              "Stratification  Start",
              paste(rep("-", 29), collapse = ""), "<<<", sep = ""), "\n")
  }

  #pve <- as.vector(sp[, ncol(sp)])
  #names(pve) <- rownames(map)
  pve <- as.data.frame(sp[, 7:9])
  rownames(pve) <- rownames(map)

  inter.s <- c(min(interval_S) : max(interval_S))
  if (verbose) {
    if (length(inter.s) == 1) {
      t1 <- as.numeric(Sys.time())
      cat("Stratified genome into", inter.s, "groups\n")
    }
    else {
      t1 <- as.numeric(Sys.time())
      cat("Optimizing stratifications by k-means...\n")
    }
  }

  M <- dim(geno)[1]
  st.loop <- function(i) {
    set.seed(i)
    part <- stats::kmeans(pve, centers = i, nstart = 50)$cluster
    genolist <- lapply(1:i, function(j) {
      snpij <- as.numeric(which(part == j))
      gij <- matrix(geno[snpij, ], nrow = length(snpij))
      colnames(gij) <- colnames(geno)
      rownames(gij) <- rownames(geno)[snpij]
      return(gij)
    })
    #wi <- pve
    wi <- as.vector(pve[, ncol(pve)])
    names(wi) <- rownames(map)

    mix <- EBV.trans(y = y, CV = CV, geno = genolist, weight = wi,
                     random = random, EMsteps = EMsteps, EM_alpha = EM_alpha,
                     EMsteps_fail = EMsteps_fail, max_iter = max_iter,
                     eps = eps, get.KR = FALSE, get.P = FALSE,
                     doTrans = FALSE, verbose = FALSE)
    hg <- mix[[1]]$vg / sum(mix[[1]]$vg, mix[[1]]$vr, mix[[1]]$ve)
    m <- unlist(lapply(genolist, nrow))
    enrich <- M * hg / m
    ll <- mix[[1]]$logL
    if (verbose) {
      cat("  [", "Num. of Stra.:", i, ";", "LogL:", ll, "]\n")
    }
    return(list(part, ll, enrich))
  }

  if (length(inter.s) == 1 || ncpus == 1) {
    loop <- lapply(inter.s, st.loop)
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      loop <- parallel::mclapply(inter.s, st.loop, mc.cores = ncpus,
                                 mc.preschedule = FALSE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      geno <- geno[, ]
      suppressMessages(snowfall::sfExport("pve", "geno", "y",
                                          "CV", "random", "M",
                                          "EMsteps", "EM_alpha",
                                          "EMsteps_fail", "max_iter",
                                          "eps", "verbose"))
      suppressMessages(snowfall::sfLibrary(stats))
      suppressMessages(snowfall::sfLibrary(CMGP))
      loop <- snowfall::sfLapply(inter.s, st.loop)
      suppressMessages(snowfall::sfStop())
    }
  }
  ll <- unlist(lapply(loop, function(l) { return(l[[2]]) }))
  index <- min(which(ll == max(ll, na.rm = TRUE)))
  ops <- inter.s[index]
  partition <- loop[[index]][[1]]
  h2_enrich <- loop[[index]][[3]]
  if (verbose) {
    t2 <- as.numeric(Sys.time())
    cat(paste("  | 100% elapsed=", time.trans(round(t2 - t1)),sep = ""), "\n")
    if (length(inter.s) != 1) {
      cat("Optimal num. of stra.:", ops ,"\n")
    }
    cat("SNPs Num. in Stras. / Enrichment h2:\n")
    nums <- as.vector(table(partition))
    for (i in 1:length(nums)){
      cat(" ", paste("[Stra.", i, "]", sep = ""), nums[i], "/", h2_enrich[i], "\n")
    }
  }
  revl <- cbind(sp, partition)
  colnames(revl) <- c(colnames(sp), "Cluster")

  if (verbose) {
    cat(paste(">>>", paste(rep("-", 3), collapse = ""),
              "Stratification  End",
              paste(rep("-", 31), collapse = ""), "<<<", sep = ""), "\n")
  }
  return(list(Cluster.map = revl, Enrichment.h2 = h2_enrich))
}
