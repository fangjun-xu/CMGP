#' Parallel calculate the LD score
#'
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param map geno information, SNP, Chr, Pos
#' @param bin region size for correlation calculation
#' @param ncpus integer value, default as 1, how many threads used
#' @param verbose logical, TRUE for print the output
#'
#' @return vertor
#' @export
#'
#' @examples
#' #\donttest{
#' #ldscore <- LDSC(geno = geno, map = map, bin = bin,
#' #ncpus = ncpus, verbose = verbose)
#' #}
LDSC <- function(geno = NULL, map = NULL, bin = 1e6,
                 ncpus = 1, verbose = TRUE) {
  if (dim(geno)[1] > ncpus && ncpus > 1) {
    pbseq <- seq(dim(geno)[1], 1, -ceiling(dim(geno)[1] / ncpus))[1:ncpus]
  }else {
    pbseq <- 1:dim(geno)[1]
  }
  if (verbose) {
    cat("Calculating LD score...\n")
    cat("Bin size for correlation:", bin / 1000, "Kb\n")
    maxpb <- length(pbseq) + 1
    pb <- pbapply::timerProgressBar(max = maxpb,
                                    width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  if (bigmemory::is.big.matrix(geno)) {
    options(bigmemory.allow.dimnames = TRUE)
  }
  rownames(geno) <- as.character(map[, 1])
  map[, 1] <- rownames(geno)
  rownames(map) <- rownames(geno)
  ldCal <- function(i) {
    chri <- as.numeric(map[i, 2])
    posi <- as.numeric(map[i, 3])
    mapLD <- map[map[, 2] == chri, ]
    mapLD <- mapLD[abs(posi - mapLD[, 3]) <= bin, ]
    index <- mapLD[, 1]
    gi <- geno[i, ]
    gLD <- geno[index, ]
    scor <- stats::cor(gi, t(gLD)) ^ 2
    adjcor <- scor - ((1 - scor) / (dim(geno)[2] - 2))
    LDscore <- sum(adjcor) - 1
    if (verbose && i %in% pbseq) {
      pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
    }
    if (LDscore < 0) {
      LDscore <- 0
    }
    return(LDscore)
  }
  if (dim(geno)[1] == 1 || ncpus == 1) {
    ldsc <- lapply(1:dim(geno)[1], FUN = ldCal)
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      ldsc <- parallel::mclapply(1:dim(geno)[1], FUN = ldCal,
                                 mc.cores = ncpus,
                                 mc.preschedule = TRUE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      geno <- geno[, ]
      suppressMessages(snowfall::sfExport("geno", "map", "bin",
                                          "pbseq", "pb", "verbose"))
      suppressMessages(snowfall::sfLibrary(pbapply))
      suppressMessages(snowfall::sfLibrary(stats))
      ldsc <- snowfall::sfLapply(1:dim(geno)[1], fun = ldCal)
      suppressMessages(snowfall::sfStop())
    }
  }
  if (verbose) {
    pbapply::setTimerProgressBar(pb, maxpb)
  }
  ldsc <- unlist(ldsc)
  return(ldsc)
}
