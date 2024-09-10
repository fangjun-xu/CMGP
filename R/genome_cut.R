#' Title
#'
#' @param geno genotype, m * n matrix/bigmatrix
#' @param map marker information, SNP, Chr, Pos
#' @param bin region size for genome dividing
#' @param ncpus integer value for number of threads used
#' @param verbose logical, should the running progress be printed
#'
#' @return divided genome list
#' @export
#'
#' @examples
#' #\donttest{
#' #genolist <- genome_cut(geno = geno, map = map, bin = bin,
#' #ncpus = ncpus, verbose = verbose)
#' #}
genome_cut <- function(geno = NULL, map = NULL, bin = 1e6,
                       ncpus = 1, verbose = TRUE) {

  if (bigmemory::is.big.matrix(geno)) {
    options(bigmemory.allow.dimnames = TRUE)
  }
  rownames(geno) <- as.character(map[, 1])
  map[, 1] <- rownames(geno)
  rownames(map) <- rownames(geno)
  chr <- unique(map[, 2])
  if (length(chr) > ncpus) {
    pbseq <- seq(length(chr), 1, -ceiling(length(chr) / ncpus))[1:ncpus]
  }else {
    pbseq <- 1:length(chr)
  }
  if (verbose) {
    cat("Dividing the genome regions...\n")
    cat("Region size:", bin / 1000, "Kb\n")
    maxpb <- min(ncpus, length(chr)) + 1
    pb <- pbapply::timerProgressBar(max = maxpb ,
                                    width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  cut_ge <- function(i) {
    chi <- map[map[, 2] == chr[i], ]
    st <- min(chi[, 3]) - 1
    ed <- max(chi[, 3])
    break_seq <- unique(c(seq(st, ed, bin), ed))
    clus <- cut(chi[, 3], breaks = break_seq, labels = FALSE)
    gi <- lapply(unique(clus), function(j) {
      ind <- chi[clus == j, 1]
      return(geno[ind, ])
    })
    if (verbose && i %in% pbseq) {
      pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% i))
    }
    return(gi)
  }
  if (ncpus == 1 || length(chr) == 1) {
    genolist <- lapply(1:length(chr), cut_ge)
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      genolist <- parallel::mclapply(1:length(chr), cut_ge,
                                     mc.cores = ncpus,
                                     mc.preschedule = FALSE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      suppressMessages(snowfall::sfExport("geno", "map", "bin", "pbseq",
                                          "chr", "pb", "verbose"))
      suppressMessages(snowfall::sfLibrary(pbapply))
      L1 <- snowfall::sfLapply(1:length(chr), cut_ge)
      suppressMessages(snowfall::sfStop())

    }
  }
  genolist <- do.call(c, genolist)
  if (verbose) {
    pbapply::setTimerProgressBar(pb, maxpb)
  }
  if (verbose) {
    cat("\n")
    cat(length(genolist), "genome regions obtained")
  }
  return(genolist)
}
