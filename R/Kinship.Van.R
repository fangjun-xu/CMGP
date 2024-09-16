#' Calculate the Kinship in VanRaden's method
#'
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param weight numeric vector, the weight of markers
#' @param ncut integer value, default as 1, how many block should the geno matrix be divided into?
#' @param ncpus integer value, default as 1, how many threads used
#'
#' @return genetic relationship matrix
#' @export
#'
#' @examples
#' #\donttest{
#' #K <- Kinship.Van(geno=geno)
#' #}
Kinship.Van <- function(geno = NULL, weight = NULL, ncut = 1, ncpus = 1) {

  #---------------------------------------------------------#
  #Calculate the genomic Kinship by VanRaden's method
  if (!bigmemory::is.big.matrix(geno) && !is.matrix(geno)) {
    stop("Wrong format of geno inputed! matrix/big.matrix is allowed.")
  }
  if (sum(is.na(geno[, ])) != 0) stop("'NA' isn't allowed in geno!")
  if (!is.null(weight)) {
    if (sum(is.na(weight)) != 0) stop("'NA' isn't allowed in weight!")
  }
  if (ncut < 1 || ncut %% 1 != 0) stop("ncut should be a positive integer!")

  n <- dim(geno)[2] #NUM OF INDIVIDUALS
  m <- dim(geno)[1] #NUM OF MARKERS
  if (m == 1) {
    maf <- mean(geno[, ]) / 2
  }else {
    maf <- rowMeans(geno[, ]) / 2
  }
  Z <- bigmemory::big.matrix(m, n, type = "double",
                             init = NULL, shared = FALSE)
  if (!is.null(weight)) {
    Z[, ] <- (geno[, ] - 2 * maf) * sqrt(weight)
  }else {
    Z[, ] <- geno[, ] - 2 * maf
  }
  rm(geno, weight)
  K <- matrix(NA, n, n)
  if (ncut == 1 || ncpus == 1) {
    ##recommend##
    if (m == 1) {
      K[, ] <- tcrossprod(Z[, ]) / (2 * sum(maf * (1 - maf)))
    }else {
      K[, ] <- crossprod(Z[, ]) / (2 * sum(maf * (1 - maf)))
    }
    rm(Z)
  }else {
    cut.lab <- cut(1:n, breaks = ncut, labels = FALSE)
    cut.list <- list()
    for (i in 1:ncut) {
      for (j in i:ncut) {
        cuti <- which(cut.lab == i)
        cutj <- which(cut.lab == j)
        cut.list <- c(cut.list, list(list(cuti, cutj)))
      }
    }
    rm(cut.lab, cuti, cutj)
    block.prod <- function(i) {
      row <- cut.list[[i]][[1]]
      col <- cut.list[[i]][[2]]
      block <- crossprod(Z[, row], Z[, col])
      return(block)
    }
    if (Sys.info()[["sysname"]] == "Linux") {
      bloc <- parallel::mclapply(1:length(cut.list), FUN = block.prod,
                                 mc.cores = ncpus,
                                 mc.preschedule = FALSE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      Z <- Z[, ]
      suppressMessages(snowfall::sfExport("Z", "cut.list"))
      bloc <- snowfall::sfLapply(1:length(cut.list), fun = block.prod)
      suppressMessages(snowfall::sfStop())
    }
    rm(Z)
    for (i in 1:length(cut.list)) {
      row <- cut.list[[i]][[1]]
      col <- cut.list[[i]][[2]]
      K[row, col] <- bloc[[i]]
      if (!identical(row, col)) {
        K[col, row] <- t(bloc[[i]])
      }
    }
    K[, ] <- K[, ] / (2 * sum(maf * (1 - maf)))
  }
  return(K)
}
