#' Solve the mixed model using AIREML and calculate marker effects by linear transform
#'
#' @param y phenotype vector, NAs refer to validation
#' @param CV covariates
#' @param geno Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size
#' @param weight numeric vector, the weight of markers, using marker names to name it
#' @param random list, non-additive genetic random effect
#' @param ncpus integer value, default as 1, how many threads used
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
#' @return list, GP and GWAS
#' @export
#'
#' @examples
#' #\donttest{
#' #mix <- EBV.trans(y = y, CV = CV, geno = genolist, weight = pve,
#' #random = random, EMsteps = EMsteps, EM_alpha = EM_alpha,
#' #EMsteps_fail = EMsteps_fail, max_iter = max_iter,
#' #eps = eps, get.KR = FALSE, get.P = FALSE,
#' #doTrans = FALSE, verbose = TRUE)
#' #}
EBV.trans <- function(y = NULL, CV = NULL, geno = NULL, weight = NULL,
                      random = NULL, ncpus = 1, EMsteps = 0, EMsteps_fail = 10,
                      EM_alpha = 1, max_iter = 50, eps = 1e-02, get.KR = FALSE,
                      get.P = FALSE, doTrans = FALSE, verbose = TRUE) {
  #---------------------------------------------------------#
  #Solve the mixed model and transform the EBV to marker effects

  if (!is.null(CV) && length(y) != dim(CV)[1]) {
    stop("Number of individuals don't match between phenotype and covariates!")
  }
  if (!is.null(random) && length(y) != dim(random)[1]) {
    stop("Number of individuals don't match between phenotype and covariates!")
  }
  if (!is.list(geno)) {
    warning("geno should be provided as list format.")
    geno <- list(geno)
  }
  for (i in 1:length(geno)) {
    gi <- geno[[i]]
    if (!bigmemory::is.big.matrix(gi) && !is.matrix(gi)) {
      stop("Wrong format in geno list inputed! matrix/big.matrix is allowed.")
    }
    if (sum(is.na(gi[, ])) != 0) stop("'NA' isn't allowed in geno!")
    if (length(y) != dim(gi)[2]) {
      stop("Number of individuals don't match between y and geno!", "\n",
           "In addition: geno should be (marker size) * (individual size).")
    }
  }
  rm(gi)
  if (sum(is.na(weight)) != 0) stop("'NA' isn't allowed in weight!")

  ng <- length(geno)
  ref <- which(!is.na(y))
  val <- which(is.na(y))
  n <- length(y)
  CV <- CovMatrix(CV = CV, n = n)

  if (verbose) {
    cat("Calculating kinship(", ng, "group )...\n")
    cat("  |")
    t1 <- as.numeric(Sys.time())
  }
  if (ng == 1 || ncpus == 1) {
    Klist <- lapply(1:ng, function(i) {
      gi <- geno[[i]]
      wi <- weight[row.names(gi)]
      if (!is.null(wi)) {
        wi <- length(wi) * wi / sum(wi)
      }
      Ki <- Kinship.Van(geno = gi, weight = wi, ncut = 1, ncpus = 1)
      rm(gi)
      if (verbose) cat("-")
      return(Ki)
    })
  }else {
    if (Sys.info()[["sysname"]] == "Linux") {
      Klist <- parallel::mclapply(1:ng, function(i) {
        gi <- geno[[i]]
        wi <- weight[row.names(gi)]
        if (!is.null(wi)) {
          wi <- length(wi) * wi / sum(wi)
        }
        Ki <- Kinship.Van(geno = gi, weight = wi, ncut = 1, ncpus = 1)
        rm(gi)
        if (verbose) cat("-")
        return(Ki)
      }, mc.cores = ncpus, mc.preschedule = FALSE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      suppressMessages(snowfall::sfExport("geno", "weight",
                                          "verbose", "Kinship.Van"))
      Klist <- snowfall::sfLapply(1:ng, function(i) {
        gi <- geno[[i]]
        wi <- weight[row.names(gi)]
        if (!is.null(wi)) {
          wi <- length(wi) * wi / sum(wi)
        }
        Ki <- Kinship.Van(geno = gi, weight = wi, ncut = 1, ncpus = 1)
        rm(gi)
        if (verbose) cat("-")
        return(Ki)
      })
      suppressMessages(snowfall::sfStop())
    }
  }
  if (verbose) {
    t2 <- as.numeric(Sys.time())
    if (round(t2 - t1) == 0) {
      tud <- "0s"
    }else {
      tud <- time.trans(round(t2 - t1))
    }
    cat(paste("| 100% elapsed=", tud, sep = ""), "\n")
  }
  Randlist <- Klist
  if (!is.null(random)) {
    if (verbose) cat(length(random), "additional random effect(s) inputed.\n")
    Krand <- lapply(1:length(random), function(i) {
      ri <- as.matrix(CovMatrix(CV = random[[i]], n = n)[, -1])
      Kri <- bigmemory::big.matrix(n, n, type = "double",
                                   init = NULL, shared = FALSE)
      Kri[, ] <- tcrossprod(ri) / ncol(ri)
      return(Kri)
    })
    Randlist <- c(Klist, Krand)
  }

  if (verbose) {
    cat("Solving the mixed model...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  Pcall <- TRUE
  fit <- gaston::lmm.aireml(y[ref],
                            X = as.matrix(CV[ref, ]),
                            K = lapply(Randlist, function(i) {
                              return(i[ref, ref])
                            }),
                            EMsteps = EMsteps,
                            EMsteps_fail = EMsteps_fail,
                            EM_alpha = EM_alpha,
                            max_iter = max_iter, eps = eps,
                            get.P = Pcall, verbose = FALSE)
  if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
  }
  if (verbose) {
    cat("Generating prediction result...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
  }
  ref.ebv <- lapply(1:length(Randlist), function(i) {
    ebv.ri <- crossprod(Randlist[[i]][ref, ref], as.matrix(fit$tau[i] * fit$Py))
    return(ebv.ri)
  })
  val.ebv <- lapply(1:length(Randlist), function(i) {
    ebv.vi <- crossprod(Randlist[[i]][ref, val], as.matrix(fit$tau[i] * fit$Py))
    return(ebv.vi)
  })
  ref.ebv <- matrix(unlist(ref.ebv), ncol = length(Randlist))
  val.ebv <- matrix(unlist(val.ebv), ncol = length(Randlist))

  ebvgr <- matrix(NA, nrow = n, ncol = length(Randlist))
  ebvgr[ref, ] <- ref.ebv
  ebvgr[val, ] <- val.ebv
  ebv.g <- as.matrix(ebvgr[, 1:ng])
  vg <- fit$tau[1:ng]
  if (!is.null(random)) {
    ebv.r <- as.matrix(ebvgr[, -c(1:ng)])
    vr <- fit$tau[-c(1:ng)]
  }else {
    ebv.r <- NULL
    vr <- NULL
  }
  ve <- fit$sigma2
  logL <- fit$logL
  beta <- fit$BLUP_beta

  K <- lapply(1:length(Randlist), function(i) {
    return(fit$tau[i] * Randlist[[i]][ref, ref])
  })
  V <- matrix(0, length(ref), length(ref))
  for (i in 1:length(K)) {
    V <- V + K[[i]]
  }
  V <- V + diag(ve, length(ref))
  P <- fit$P
  V_VPV <- V - tcrossprod(crossprod(V, P), V)
  rm(K, V, P)
  X <- as.matrix(CV[ref, ])
  XXiX <- try(tcrossprod(solve(crossprod(X)), X), silent = TRUE)
  if (inherits(XXiX, "try-error")) {
    warning("The CV matrix is singular!\nUsing general inverse insteaded")
    XXiX <- tcrossprod(MASS::ginv(crossprod(X)), X)
  }
  rm(X)
  varbeta <- tcrossprod(tcrossprod(XXiX, V_VPV), XXiX)
  rm(XXiX, V_VPV)

  xb <- CV %*% as.matrix(beta)
  y_pre <- as.vector(xb + apply(ebvgr, 1, sum))
  if (get.P) {
    P <- fit$P
  }else {
    P <- NULL
  }
  Py <- fit$Py
  if (get.KR) {
    KR <- Randlist
  }else {
    KR <- NULL
  }

  GP <- list(y = y,
             y_predict = y_pre,
             beta = as.vector(beta),
             varbeta = varbeta,
             gebv = ebv.g,
             rebv = ebv.r,
             vg = as.vector(vg),
             vr = as.vector(vr),
             ve = as.vector(ve),
             logL = as.vector(logL),
             P = P,
             Py = Py,
             Randomship = KR)
  if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
  }

  GWAS <- NULL
  if (doTrans) {
    if (verbose) {
      cat("Transforming the EBV to marker effects...\n")
      pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
      on.exit(close(pb))
    }
    P <- fit$P
    DF <- n - ncol(CV) - length(Randlist)
    if (ng == 1 || ncpus == 1) {
      gwas <- lapply(1:ng, function(i) {
        genoi <- geno[[i]]
        wi <- weight[row.names(genoi)]
        mafi <- rowMeans(genoi[, ]) / 2
        scalei <- 1 / (2 * sum(mafi * (1 - mafi)))
        if (is.null(wi)) {
          DZ <- (genoi[, ref] - mafi)
        }else {
          wi <- length(wi) * wi / sum(wi)
          DZ <- (genoi[, ref] - mafi) * sqrt(wi)
        }
        Ki <- Klist[[i]]
        K_inv <- try(solve(Ki[, ] + diag(1, ncol(Ki)) * (1e-10)), silent = TRUE)[ref, ref]
        if (inherits(K_inv, "try-error")) {
          warning("The GRM of geno component ", i, " is singular!\nUsing general inverse insteaded")
          K_inv <- MASS::ginv(Ki[, ])[ref, ref]
        }
        ef <- colSums(tcrossprod(K_inv, DZ) * ebv.g[, i])
        Effect <- scalei * ef
        pch <-try(chol(P), silent = TRUE)
        if (inherits(pch, "try-error")) {
          warning("The ", i, "th P matrix is not positive definite\nTry spectral decomposition instead of cholesky decomposition")
          eig <- eigen(P, symmetric = TRUE)
          vec <- eig$vectors
          val <- eig$values
          DZS <- tcrossprod(DZ, t(vec))
          diagv <- apply(DZS, 1, function(i) {
            return(sum(val * i^2))
          })
          SE <- vg[i] * scalei * sqrt(diagv)
          rm(ef, eig, vec, val, DZS, pch, diagv)
        }else {
          ZP <- tcrossprod(DZ, pch)
          n <- dim(pch)[2]
          diagv <- (n - 1) * Rfast::rowVars(ZP) + n * rowMeans(ZP) ^ 2
          SE <- vg[i] * scalei * sqrt(diagv)
          rm(ef, ZP, pch, diagv)
        }
        Pvalue <- 2 * stats::pt(abs(Effect / SE), DF, lower.tail = FALSE)
        if (sum(is.na(Pvalue)) != 0) {
          warning("Non-positive variance appeared at ", which(is.na(Pvalue)), " markers in geno group ", i)
          Pvalue[which(is.na(Pvalue))] <- 1
        }
        rm(genoi, wi, mafi, Ki, scalei, DZ, K_inv)
        return(cbind(Effect, SE, Pvalue))
      })
    }else {
      if (Sys.info()[["sysname"]] == "Linux") {
        gwas <- parallel::mclapply(1:ng, function(i) {
          genoi <- geno[[i]]
          wi <- weight[row.names(genoi)]
          mafi <- rowMeans(genoi[, ]) / 2
          scalei <- 1 / (2 * sum(mafi * (1 - mafi)))
          if (is.null(wi)) {
            DZ <- (genoi[, ref] - mafi)
          }else {
            wi <- length(wi) * wi / sum(wi)
            DZ <- (genoi[, ref] - mafi) * sqrt(wi)
          }
          Ki <- Klist[[i]]
          K_inv <- try(solve(Ki[, ] + diag(1, ncol(Ki)) * (1e-10)), silent = TRUE)[ref, ref]
          if (inherits(K_inv, "try-error")) {
            warning("The GRM of geno component ", i, " is singular!\nUsing general inverse insteaded")
            K_inv <- MASS::ginv(Ki[, ])[ref, ref]
          }
          ef <- colSums(tcrossprod(K_inv, DZ) * ebv.g[, i])
          Effect <- scalei * ef
          pch <-try(chol(P), silent = TRUE)
          if (inherits(pch, "try-error")) {
            warning("The ", i, "th P matrix is not positive definite\nTry spectral decomposition instead of cholesky decomposition")
            eig <- eigen(P, symmetric = TRUE)
            vec <- eig$vectors
            val <- eig$values
            DZS <- tcrossprod(DZ, t(vec))
            diagv <- apply(DZS, 1, function(i) {
              return(sum(val * i^2))
            })
            SE <- vg[i] * scalei * sqrt(diagv)
            rm(ef, eig, vec, val, DZS, pch, diagv)
          }else {
            ZP <- tcrossprod(DZ, pch)
            n <- dim(pch)[2]
            diagv <- (n - 1) * Rfast::rowVars(ZP) + n * rowMeans(ZP) ^ 2
            SE <- vg[i] * scalei * sqrt(diagv)
            rm(ef, ZP, pch, diagv)
          }
          Pvalue <- 2 * stats::pt(abs(Effect / SE), DF, lower.tail = FALSE)
          if (sum(is.na(Pvalue)) != 0) {
            warning("Non-positive variance appeared at ", which(is.na(Pvalue)), " markers in geno group ", i)
            Pvalue[which(is.na(Pvalue))] <- 1
          }
          rm(genoi, wi, mafi, Ki, scalei, DZ, K_inv)
          return(cbind(Effect, SE, Pvalue))
        }, mc.cores = ncpus, mc.preschedule = FALSE)
      }else {
        suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
        suppressMessages(snowfall::sfLibrary(MASS))
        suppressMessages(snowfall::sfLibrary(stats))
        suppressMessages(snowfall::sfExport("geno", "weight", "Klist",
                                            "ebv.g", "vg", "P", "DF"))
        gwas <- snowfall::sfLapply(1:ng, function(i) {
          genoi <- geno[[i]]
          wi <- weight[row.names(genoi)]
          mafi <- rowMeans(genoi[, ]) / 2
          scalei <- 1 / (2 * sum(mafi * (1 - mafi)))
          if (is.null(wi)) {
            DZ <- (genoi[, ref] - mafi)
          }else {
            wi <- length(wi) * wi / sum(wi)
            DZ <- (genoi[, ref] - mafi) * sqrt(wi)
          }
          Ki <- Klist[[i]]
          K_inv <- try(solve(Ki[, ] + diag(1, ncol(Ki)) * (1e-10)), silent = TRUE)[ref, ref]
          if (inherits(K_inv, "try-error")) {
            warning("The GRM of geno component ", i, " is singular!\nUsing general inverse insteaded")
            K_inv <- MASS::ginv(Ki[, ])[ref, ref]
          }
          ef <- colSums(tcrossprod(K_inv, DZ) * ebv.g[, i])
          Effect <- scalei * ef
          pch <-try(chol(P), silent = TRUE)
          if (inherits(pch, "try-error")) {
            warning("The ", i, "th P matrix is not positive definite\nTry spectral decomposition instead of cholesky decomposition")
            eig <- eigen(P, symmetric = TRUE)
            vec <- eig$vectors
            val <- eig$values
            DZS <- tcrossprod(DZ, t(vec))
            diagv <- apply(DZS, 1, function(i) {
              return(sum(val * i^2))
            })
            SE <- vg[i] * scalei * sqrt(diagv)
            rm(ef, eig, vec, val, DZS, pch, diagv)
          }else {
            ZP <- tcrossprod(DZ, pch)
            n <- dim(pch)[2]
            diagv <- (n - 1) * Rfast::rowVars(ZP) + n * rowMeans(ZP) ^ 2
            SE <- vg[i] * scalei * sqrt(diagv)
            rm(ef, ZP, pch, diagv)
          }
          Pvalue <- 2 * stats::pt(abs(Effect / SE), DF, lower.tail = FALSE)
          if (sum(is.na(Pvalue)) != 0) {
            warning("Non-positive variance appeared at ", which(is.na(Pvalue)), " markers in geno group ", i)
            Pvalue[which(is.na(Pvalue))] <- 1
          }
          rm(genoi, wi, mafi, Ki, scalei, DZ, K_inv)
          return(cbind(Effect, SE, Pvalue))
        })
        suppressMessages(snowfall::sfStop())
      }
    }
    gwas <- do.call(rbind, gwas)
    colnames(gwas) <- c("Effect", "SE", "Pvalue")
    rname <- lapply(1:ng, function(i) {
      return(row.names(geno[[i]]))
    })
    rownames(gwas) <- unlist(rname)
    rm(rname)
    GWAS <- gwas
    if (verbose) {
      pbapply::setTimerProgressBar(pb, 1)
    }
  }
  revl <- list(GP = GP, GWAS = GWAS)
  return(revl)
}
