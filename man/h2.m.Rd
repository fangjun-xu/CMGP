% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h2.m.R
\name{h2.m}
\alias{h2.m}
\title{Calculate the PVE}
\usage{
h2.m(
  y = NULL,
  CV = NULL,
  geno = NULL,
  map = NULL,
  random = NULL,
  alpha = 0,
  ldscore = NULL,
  bin = 1e+06,
  EMsteps = 0,
  EMsteps_fail = 10,
  EM_alpha = 1,
  max_iter = 50,
  eps = 0.01,
  ncpus = 1,
  verbose = TRUE
)
}
\arguments{
\item{y}{phenotype vector, NAs refer to validation}

\item{CV}{covariates}

\item{geno}{Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size}

\item{map}{geno information, SNP, Chr, Pos}

\item{random}{list, non-additive genetic random effect}

\item{alpha}{genetype scaling value}

\item{ldscore}{LD score}

\item{bin}{region size for correlation calculation}

\item{EMsteps}{see ?gaston::lmm.aireml}

\item{EMsteps_fail}{see ?gaston::lmm.aireml}

\item{EM_alpha}{see ?gaston::lmm.aireml}

\item{max_iter}{see ?gaston::lmm.aireml}

\item{eps}{see ?gaston::lmm.aireml}

\item{ncpus}{integer value, default as 1, how many threads used}

\item{verbose}{logical, TRUE for print the output}
}
\value{
matrix, Beta, SE, Pvalue
}
\description{
Calculate the PVE
}
\examples{
#\donttest{
#revl <- h2.m(y = y, CV = CV, geno = geno, map = map,
#random = random, alpha = alpha, ldscore = ldscore,
#bin = bin, EMsteps = EMsteps, EM_alpha = EM_alpha,
#EMsteps_fail = EMsteps_fail, max_iter = max_iter,
#eps = eps, ncpus = ncpus, verbose = verbose)
#}
}
