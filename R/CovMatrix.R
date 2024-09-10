#' Transform the covariates as model matrix
#' @description
#' The quantitative cov. were inputed as numeric format and remained unchanged,
#' the discrete cov. were transformed as one hot model matrix
#'
#' @param CV a matrix or data.frame for covariates
#' @param n integer, number of individuals
#'
#' @return a model matrix
#' @export
#'
#' @examples
#' cv1 <- rnorm(100)
#' cv2 <- rep(c("m", "f"), 50)
#' n <- 100
#' CV <- CovMatrix(cbind(cv1,cv2), 100)
CovMatrix <- function(CV = NULL, n = NULL) {
  Cov <- matrix(1, n)
  if (!is.null(CV)) {
    if (sum(is.na(CV)) != 0) stop("'NA' isn't allowed in CV!")
    for (i in 1:ncol(CV)){
      if (is.numeric(CV[, i])) {
        #quantitative
        Cov <- cbind(Cov, CV[, i])
      }else {
        #discrete
        if (length(unique(CV[, i])) == 1) {
          stop("No groups in column ", i, " of CV!")
        }
        ci <- nnet::class.ind(CV[, i])
        Cov <- cbind(Cov, ci[, -1])
      }
    }
  }
  return(Cov)
}
