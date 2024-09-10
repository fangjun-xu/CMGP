#' Transform the time from second format
#'
#' @param x numeric, elapsed time in second
#'
#' @return a character
#' @export
#'
#' @examples
#' #time.trans(100)
time.trans <- function(x) {
  h <- x %/% 3600
  m <- (x %% 3600) %/% 60
  s <- ((x %% 3600) %% 60)
  index <- which(c(h, m, s) != 0)
  num <- c(h, m, s)[index]
  char <- c("h", "m", "s")[index]
  return(paste(round(num), char, sep = "", collapse = ""))
}
