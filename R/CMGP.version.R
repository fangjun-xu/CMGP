#' CMGP.version: Print the version of CMGP
#'
#' @param verbose logical, whether print
#'
#' @return NULL
#'
#' @examples
#' #CMGP.version(TRUE)
CMGP.version <-  function(verbose = TRUE) {
  if (verbose) {
    cat(paste("+",
              paste(rep("=", 19), collapse = ""),
              "< Welcome to CMGP >",
              paste(rep("=", 19), collapse = ""),
              "+", sep = ""), "\n")
    cat("|      A Comprehensive Mixed-linear-model-based Tool      |\n")
    cat("|                 for Genomic Prediction                  |\n")
    cat("|          _ _ _    _      _     _ _ _    _ _ _           |\n")
    cat("|         / _ _ \\  |  \\  /  |   / _ _ \\  |  _ _ \\         |\n")
    cat("|        / /   \\_| | \\ \\/ / |  / /   \\_| | /   \\ |        |\n")
    cat("|       | |        | |\\  /| | | |        | \\_ _/ |        |\n")
    cat("|       | |        | | \\/ | | | |   _ _  |  _ _ /         |\n")
    cat("|       | |     _  | |    | | | |  |__ | | |              |\n")
    cat("|        \\ \\_ _/ | | |    | |  \\ \\_ _/ | | |              |\n")
    cat("|         \\_ _ _/  |_|    |_|   \\_ _ _/  |_|     V1.0.0   |\n")
    cat("|                                                         |\n")
    cat("|    E-mail: fjxu@webmail.hzau.edu.cn                     |\n")
    cat(paste("+", paste(rep("=", 57), collapse = ""), "+", sep = ""), "\n")
  }
}
