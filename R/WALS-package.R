## Specify global imports
#' @importFrom Rdpack reprompt
#' @import stats
#' @import methods
NULL

.onLoad <- function(lib, pkg){
  Rdpack::Rdpack_bibstyles(package = pkg, authors = "LongNames")
  invisible(NULL)
}
