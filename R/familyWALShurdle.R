#' Constructor of familyWALShurdle
#' 
#' Uses a zerotruncated family and a family for the zero-part to construct
#' a joint density function used in predict.walsHurdle
#' @export
familyWALShurdle <- function(count, zero) {
  out <- list(count = count, zero = zero)
  
  out$density <- function(x, etaCount, etaZero, log = FALSE, ...) {
    # init vector
    x1 <- x > 0
    rval <- as.numeric(rep(NA, length(x)))
    
    zerotruncCount <- count$zerotrunc$density(x[x1], etaCount[x1], log = TRUE)
    probLargerZero <- zero$density(1L, etaZero, log = TRUE, size = 1L)
    rval[x1] <- probLargerZero[x1] + zerotruncCount
    rval[!x1] <- log(1 - exp(probLargerZero[!x1]))
    
    if (log) return(rval) else return(exp(rval))
  }
  
  class(out) <- "familyWALShurdle"
  return(out)
}
