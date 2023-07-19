predictCounts <- function(x, ...) UseMethod("predictCounts", x)

predictCounts.default <- function(family, yUnique, rowNames, eta, ...) {
  # ... to insert other parameters, e.g. for negbin also need scale
  nUnique <- length(yUnique)
  rval <- matrix(NA, nrow = length(eta), ncol = nUnique)
  dimnames(rval) <- list(rowNames, yUnique)

  for (i in 1:nUnique) rval[,i] <- family$density(yUnique[i], eta = eta, ...)

  return(rval)
}
