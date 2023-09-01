predictCounts <- function(x, ...) UseMethod("predictCounts", x)

predictCounts.familyWALScount <- function(family, yUnique, rowNames, eta, ...) {
  # ... to insert other parameters, e.g. for negbin also need scale
  # Inspired by code in countreg (Kleiber & Zeileis)
  nUnique <- length(yUnique)
  rval <- matrix(NA, nrow = length(eta), ncol = nUnique)
  dimnames(rval) <- list(rowNames, yUnique)

  for (i in 1:nUnique) rval[,i] <- family$density(yUnique[i], eta = eta, ...)

  return(rval)
}

predictCounts.family <- function(x, ...) {
  stop("Probability predictions of counts not supported for ",
       sQuote(x$family), " family.")
}

predictCounts.default <- function(x, ...) {
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}
