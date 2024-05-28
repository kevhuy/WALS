#' Internal methods: Predict probability for counts
#'
#' Predicts the probability of counts given a family object of class
#' \code{"\link{familyWALScount}"}. Only works for count data models.
#'
#' @param x object of class \code{"\link[WALS]{familyWALScount}"}.
#'
#' @details
#' \code{"\link{familyWALScount}"} objects are used in the fitting methods
#' \code{\link[WALS]{walsNB}}, \code{\link[WALS]{walsNBmatrix}},
#' \code{\link[WALS]{walsGLM}} or \code{\link[WALS]{walsGLMmatrix}}. For the
#' latter two, only the family \code{\link[WALS]{poissonWALS}} is currently
#' supported.
#'
#' \code{predictCounts()} is not available for objects of any class except for
#' \code{"\link{familyWALScount}"}.
#'
#' @references
#' \insertAllCited{}
#'
predictCounts <- function(x, ...) UseMethod("predictCounts", x)

#' @rdname predictCounts
#'
#' @param yUnique vector. The counts (larger or equal to zero) which to predict
#' probabilities for.
#' @param rowNames vector. The names of the observations.
#' @param eta vector. The fitted linear link \eqn{\hat{\eta}} of the model.
#' @param ... Further parameters passed to \code{density()} function in
#' \code{family}.
#'
#' @details
#' The \code{predictCounts.familyWALScount()} method is a modified version of the
#' \code{predict.hurdle()} method from the \code{countreg} package
#' version 0.2-1 (2023-06-13) \insertCite{countreg,countreghurdle}{WALS} using the argument
#' \code{type = "prob"}.
#'
#' @returns Returns a matrix of dimension \code{length(eta)} times
#' \code{length{yUnique}} with the predicted probabilities of the counts given
#' in \code{yUnique} for every observation in \code{eta}.
#'
#' @export
predictCounts.familyWALScount <- function(x, yUnique, rowNames, eta, ...) {
  # ... to insert other parameters, e.g. for negbin also need scale
  # Inspired by code in countreg (Kleiber & Zeileis)
  nUnique <- length(yUnique)
  rval <- matrix(NA, nrow = length(eta), ncol = nUnique)
  dimnames(rval) <- list(rowNames, yUnique)

  for (i in 1:nUnique) rval[,i] <- x$density(yUnique[i], eta = eta, ...)

  return(rval)
}

#' @export
predictCounts.family <- function(x, ...) {
  stop("Probability predictions of counts not supported for ",
       sQuote(x$family), " family.")
}

#' @export
predictCounts.default <- function(x, ...) {
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}
