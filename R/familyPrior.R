#' Family Objects for Prior Distributions in WALS
#'
#' \code{"familyPrior"} objects provide a convenient way to specify the prior
#' distribution used for the Bayesian posterior mean estimation of the WALS
#' estimators in \code{\link[WALS]{wals}}, \code{\link[WALS]{walsGLM}} and
#' \code{\link[WALS]{walsNB}}
#'
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution. See below for details.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution. See below for details.
#' @param object,x Object of of class \code{"familyPrior"} or \code{"\link[WALS]{wals}"}.
#' The function \code{familyPrior()} accesses the \code{"familyPrior"} objects
#' that are stored in objects of class \code{"\link[WALS]{wals}"}.
#' @param digits The number of significant digits to display.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' \code{familyPrior()} is a generic function that extracts the family used in
#' \code{"\link[WALS]{wals}"} objects.
#'
#' The density function of the reflected generalized gamma distribution is
#' \deqn{\pi(x) = \frac{q c^{(1 - \alpha)/q}}{2 \Gamma((1 - \alpha)/q)}
#'                |x|^{-\alpha} \exp(-c |x|^{q}).}
#'
#'
#' The double (reflected) Weibull, Subbotin and Laplace distributions are all
#' special cases of the reflected generalized gamma distribution. The Laplace
#' distribution is also a special case of the double Weibull and of the Subbotin
#' distribution.
#'
#' The double (reflected) Weibull density sets \eqn{q = 1 - \alpha}, the Subbotin
#' density sets \eqn{\alpha = 0} and the Laplace density sets \eqn{\alpha = 0}
#' and \eqn{q = 1}.
#'
#' The default values for the parameters \code{q} and \code{b} are minimax regret
#' solutions for the corresponding priors. The double (reflected) Weibull and
#' Subbotin prior are both neutral and robust. In contrast, the Laplace prior
#' is only neutral but not robust. See section 9 "Enter Bayes: Neutrality and
#' Robustness" of \insertCite{magnus2016wals;textual}{WALS} for details and
#' Table 1 for the optimal parameter values.
#'
#' @returns An object of class \code{"familyPrior"} to be used in
#' \code{\link[WALS]{wals}}, \code{\link[WALS]{walsGLM}} and \code{\link[WALS]{walsNB}}.
#' This is a list with the elements
#' \item{q}{Parameter \eqn{q}.}
#' \item{alpha}{Parameter \eqn{\alpha} (of the reflected generalized gamma
#' distribution).}
#' \item{b}{Parameter \eqn{c}.}
#' \item{delta}{Parameter \eqn{\delta = (1 - \alpha)/q}.}
#' \item{printPars}{vector. Contains the parameters that are shown in printing
#' functions, e.g. \code{print.familyPrior()}.}
#' \item{prior}{String with the name of the prior distribution.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [wals], [walsGLM], [walsNB], [computePosterior], [ddweibull],
#' [dsubbotin], [dlaplace].
#'
#' @examples
#' ## Use in wals():
#' fit <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
#'             law + tropics + avelf + confucian, data = GrowthMPP,
#'             prior = weibull(q = 0.8, b = log(1.8)))
#' summary(fit)
#'
#' @export
familyPrior <- function(object, ...) UseMethod("familyPrior", object)


#' Double (reflected) Weibull prior
#' @rdname familyPrior
#' @export
weibull <- function(q = 0.887630085544086, b = log(2.0)) {
  alpha <- 1 - q
  delta <- (1 - alpha)/q

  # fix parameters of distribution
  dens <- function(x, log = FALSE){
    ddweibull(x, q = q, b = b, log = log)
  }

  out <- list(q = q, alpha = alpha, b = b, delta = delta, density = dens,
              printPars = c(q = q, b = b), prior = "weibull")
  class(out) <- "familyPrior"
  return(out)
}


#' Subbotin prior
#' @rdname familyPrior
#' @export
subbotin <- function(q = 0.799512530172489, b = 0.937673273794677) {
  alpha <- 0.0
  delta <- (1.0 - alpha)/q

  # fix parameters of distribution
  dens <- function(x, log = FALSE) {
    dsubbotin(x, q = q, b = b, log = log)
  }

  out <- list(q = q, alpha = alpha, b = b, delta = delta, density = dens,
              printPars = c(q = q, b = b), prior = "subbotin")
  class(out) <- "familyPrior"
  return(out)
}


#' Laplace prior
#' @rdname familyPrior
#' @aliases familyPrior_laplace
#' @returns \code{laplace()} returns an object of the specialized class
#' \code{"familyPrior_laplace"} that inherits from \code{"familyPrior"}.
#' This allows separate processing of the Laplace prior in the estimation
#' functions as closed-form formulas exists for its posterior mean and variance.
#' The list elements are the same as for objects of class \code{"familyPrior"}.
#'
#' @export
laplace <- function(b = log(2.0)) {
  alpha <- 0.0
  q <- 1.0
  delta <- (1.0 - alpha)/q

  # fix parameters of distribution
  dens <- function(x, log = FALSE) {
    dlaplace(x, b = b, log = log)
  }

  out <- list(q = q, alpha = alpha, b = b, delta = delta, density = dens,
              printPars = c(b = b), prior = "laplace")
  class(out) <- c("familyPrior_laplace", "familyPrior")
  return(out)
}

## Class methods ---------------------------------------------------------------
#' @rdname familyPrior
#' @export
print.familyPrior <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  priorPars <- paste(names(x$printPars), signif(x$printPars, digits),
                     sep = " = ", collapse = ", ")
  cat(paste0("\nPrior: ", x$prior, "(", priorPars, ")\n\n"))
  invisible(x)
}

## Densities of prior distributions --------------------------------------------

#' Internal function: double (reflected) Weibull density
#'
#' Wrapper around \code{\link[stats]{dweibull}} to use the parametrization on
#' pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param x vector of quantiles.
#' @inheritParams familyPrior
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @returns Gives the (log-)density.
#'
#' @details
#' The density function is
#' \deqn{\pi(x) = \frac{q c}{2} |x|^{q - 1} \exp(-c |x|^{q}).}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [weibull], [dweibull].
#'
ddweibull <- function(x, q, b, log = FALSE) {
  shape <- q
  scale <- exp(-log(b)/shape)

  if (log) {
    return(dweibull(abs(x), shape = shape, scale = scale, log = TRUE) - log(2.0))
  } else {
    return(dweibull(abs(x), shape = shape, scale = scale, log = FALSE) / 2.0)
  }
}


#' Internal function: Subbotin density
#'
#' Subbotin density, uses the parametrization on pp. 131 of
#' \insertCite{magnus2016wals;textual}{WALS}. Also called generalized normal
#' distribution.
#'
#' @inheritParams ddweibull
#'
#' @returns Gives the (log-)density.
#'
#' @details
#' The density function is
#' \deqn{\pi(x) = \frac{q c^{1/q}}{2 \Gamma(1/q)} \exp(-c |x|^{q}).}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [subbotin].
#'
dsubbotin <- function(x, q, b, log = FALSE) {
  delta <- 1.0 / q
  loglik <- -(b*(abs(x)^(q))) + delta * log(b) - log(delta) - log(2.0) - lgamma(delta)
  if (log) return(loglik) else return(exp(loglik))
}


#' Internal function: Laplace density
#'
#' Wrapper around \code{\link[WALS]{dsubbotin}} with fixed \code{q = 1}. Uses
#' the parametrization on pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @inheritParams ddweibull
#'
#' @details
#' The density function is
#' \deqn{\pi(x) = \frac{c}{2} \exp(-c |x|).}
#'
#' @returns Gives the (log-)density.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [laplace], [dsubbotin].
#'
dlaplace <- function(x, b, log = FALSE) {
  return(dsubbotin(x, q = 1, b = b, log = log))
}
