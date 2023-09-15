#' Family Objects for Prior Distributions in WALS
#'
#' \code{familyPrior} objects provide a convenient way to specify the prior
#' distribution usedfor the Bayesian posterior mean estimation of the WALS
#' estimators in \link{wals}, \link{walsGLM} and \link{walsNB}.
#'
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution
#' (\code{\link[WALS]{ddgengamma}}).
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution
#' (\code{\link[WALS]{ddgengamma}}).
#' @param object Object of of class \code{familyPrior}.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' \code{familyPrior} is a generic function that extracts the family used in
#' \code{wals} objects.
#'
#' The double (reflected) Weibull (\code{\link[WALS]{ddweibull}}), Subbotin
#' (\code{\link[WALS]{dsubbotin}}) and Laplace distributions
#' (\code{\link[WALS]{dlaplace}}) are all special cases of the reflected
#' generalized gamma distribution (\code{\link[WALS]{ddgengamma}}). The Laplace
#' distribution is also a special case of the double Weibull and of the
#' Subbotin distribution.
#'
#' The default values for the parameters \code{q} and \code{b} are minimax regret
#' solutions for the corresponding priors. The double (reflected) Weibull and
#' Subbotin prior are both neutral and robust. In contrast, the Laplace prior
#' is only neutral but not robust. See section 9 "Enter Bayes: Neutrality and
#' Robustness" of \insertCite{magnus2016wals;textual}{WALS} for details.
#'
#' @returns An object of class \code{familyPrior}. This is a list with the
#' elements
#' \item{q}{Parameter \eqn{q}.}
#' \item{alpha}{Parameter \eqn{\alpha} (of the reflected generalized gamma
#' distribution).}
#' \item{b}{Parameter \eqn{c}.}
#' \item{delta}{Parameter \eqn{\delta = (1 - \alpha)/q}.}
#' \item{printPars}{vector. Contains the parameters that are shown in printing
#' functions, e.g. \code{print.familyPrior}.}.
#' \item{prior}{String with the name of the prior distribution.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [wals], [walsGLM], [walsNB], [computePosterior], [ddweibull],
#' [dsubbotin], [dlaplace], [ddgengamma].
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
#' @returns \code{laplace} returns an object of the specialized class
#' \code{familyPrior_laplace} that inherits from \code{familyPrior}.
#' This allows separate processing of the Laplace prior in the estimation
#' functions as closed-form formulas exists for its posterior mean and variance.
#' The list elements are the same as for objects of class \code{familyPrior}.
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

#' @export
print.familyPrior <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  priorPars <- paste(names(x$printPars), signif(x$printPars, digits),
                     sep = " = ", collapse = ", ")
  cat(paste0("\nPrior: ", x$prior, "(", priorPars, ")\n\n"))
  invisible(x)
}

## Densities of distributions --------------------------------------------------

#' Double (reflected) weibull density
#'
#' Wrapper around \link[VGAM]{dgengamma.stacy} of \link[VGAM]{VGAM} to use the
#' parametrization on pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param x vector of quantiles.
#' @inheritParams familyPrior
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @returns Gives the (log-)density.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [weibull], [ddgengamma].
#'
#' @export
ddweibull <- function(x, q, b, log = FALSE) {
  ddgengamma(x, q = q, alpha = 1 - q, b = b, log = log)
}


#' Subbotin density
#'
#' Wrapper around \link[VGAM]{dgengamma.stacy} of \link[VGAM]{VGAM} to use the
#' parametrization on pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @inheritParams ddweibull
#'
#' @returns Gives the (log-)density.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [subbotin], [ddgengamma].
#'
#' @export
dsubbotin <- function(x, q, b, log = FALSE) {
  ddgengamma(x, q = q, alpha = 0.0, b = b, log = log)
}


#' Laplace density
#'
#' Wrapper around \link[VGAM]{dgengamma.stacy} of \link[VGAM]{VGAM} to use the
#' parametrization on pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @inheritParams ddweibull
#'
#' @returns Gives the (log-)density.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [laplace], [ddgengamma].
#'
#' @export
dlaplace <- function(x, b, log = FALSE) {
  ddgengamma(x, q = 1.0, alpha = 0.0, b = b, log = log)
}


# TODO: consider replacing dgengamma.stacy of VGAM with generalized normal
# because that's what double generalized gamma is but with restrictions on
# its parameters. gnorm packages provides it.

#' Double generalized gamma density
#'
#' Wrapper around \code{\link[VGAM]{dgengamma.stacy}} of \link[VGAM]{VGAM} to
#' use the parametrization on pp. 131 of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param x vector of quantiles.
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param alpha \eqn{\alpha} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#' The density function is
#' \deqn{\pi(x) = \frac{q c^{(1 - \alpha)/q}}{2 \Gamma((1 - \alpha)/q)} |x|^{-\alpha} \exp(-c |x|^{q}).}
#'
#' The function uses \code{\link[VGAM]{dgengamma.stacy}} internally from
#' \link[VGAM]{VGAM} that uses the parametrization in table 12.13, p.369 of
#' \insertCite{yee2015vgam;textual}{WALS}.
#'
#' @returns Gives the (log-)density.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
ddgengamma <- function(x, q, alpha, b, log = FALSE) {
  d <- q
  scale <- b^(-(1/d))
  k <- (1 - alpha)/q # delta

  if (log) {
    return(VGAM::dgengamma.stacy(abs(x), scale = scale, d = d, k = k, log = TRUE) - log(2.0))
  } else {
    return(VGAM::dgengamma.stacy(abs(x), scale = scale, d = d, k = k, log = FALSE) / 2.0)
  }
}
