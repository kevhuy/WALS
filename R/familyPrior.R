## Contains distribution families for priors in WALS

#' Double (reflected) Weibull prior
#'
#' Double (reflected) Weibull prior for Bayesian posterior mean estimation of
#' WALS estimators. Special case of the reflected generalized gamma distribution.
#'
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#'
#' @details
#' The default values for the parameters \code{q} and \code{b} are minimax regret
#' solutions. The double weibull prior is both neutral and robust.
#' See section 9 "Enter Bayes: Neutrality and Robustness" of
#' \insertCite{magnus2016wals;textual}{WALS} for details.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @export
weibull <- function(q = 0.887630085544086, b = log(2.0)) {
  alpha <- 1 - q
  delta <- (1 - alpha)/q

  # fix parameters of distribution
  dens <- function(x, log = FALSE){
    ddweibull(x, q = q, b = b, log = log)
  }

  out <- list(q = q, alpha = alpha, b = b, delta = delta,
              density = dens, prior = "weibull")
  class(out) <- "familyPrior"
  return(out)
}


#' Subbotin prior
#'
#' Subbotin prior for Bayesian posterior mean estimation of WALS estimators.
#' Special case of the reflected generalized gamma distribution.
#'
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#'
#' @details
#' The default values for the parameters \code{q} and \code{b} are minimax regret
#' solutions. The Subbotin prior is both neutral and robust.
#' See section 9 "Enter Bayes: Neutrality and Robustness" of
#' \insertCite{magnus2016wals;textual}{WALS} for details.
#'
#' @references
#' \insertAllCited{}
#'
#'
#'
#' @export
subbotin <- function(q = 0.799512530172489, b=0.937673273794677) {
  alpha <- 0.0
  delta <- (1.0 - alpha)/q

  # fix parameters of distribution
  dens <- function(x, log = FALSE) {
    ddsubbotin(x, q = q, b = b, log = log)
  }

  out <- list(q = q, alpha = alpha, b = b, delta = delta,
              density = dens, prior = "subbotin")
  class(out) <- "familyPrior"
  return(out)
}


#' Laplace prior
#'
#' Laplace prior for Bayesian posterior mean estimation of WALS estimators.
#' Special case of the reflected generalized gamma distribution and of
#' the Weibull and Subbotin prior.
#'
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#'
#' @details
#' The default value for the parameter \code{b} is the minimax regret
#' solution. The Laplace prior is only neutral but not robust.
#' See section 9 "Enter Bayes: Neutrality and Robustness" of
#' \insertCite{magnus2016wals;textual}{WALS} for details.
#'
#' @references
#' \insertAllCited{}
#'
#'
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

  out <- list(q = q, alpha = alpha, b = b, delta = delta,
              density = dens, prior = "laplace")
  class(out) <- c("familyPrior_laplace", "familyPrior")
  return(out)
}


#' Double (reflected) weibull density
#'
#' Wrapper around dgengamma.stacy() of VGAM to use the parametrization on pp. 131
#' of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param x vector of quantiles.
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @references
#' \insertAllCited{}
#'
#' @export
ddweibull <- function(x, q, b, log = FALSE) {
  ddgengamma(x, q = q, alpha = 1 - q, b = b, log = log)
}


#' Subbotin density
#'
#' Wrapper around dgengamma.stacy() of VGAM to use the parametrization on pp. 131
#' of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param q \eqn{q} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @references
#' \insertAllCited{}
#'
#' @export
ddsubbotin <- function(x, q, b, log = FALSE) {
  ddgengamma(x, q = q, alpha = 0.0, b = b, log = log)
}


#' Laplace density
#'
#' Wrapper around dgengamma.stacy() of VGAM to use the parametrization on pp. 131
#' of \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @param x vector of quantiles.
#' @param b \eqn{c} in \insertCite{magnus2016wals;textual}{WALS}.
#' Parameter of reflected generalized gamma distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @references
#' \insertAllCited{}
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
#' Wrapper around dgengamma.stacy() of VGAM to use the parametrization on pp. 131
#' of \insertCite{magnus2016wals;textual}{WALS}.
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
#' \deqn{\pi(x) = \frac{q c^{(1 - \alpha)/q}}{2 \Gamma((1 - \alpha)/q)} |x|^{-\alpha} \exp(-c |x|^{q})}.
#'
#' The function uses \code{dgengamma.stacy} internally from \code{VGAM} that uses
#' the parametrization in table 12.13, p.369 of \insertCite{yee2015vgam;textual}{WALS}.
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
