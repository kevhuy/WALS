#' Compute posterior mean and variance of normal location problem
#'
#' Computes the posterior mean and variance of the normal location problem with
#' fixed variance to 1, i.e. \eqn{x | \gamma \sim N(\gamma, 1)}.
#' The priors for \eqn{\gamma} are either \link[WALS]{weibull},
#' \link[WALS]{subbotin} or \link[WALS]{laplace}. Their properties are briefly
#' discussed in \insertCite{magnus2016wals;textual}{WALS}.
#' Currently only the \link[WALS]{weibull} and \link[WALS]{laplace} prior for
#' \eqn{\gamma} are supported.
#' Default method of computePosterior uses numerical integration. This is used
#' for the \link[WALS]{weibull}, \link[WALS]{subbotin} priors in my case.
#' For the \link[WALS]{laplace} prior closed form expressions exist for the integrals.
#' In the original MATLAB code, the Gauss-Kronrod quadrature was used for
#' numerical integration. Here we use the default \code{\link[stats]{integrate}} which
#' combines Gauss-Kronrod with Wynn's Epsilon algorithm for extrapolation.
#'
#' @param x Observed values, i.e. in WALS these are the regression coefficients
#' of the transformed regressor Z2 standardized by the standard deviation
#' \eqn{\gamma_{2u} / s}.
#' @param prior Object of class \code{familyPrior}, e.g. \code{\link[WALS]{weibull}},
#' should contain all necessary parameters needed for the posterior.
#'
#'
#' @details
#' See section "Numerical integration in Bayesian estimation step"
#' in the appendix of \insertCite{huynhwals;textual}{WALS} for details.
#'
#'
#' \code{computePosterior.familyPrior_laplace} is the specialized method for the
#' S3 class \code{familyPrior_laplace} and computes the posterior first and
#' second moments of the normal location problem with a
#' Laplace prior using the analytical formula (without numerical integration).
#' For more details, see \insertCite{deluca2020laplace;textual}{WALS} and the
#' original code of Magnus and De Luca.
#'
#' @references
#' \insertAllCited{}
#'
#' @references Original MATLAB code on Jan Magnus' website.
#' \url{https://www.janmagnus.nl/items/WALS.pdf}
#'
computePosterior <- function(x, ...) UseMethod("computePosterior", x)


#' @rdname computePosterior
computePosterior.default <- function(prior, x) {
  # Use numerical integration for Weibull and Subbotin priors.

  # preallocate
  postMean <- rep(NA, length(x))
  postVariance <- rep(NA, length(x))

  # Numerical integration and posterior mean and var. computations
  for (i in 1:length(x)) {
    A0 <- integrate(A0f, lower = 0, upper = Inf, x = x[i],
                    priorDensity = prior$density)
    if (A0$message != "OK") warning(paste0("Warning in integrate A0f: ", A0$message))

    A1 <- integrate(A1f, lower = 0, upper = Inf, x = x[i],
                    priorDensity = prior$density)
    if (A1$message != "OK") warning(paste0("Warning in integrate A1f: ", A1$message))

    A2 <- integrate(Ajf, lower = 0, upper = Inf, x = x[i],
                    priorDensity = prior$density, j = 2.0)
    if (A2$message != "OK") warning(paste0("Warning in integrate A2f: ", A2$message))

    postMean[i] <- x[i] - A1$value/A0$value
    postVariance[i] <- A2$value/A0$value - (A1$value/A0$value)^2.0
  }

  return(list(postMean = postMean, postVariance = postVariance))
}

#' Compute posterior moments for Laplace prior
#'
#' @rdname computePosterior
computePosterior.familyPrior_laplace <- function(prior, x) {
  signx <- sign(x)
  absx <- abs(x)

  g0 <- pnorm(-absx - prior$b)
  g1 <- pnorm(absx - prior$b)
  g2 <- dnorm(absx - prior$b)
  # psi0 <- g0 / g1
  psi1 <- g2 / g1
  # psi2 <- exp(2.0 * prior$b * absx) * psi0
  psi2 <- exp(2.0 * prior$b * absx + log(g0) - log(g1))
  hratio <- (1.0 - psi2) / (1.0 + psi2)

  postMean <- signx * (absx - prior$b * hratio)
  postVariance <- (1.0 + (prior$b^2.0) * (1.0 - (hratio^2.0))
                   - prior$b * (1.0 + hratio) * psi1)

  return(list(postMean = postMean, postVariance = postVariance))
}


# helper functions for computePosterior, e.g. integrands for A0, A1 and A2

# Integrand of A0, eq.25 in Magnus and De Luca (2016, pp. 134), WALS survey
A0f <- function(gam, x, priorDensity) {
  (dnorm(x - gam) + dnorm(x + gam)) * priorDensity(gam)
}

# Integrand of A1
A1f <- function(gam, x, priorDensity) {
  ( (x - gam)*dnorm(x - gam) + (x + gam)*dnorm(x + gam) ) * priorDensity(gam)
}


# Integrand of Aj with general power j
Ajf <- function(gam, x, priorDensity, j) {
  ( ((x - gam)^j)*dnorm(x - gam) + ((x + gam)^j)*dnorm(x + gam) ) * priorDensity(gam)
}
