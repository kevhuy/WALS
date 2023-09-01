#' Poisson family with additional functions
#'
#' Extends \code{poisson} family from \code{stats} with transformations required
#' for \code{walsGLM}.
#'
#' @param ... see \link[stats]{family} for details.
#' @export
poissonWALS <- function(...) {
  # copy default family
  fam <- poisson(...)


  if (fam$link == "log") {
    # Canonical link

    # d theta/ d eta as function of linear link eta
    fam$theta.eta <- function(etaStart) 1.0
    fam$psi <- function(etaStart) fam$variance(fam$linkinv(etaStart))

    # could insert specialized versions of transformX and transformY here
    # if they are numerically better

  } else stop(sprintf("%s link not implemented yet", fam$link))

  fam$initializeY <- function(y) {
    # only checks y
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    return(y)
  }

  ## functions for Xbar_1, Xbar_2 and ybar
  # does nothing with y, just passes it
  fam$transformX <- function(X, etaStart, y = NULL) transformX(X, etaStart,
                                                               fam$psi(etaStart))
  fam$transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart) {
    y <- fam$initializeY(y)
    transformY(y, X1bar, X2bar, beta1, beta2, etaStart,
               vBar = fam$theta.eta(etaStart), muBar = fam$linkinv(etaStart),
               psiBar = fam$psi(etaStart))
  }

  fam$density <- function(x, eta, log = FALSE, ...) {
    mu <- fam$linkinv(eta)
    return(dpois(x, lambda = mu, log = log))
  }

  class(fam) <- c("familyWALScount", "familyWALS", class(fam))
  return(fam)
}

#' Binomial family with additional functions
#'
#' Extends \code{binomial} family from \code{stats} with transformations required
#' for \code{walsGLM}.
#'
#' @param ... see \link[stats]{family} for details.
#' @export
binomialWALS <- function(...) {
  fam <- binomial(...)

  if (fam$link == "logit") {
    # Canonical link

    fam$theta.eta <- function(etaStart) 1.0
    fam$psi <- function(etaStart) fam$variance(fam$linkinv(etaStart))

  } else stop(sprintf("%s link not implemented yet", fam$link))

  fam$initializeY <- function(y) {
    # convert factor to numeric
    if (is.factor(y)) y <- y != levels(y)[1L]
    return(y)
  }

  ## functions for Xbar_1, Xbar_2 and ybar
  # does nothing with y, just passes it
  fam$transformX <- function(X, etaStart, y = NULL) transformX(X, etaStart,
                                                               fam$psi(etaStart))
  fam$transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart) {
    y <- fam$initializeY(y)
    transformY(y, X1bar, X2bar, beta1, beta2, etaStart,
               vBar = fam$theta.eta(etaStart), muBar = fam$linkinv(etaStart),
               psiBar = fam$psi(etaStart))
  }

  fam$density <- function(x, eta, log = FALSE, ...) {
    mu <- fam$linkinv(eta)
    return(dbinom(x, size = 1, prob = mu, log = log, ...))
  }

  class(fam) <- c("familyWALS", class(fam))
  return(fam)
}


#' Negative binomial family
#'
#' Reconstruct family object for negative binomial type 2 (NB2) with fixed
#' scale parameter theta. Analogous to \link[MASS]{negative.binomial} in
#' \code{MASS} but \code{MASS} uses non-canonical link.
#'
#' @param theta dispersion parameter of NB2, always larger than 0.
#' @param link specifies link function, currently only "log" and "canonical"
#' are supported.
#'
#' @export
negativeBinomial <- function(theta, link = "log") {
  if (link == "canonical") {
    # implementing functions required in glm.fit
    linkfun <- function(mu) log(mu) - log(mu + theta)
    variance <- function(mu) mu*(1 + mu/theta)
    linkinv <- function(eta) theta/(exp(-eta) - 1)
    mu.eta <- function(eta) theta * (exp(-eta) / ((exp(-eta) - 1)^2))
    ## alternatively (maybe more numerically stable)
    # mu.eta <- function(eta) theta / (exp(eta) - 2 + exp(-eta))
    valideta <- function(eta) TRUE
    validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)

    ## copied from MASS
    aic <- function(y, n, mu, wt, dev) {
      term <- ( (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1)
                - theta * log(theta) + lgamma(theta) - lgamma(theta + y) )
      return(2 * sum(term * wt))
    }

    dev.resids <- function(y, mu, wt) {
      return(2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta))))
    }

    # initializes mustart and n for IRLS in glm.fit
    initialize <- expression({
      if (any(y < 0))
        stop("negative values not allowed for the negative binomial family")
      n <- rep(1, nobs)
      mustart <- y + (y == 0)/6
    })


    famname <- paste("Negative Binomial(", format(round(theta, 4)), ")",
                     sep = "")
    return(structure(list(family = famname, link = link, linkfun = linkfun,
                   linkinv = linkinv, variance = variance,
                   dev.resids = dev.resids, aic = aic, mu.eta = mu.eta,
                   initialize = initialize, validmu = validmu,
                   valideta = valideta, simulate = NULL),
              class = "family"))

  } else {
    return(MASS::negative.binomial(theta = theta, link = link))
  }
}



#' NB2 family with fixed dispersion parameter and additional functions
#'
#' Family object for negative binomial distribution type 2 (NB2) with fixed
#' dispersion parameter. Extends \link{negativeBinomial} with functions
#' required in \code{walsGLM}.
#'
#' @param scale dispersion parameter of NB2, always larger than 0.
#' @param link specifies link function, currently only "log" and "canonical"
#' are supported.
#'
#' @export
negbinFixedWALS <- function(scale, link) {
  fam <- negativeBinomial(theta = scale, link = link)

  fam$initializeY <- function(y) {
    # only checks for valid values
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    return(y)
  }

  if (fam$link == "canonical") {
    # Canonical link

    # d theta/ d eta as function of linear link eta
    fam$theta.eta <- function(etaStart) 1.0
    fam$psi <- function(etaStart) scale * (exp(-etaStart) / ((exp(-etaStart) - 1)^2))

    ## alternatively (maybe more numerically stable)
    # fam$psi <- function(etaStart) scale / (exp(etaStart) - 2 + exp(-etaStart))

    fam$transformX <- function(X, etaStart, y = NULL) {
      # does nothing with y, just passes it
      # as.vector to allow broadcasting along columns
      return(as.vector( (sqrt(scale)*exp(0.5*etaStart)) / (1 - exp(etaStart)) ) * X)
    }
    fam$transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart) {
      y <- fam$initializeY(y)
      return( X1bar %*% beta1 + X2bar %*% beta2
              + (exp(-0.5)/sqrt(scale)) * (y*(1 - exp(etaStart)) - scale*exp(etaStart)) )
    }
  } else if (fam$link == "log") {

    # d theta / d eta as function of linear link eta and y
    fam$theta.eta <- function(eta) scale / (exp(eta) + scale)
    fam$psi <- function(eta, y) {
      return(exp( eta + log(scale) + ( log(y + scale) - 2.0*log(exp(eta) + scale) ) ))
    }

    fam$transformX <- function(X, eta, y) {
      return(as.vector(  exp(0.5*eta) * sqrt(scale * (scale + y))/(exp(eta) + scale)) * X)
    }
    fam$transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart) {
      y <- fam$initializeY(y)
      return( X1bar %*% beta1 + X2bar %*% beta2
              + exp(-0.5*etaStart) * sqrt( scale / (scale + y))
              * (y - exp(etaStart)) )
    }

  } else stop(sprintf("%s link not implemented yet", fam$link))


  fam$density <- function(x, eta, log = FALSE, ...) {
    mu <- fam$linkinv(eta)
    return(dnbinom(x, size = scale, mu = mu, log = log))
  }

  class(fam) <- c("familyWALScount", "familyWALS", class(fam))
  return(fam)
}

#' NB2 family with variable scale parameter and additional functions
#'
#' Family object for negative binomial distribution type 2 (NB2) with variable
#' dispersion parameter. Extends \link{negativeBinomial} with functions
#' required in \link{walsNB}.
#'
#' @param scale dispersion parameter of NB2 to be used.
#' @param link specifies the link function, currently only "log" is supported.
#'
#' @export
negbinWALS <- function(scale, link) {
  alpha <- log(scale)
  out <- negbinFixedWALS(scale, link)

  if (link == "log") {
    out$q <- function(eta, y) { # q_i = C'(y-muStart)
      return(exp(eta - 2*log(exp(eta) + exp(alpha))) * (y - exp(eta)))
    }

    out$g <- function() return(scale)

    out$transformX <- function(X, eta, y) { # transform X to Xbar
      return(as.numeric( (exp(0.5*(eta + alpha)) * sqrt(y + scale)) /
        (exp(eta) + scale) ) * X)
    }

    out$transformY0 <- function(y, X1bar, X2bar, beta1, beta2, eta) {
      y <- out$initializeY(y)
      mu <- exp(eta)
      uplus <- ((y - mu) * ( exp(0.5*(alpha - eta - log(y + scale))) *
        ((mu*(1.0 - alpha) + scale)/(mu + scale))  ))

      return(X1bar %*% beta1 + X2bar %*% beta2 + uplus)
    }

    out$t <- function(eta, y) {
      n <- length(y)
      mu <- exp(eta)
      sum1scale <- sum((y - mu) / (exp(2.0*eta - alpha) + 2.0*exp(eta) + exp(alpha)))
      sum2scale <- sum(mu / (exp(eta - alpha) + 1))
      g2ksum <- (sum1scale + sum2scale +
                   scale * (sum(trigamma(y + scale)) - n*trigamma(scale)))

      return(scale * ( out$kappaSum(eta, y)*(1.0 - alpha) - sum(out$q(eta,y)*eta)) -
               alpha*(g2ksum))
    }

    out$epsilon <- function(eta, y) {
      return(1.0 / out$epsiloninv(eta, y))
    }

    out$epsiloninv <- function(eta, y) {
      n <- length(y)
      mu <- exp(eta)
      sum1scale <- sum((y - mu) / (exp(2.0*eta - alpha) + 2.0*exp(eta) + exp(alpha)))
      sum2 <- sum((2.0*mu - y) / (mu + scale))

      return(sum1scale + sum2 + (sum(digamma(y + scale)) - n*digamma(scale)) -
        sum(log(mu + scale)) +
        scale*(sum(trigamma(y + scale)) - n*trigamma(scale)) + n * alpha)
    }

    out$kappaSum <- function(eta, y){
      n <- length(y)
      mu <- exp(eta)
      return(sum((mu - y)/(mu + scale)) - (sum(log(mu + scale)) - n*alpha) +
        (sum(digamma(y + scale)) - n * digamma(scale)))
    }

    out$psi <- function(eta, y) {
      return(exp( eta + alpha + ( log(y + scale) - 2.0*log(exp(eta) + scale) ) ))
    }

    out$computeAlpha <- function(gamma1, gamma2, Z1, Z2, y, eta, q) {
      kappaSum <- out$kappaSum(eta, y)
      sum2 <- out$epsilon(eta, y) * sum(q * (Z1 %*% gamma1 + Z2 %*% gamma2))

      return((sum(q*eta) - kappaSum) / out$epsiloninv(eta, y) + alpha - sum2)
    }

  } else {
    stop("link not implemented")
  }


  class(out) <- c("familyNBWALS", class(out))
  return(out)
}
