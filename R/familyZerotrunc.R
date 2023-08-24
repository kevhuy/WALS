#' @export
zerotruncPoissonWALS <- function(link = "log") {
  fam <- list(untrunc = poisson(link), zerotrunc = NULL)
  
  # needed for predictions of type = "zero" and "count"
  fam$untrunc$distFun <- function(x, eta, ...) {
    mu <- fam$untrunc$linkinv(eta)
    ppois(x, lambda = mu, ...)
  }
  
  fam$zerotrunc$density <- function(x, eta, log = FALSE, ...) {
    # eta is converted to mu inside distFun and also here...redundant...
    # can we remove one redundancy?
    mu <- fam$untrunc$linkinv(eta)
    # logOneminusp0 <- fam$untrunc$distFun(0, eta = eta, lower.tail = FALSE,
    #                                      log.p = TRUE)
    logOneminusp0 <- ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)
    logDens <- dpois(x, lambda = mu, log = TRUE) - logOneminusp0
    if (log) {
      return(logDens)
    } else {
      return(exp(logDens))
    }
    
  }
  
  if (link == "log") {
    # Canonical link
    fam$zerotrunc$theta.eta <- function(etaStart) 1.0
    fam$zerotrunc$psi <- function(etaStart) {
      thetaStart <- etaStart
      muStart <- fam$untrunc$linkinv(etaStart)
      
      exp(thetaStart + ppois(1.0, lambda = muStart, lower.tail = FALSE, log.p = TRUE)
          - 2.0*ppois(0, lambda = muStart, lower.tail = FALSE, log.p = TRUE))
    }
    
    # get mean of zerotruncated distribution
    fam$zerotrunc$linkinv <- function(eta) {
      theta <- eta
      mu <- fam$untrunc$linkinv(eta) # mean of untruncated
      pmax(exp(theta - ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE)),
           .Machine$double.eps)
    }
    
    fam$zerotrunc$variance <- function(mu) {
      eta <- log(mu)
      fam$zerotrunc$psi(eta) # for canonical link, psi = sigma^2
    }
    
    # functions for zerotruncBetaStart
    # negative loglik to min.
    fam$zerotrunc$nllfun <- function(parms, X, Y, linkinv, weights, offset) {
      # does not do anything with linkinv, just passes it
      # needed later for more general links
      mu <- as.vector(exp(X %*% parms + offset))
      -sum(weights * (
        dpois(Y, lambda = mu, log = TRUE) - ppois(0, lambda = mu, lower.tail = FALSE,
                                                  log.p = TRUE)
      ))
    }
    
    # negative gradient corresponding to negative loglik
    fam$zerotrunc$ngrad <- function(parms, X, Y, linkinv, weights, offset) {
      # does not do anything with linkinv
      eta <- as.vector(X %*% parms + offset)
      mu <- exp(eta)
      -colSums(((Y - mu) - exp(ppois(0, lambda = mu, log.p = TRUE) -
              ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)) * weights * X)
    }
    
    # transformations for WALS
    fam$zerotrunc$transformX <- function(X, etaStart, y = NULL) {
      transformX(X, etaStart, fam$zerotrunc$psi(etaStart))
    } 
    
    fam$zerotrunc$transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart) {
      thetaStart <- etaStart
      muStart <- fam$untrunc$linkinv(etaStart)
      muZeroStart <- fam$zerotrunc$linkinv(etaStart)
      
      psisqrtInv <- exp(-0.5*(thetaStart 
                        - ppois(1.0, lambda = muStart, lower.tail = FALSE,
                                    log.p = TRUE))
                        + ppois(0, lambda = muStart, lower.tail = FALSE, log.p = TRUE))
      
      X1bar %*% beta1 + X2bar %*% beta2 + psisqrtInv * (y - muZeroStart)
    }

  } else stop(sprintf("%s link not implemented yet", link))
  
  class(fam) <- c("familyWALSzerotruncPoisson","familyWALSzerotrunc")
  return(fam)
}
