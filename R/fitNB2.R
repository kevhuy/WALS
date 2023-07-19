#' Fits a NB2 regression via maximum likelihood with log link for mean and
#' dispersion parameter.
#' 
#' Only works with log-link so far, no other links tested.
#'
#'
#' @param X Design matrix.
#' @param Y Count response vector.
#' @param family Object of class \code{\link[stats]{family}}.
#' @param control List of parameters for controlling the optimization process.
#' Use \link[WALS]{controlNB} to generate the list.
#'
#' @details
#' The available parameters for controlling the optimization are documented in
#' \link[WALS]{controlNB}.
#'
#'
fitNB2 <- function(X, Y, family, control = controlNB()) {

  # ignore starting values in control, use MASS::glm.nb
  if (control$initMASS) { # reduce memory, limit output by model = FALSE, x = FALSE ...
    fit <- MASS::glm.nb(Y ~ -1 + X, model = FALSE, x = FALSE, y = FALSE,
                        control = glm.control(maxit = control$controlOptim$maxit,
                                              epsilon = control$epsilonMASS))
    out <- list(coefficients = coef(fit), theta = fit$theta,
                # generate error code 99 if IWLS algo in glm.nb failed
                convergence = if(fit$converged) 0 else 99,
                ll = logLik(fit), message = NULL, start = NULL)
    return(out)
  }

  ## starting values
  startBeta <- control$start$mu
  startLogTheta <- control$start$logTheta

  ## default starting values

  # use poisson regression for regression coefs
  if (is.null(startBeta)) startBeta <- glm.fit(X, Y, family = poisson(family$link))$coefficient

  # maximize loglik wrt theta given regression coefs
  if (is.null(startLogTheta)) {
    if (control$initThetaMASS) {
      mu <- family$linkinv(X %*% startBeta)
      startLogTheta <- log(MASS::theta.ml(Y, mu, n = length(Y),
                                          limit = control$controlOptim$maxit,
                                          eps = control$eps))
    } else {
      startLogTheta <- 0
    }
  }

  start <- c(startBeta, startLogTheta)


  ## loglik and gradient
  nllfun <- function(parms, k) {
    mu <- family$linkinv(X %*% parms[1:k])
    theta <- exp(parms[k + 1]) # optimize wrt log theta --> unconstrained

    return(-sum(dnbinom(Y, size = theta, mu = mu, log = TRUE))) # negative loglik
  }

  ngrad <- function(parms, k) { # negative gradient because minimize
    eta <- X %*% parms[1:k]
    mu <- family$linkinv(eta)
    theta <- exp(parms[k + 1])
    gr <- countreg::snbinom(Y, mu = mu, size = theta)
    gr <- cbind(gr[, 1] * family$mu.eta(eta)[,,drop = TRUE] * X, gr[, 2] * theta)
    return(-colSums(gr))
  }

  k <- ncol(X)

  ## ML estimation
  fit <- optim(fn = nllfun, gr = ngrad, par = start, k = k,
               method = control$method, control = control$controlOptim)
  out <- list(coefficients = fit$par[1:k], theta = exp(fit$par[k + 1]),
              convergence = fit$convergence, ll = -fit$value, message = fit$message,
              start = start)
  return(out)
}
