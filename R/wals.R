#' Weighted Average Least Squares for linear regression models
#'
#' R port of MATLAB code wals.m (version 2.0, revision 18 december 2013)
#' by J.R. Magnus and G.De Luca, available from https://www.janmagnus.nl/items/WALS.pdf.
#' Calculates WALS estimates and variances when some regressors (X1) are present
#' in all models and model selection takes place over the rest (X2).
#'
#' @references
#' \insertRef{deluca2011stata}{WALS}
#'
#' \insertRef{kumar2013normallocation}{WALS}
#'
#' \insertRef{magnus2010growth}{WALS}
#'
#' \insertRef{magnus2016wals}{WALS}
#'
#' @export
wals <- function(x, ...) UseMethod("wals", x)

#' @export
wals.formula <- function(formula, data, subset, na.action, weights, offset,
                         family, prior = weibull(), model = TRUE, keepY = TRUE,
                         keepX = FALSE, sigma = NULL, ...) {
  ## call
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] < 2L) {
    # TODO: Implement what happens when only one part is specified in formula
    # Either interpret all as focus regressors (same as ML estimator basically)
    # or include all as auxiliary regressors and keep only constant as
    # focus regressor.

    stop("One part formula not implemented yet")
    # formula <- as.Formula(formula(formula), ~ 1)
    # simpleFormula <- TRUE
  } else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simpleFormula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  # extract terms, model matrix, response
  mm <- extractModel(formula, mf, data)

  # extract objects from mm
  Y <- mm$Y; X1 <- mm$X1; X2 <- mm$X2; mt <- mm$mt; mtX1 <- mm$mtX1; mtX2 <- mm$mtX2
  cont <- mm$cont
  n <- length(Y)

  rm(mm) # save RAM, do not need it anymore

  # check if X1 and X2 contain the same variables
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")

  ## weights (not used yet)
  weights <- processWeights(weights, mf, n)

  ## offsets (not used yet)
  offset <- getOffset(formula, mf, cl, n)

  ## Fit model
  out <- wals.fit(X1, X2, Y, sigma, prior, ...)

  # add more elements
  if (keepY) out$y <- Y
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$call <- cl
  out$formula <- oformula
  out$terms <- list(focus = mtX1, aux = mtX2, full = mt)
  out$levels <- list(focus = .getXlevels(mtX1, mf), aux = .getXlevels(mtX2, mf),
                     full = .getXlevels(mt, mf))
  out$contrasts <- cont
  if (model) out$model <- mf

  class(out) <- "wals"
  return(out)
}

#' Fitter function for Weighted Average Least Squares estimation
#'
#' Workhorse function behind \link{wals} and \link{walsGLM}.
#'
#' @usage
#' wals.fit(
#'  X1,
#'  X2,
#'  y,
#'  sigma = NULL,
#'  prior = weibull(),
#'  method = "original",
#'  svdTol = .Machine$double.eps,
#'  svdRtol = 1e-6,
#'  keepUn = FALSE,
#'  eigenSVD = TRUE,
#'  prescale = TRUE,
#'  ...
#'  )
#'
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1's) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param y Count response as vector
#' @param sigma if NULL (default), then the variance of the error term is estimated,
#' see p.136 of \insertCite{magnus2016wals;textual}{WALS}. If sigma is specified,
#' then the unrestricted estimator is divided by sigma before performing the
#' Bayesian posterior mean estimation.
#' @param prior Object of class \code{familyPrior}. For example \link[WALS]{weibull}
#' or \link[WALS]{laplace}. Not tested with other priors.
#' @param method Specifies method used. Available methods are
#' \code{"original"} (default) or \code{"svd"}.
#' @param svdTol Tolerance for rank of matrix \eqn{\bar{Z}_{1}}
#' Only used if \code{method = "svd"}.
#' Checks if smallest eigenvalue in SVD of \eqn{\bar{Z}_1} and \eqn{\bar{Z}}
#' is larger than \code{svdTol}, otherwise reports a rank deficiency.
#' @param svdRtol Relative tolerance for rank of matrix \eqn{\bar{Z}_{1}}.
#' Only used if \code{method = "svd"}. Checks if ratio of largest to smallest
#' eigenvalue in SVD of \eqn{\bar{Z}_1} is larger than  \code{svdRtol},
#' otherwise reports a rank deficiency.
#' @param keepUn If \code{TRUE}, keeps the estimators of the unrestricted model,
#' i.e. \eqn{\tilde{\gamma}_{u}}.
#' @param eigenSVD If \code{TRUE}, then \code{semiorthogonalize()} uses \code{svd()}
#' to compute the eigendecomposition of \eqn{\bar{\Xi}} instead of \code{eigen()}.
#' In this case, the tolerances of \code{svdTol} and \code{svdRtol} are used to
#' determine whether \eqn{\bar{\Xi}} is of full rank (need it for \eqn{\bar{\Xi}^{-1/2}}).
#' @param prescale If \code{TRUE} (default), prescales the regressors X1 and X2 with
#' \eqn{\Delta_1} and \eqn{\Delta_2}, respectively, to improve numerical stability
#' and make the coefficients of the auxiliary regressors scale equivariant.
#' See \insertCite{deluca2011stata;textual}{WALS} for more details.
#' \strong{WARNING: It is not recommended to set \code{prescale = FALSE}.}
#' The option \code{prescale = FALSE} only exists for historical reasons.
#' @param ... Arguments for internal function \code{computePosterior()}.
#'
#' @references
#' \insertAllCited{}
#'
#' @export wals.fit
wals.fit <- function(X1, X2, y, sigma = NULL, prior = weibull(),
                     method = "original", svdTol = .Machine$double.eps,
                     svdRtol = 1e-6, keepUn = FALSE, eigenSVD = TRUE,
                     prescale = TRUE, ...) {
  ### TODO: Implement using SVD but also take advantage of orthogonality of
  ### U and V. Currently only use it for inverting (X1bar'X1bar)
  ## Sanity checks

  # number of obs.
  n <- length(y)
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  stopifnot(n == n1, n == n2, n1 == n2)

  # number of focus and auxiliary regressors
  k1 <- ncol(X1)
  k2 <- ncol(X2)
  k <- k1 + k2
  stopifnot(k <= n)

  # store names of regressors for later
  X1names <- colnames(X1)
  X2names <- colnames(X2)
  Xnames <- c(X1names, X2names)


  ## Step 2.a: Scaling X1 so all diag elements of (X1*Delta1)'X1*Delta1 = 1 ----
  ## Corresponds to equation (8) of De Luca et al. 2018 except division by n
  ## missing

  Delta1 <- if (prescale) 1.0 / sqrt(colSums(X1^2.0)) else rep(1.0, ncol(X1))

  # multiply each row with Delta1
  Z1 <- multAllRows(X1, Delta1)
  # Idea: could we define Z1 by standardization of X1 by its standard deviation
  # instead? Could use scale(X1, ...) then.


  ## Step 2.b: Scaling X2 so diag els. of (X2*Delta2)'M1*X2*Delta2 = 1 ---------

  if (method == "original") {

    Z1inv <- solve(crossprod(Z1, Z1))

    # empty entries for svd
    svdZ1 <- NULL

  } else if (method == "svd") {
    svdZ1 <- svd(Z1) # only need D and V of the SVD, but need U later
    singularVals <- svdZ1$d

    checkSingularitySVD(singularVals, tol = svdTol, rtol = svdRtol)

    Z1inv <- svdZ1$v %*% ( (1.0 / ((singularVals)^2.0)) * t(svdZ1$v) )

    ## TODO:
    ## can do this more elegantly by applying SVD to Z1 too, not only for
    ## inverting but only in following operations below with V12 etc.
    ## can use the u of svdZ1...


  } else stop(paste("method", method, " not implemented"))

  VV12 <- crossprod(Z1, X2)
  ZZ <- crossprod(VV12, Z1inv %*% VV12)
  Z2d <- crossprod(X2, X2) - ZZ

  Delta2 <- if (prescale) 1.0 / sqrt(diag(Z2d)) else rep(1.0, ncol(Z2d))
  Z2s <- multAllRows(Delta2*Z2d, Delta2)


  ## Step 3: Semi-orthogonalization of Z2s -------------------------------------

  outSemiOrt <- semiorthogonalize(Z2s, X2, Delta2, eigenSVD)


  ## Step 4: OLS of unrestricted model -----------------------------------------
  Z <- cbind(Z1, outSemiOrt$Z2)
  lmUnrestricted <- lm.fit(Z, y, singular.ok = FALSE)

  gammaUnrestricted2 <- lmUnrestricted$coefficients[(k1 + 1):k]

  # estimate sigma if unknown
  if (is.null(sigma)) {
    sigma <- sqrt( sum(lmUnrestricted$residuals^2.0) / (n - k) )
  }

  x <- gammaUnrestricted2 / sigma

  ## Step 5: Compute posterior mean and variance of x ~ N(\gamma, 1) -----------

  outPosterior <- computePosterior(prior, x, ...)

  ## Step 6 & 7: WALS estimates and precision ----------------------------------
  walsEstimates <- gammaToBeta(outPosterior, y, Z1, outSemiOrt$Z2, Delta1,
                               outSemiOrt$D2, sigma, Z1inv, method = method,
                               svdZ1 = svdZ1)
  walsEstimates$sigma <- sigma
  walsEstimates$prior <- prior
  walsEstimates$method <- method

  if (keepUn) {
    walsEstimates$gammaUn1 <- lmUnrestricted$coefficients[1:k1]
    walsEstimates$gammaUn2 <- gammaUnrestricted2
    names(walsEstimates$gammaUn2) <- X2names
  }

  # store names of focus and auxiliary regressors
  walsEstimates$X1names <- X1names
  walsEstimates$X2names <- X2names

  # reassign names to variables
  names(walsEstimates$coef) <- Xnames
  names(walsEstimates$beta2) <- X2names
  names(walsEstimates$gamma2) <- X2names
  colnames(walsEstimates$vcovBeta) <- Xnames
  row.names(walsEstimates$vcovBeta) <- Xnames
  colnames(walsEstimates$vcovGamma) <- Xnames
  row.names(walsEstimates$vcovGamma) <- Xnames

  # class(walsEstimates) <- "wals"
  return(walsEstimates)
}

#' @rdname wals
#'
#' @param y response vector.
#' @param X1 matrix of focus regressors.
#' @param X2 matrix of auxiliary regressors.
#'
#' @examples
#' X <- model.matrix(mpg ~ disp + hp + wt + vs + am + carb, data = mtcars)
#' X1 <- X[,c("(Intercept)", "disp", "hp", "wt")] # focus
#' X2 <- X[,c("vs", "am", "carb")] # auxiliary
#' y <- mtcars$mpg
#'
#' wals(y, X1, X2, prior = weibull())
#'
#' @export
wals.default <- function(y, X1, X2, ...) {

  out <- wals.fit(X1 = X1, X2 = X2, y = y, ...)
  class(out) <- "wals"
  return(out)
}
