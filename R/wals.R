#' Weighted-Average Least Squares for linear regression models
#'
#' Performs model averaging for linear regression models using the
#' Weighted-Average Least Squares method by
#' \insertCite{magnus2010growth;textual}{WALS}. See also
#' \insertCite{deluca2011stata;textual}{WALS},
#' \insertCite{kumar2013normallocation;textual}{WALS} and
#' \insertCite{magnus2016wals;textual}{WALS}.
#'
#' @details
#' R port of MATLAB code wals.m (version 2.0, revision 18 December 2013)
#' by J.R. Magnus and G. De Luca, available from
#' \url{https://www.janmagnus.nl/items/WALS.pdf}.
#' Calculates WALS estimates when focus regressors (X1) are present in all
#' submodels and model averaging takes place over the auxiliary regressors (X2).
#'
#' @references
#' \insertAllCited{}
#'
#' @export
wals <- function(x, ...) UseMethod("wals", x)

#' \code{wals.formula()} uses formulas to specify the design matrix.
#' @rdname wals
#'
#' @param formula an object of class \code{"\link[Formula]{Formula}"}
#' (or one that can be coerced to that class, e.g. \code{"\link[stats]{formula}"}):
#' a symbolic description of the model to be fitted.
#' The details of model specification are given under ‘Details’.
#' @param data an optional data frame, list or environment
#' (or object coercible by \code{\link[base]{as.data.frame}} to a data frame)
#' containing the variables in the model. If not found in \code{data}, the variables
#' are taken from \code{environment(formula)}, typically the environment which
#' the function is called from.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights **not implemented yet.**
#' @param offset **not implemented yet.**
#' @param na.action **not implemented yet.**
#' @param prior Object of class \code{"\link[WALS]{familyPrior}"}. For example
#' \code{\link[WALS]{weibull}} or \code{\link[WALS]{laplace}}.
#' @param model if \code{TRUE} (default), then the model.frame is stored in
#' the return.
#' @param keepY if \code{TRUE} (default), then the response is stored in
#' the return.
#' @param keepX if \code{TRUE}, then the model matrices are stored in the return.
#' the return.
#' @param sigma if NULL (default), then the variance of the error term is
#' estimated. See \code{\link[WALS]{walsFit}} for more details.
#' @param ... Arguments for workhorse \code{\link[WALS]{walsFit}}.
#'
#' @details
#' Formulas typically contain two parts, i.e. they are of the form
#' "\code{y ~ X11 + X12 | X21 + X22}", where the variables before "\code{|}" are
#' the focus regressors (includes a constant by default) and the ones after
#' "\code{|}" are the auxiliary regressors. If only a one-part formula is
#' specified, then all regressors are considered as auxiliary regressors and only
#' a constant is employed as focus regressor, i.e.
#' "\code{y ~ X1 + X2}" is equivalent to "\code{y ~ 1 | X1 + X2}".
#'
#' **WARNING:** Interactions in formula do not work properly yet.
#' It is recommended to manually create the interactions beforehand and then
#' to insert them as 'linear terms' in the formula.
#'
#' @returns \code{wals.formula()} returns an object of class
#' \code{"wals"}. This is a list that contains all elements returned from
#' \code{\link[WALS]{walsFit}} and additionally
#' \item{y}{If \code{keepY = TRUE}, contains the response vector.}
#' \item{x}{list. If \code{keepX = TRUE}, then it is a list with elements
#' \code{x1} and \code{x2} containing the design matrices of the focus and
#' auxiliary regressors, respectively.}
#' \item{weights}{returns the argument \code{weights}.}
#' \item{offset}{returns the argument \code{offset}.}
#' \item{cl}{Call of the function.}
#' \item{formula}{\code{formula} used.}
#' \item{terms}{List containing the model terms of the focus and auxiliary
#' regressors separately, as well as for the full model.}
#' \item{levels}{List containing the levels of the focus and auxiliary
#' regressors separately, as well as for the full model.}
#' \item{contrasts}{List containing the contrasts of the design matrices of
#' focus and auxiliary regressors.}
#' \item{model}{If \code{model = TRUE}, contains the model frame.}
#'
#' See returns of \code{\link[WALS]{walsFit}} for more details.
#'
#' @examples
#' ## Replicate table on p. 534 of De Luca & Magnus (2011)
#' fitDM <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
#'                   law + tropics + avelf + confucian, data = GrowthMPP,
#'                 prior = laplace())
#' tableDM <- cbind("coef" = coef(fitDM), "se" = sqrt(diag(vcov(fitDM))))
#' print(round(tableDM, 7))
#'
#'
#' ## Replicate first panel of Table I in Amini & Parmeter (2012)
#' data("datafls", package = "BMS")
#'
#' # NOTE: Authors manually scale data, then rescale the resulting coefs and se.
#' X <- model.matrix(y ~ ., data = datafls)
#' Xscaled <- apply(X, MARGIN = 2, function(x) x/max(x))
#' Xscaled <- Xscaled[,-1]
#' scaleVector <- apply(X, MARGIN = 2, function(x) max(x))
#' flsScaled <- as.data.frame(cbind(y = datafls$y, Xscaled))
#'
#' # NOTE: prescale = FALSE, still used old version of WALS in Magnus et al. (2010).
#' # Not recommended anymore!
#' fitFLS <- wals(y ~ 1 | ., data = flsScaled, prescale = FALSE, eigenSVD = FALSE,
#'                prior = laplace())
#' tableFLS <- cbind('coef' = coef(fitFLS)/scaleVector,
#'                   'se' = sqrt(diag(vcov(fitFLS)))/scaleVector)
#' printVars <- c("(Intercept)", "GDP60", "Confucian", "LifeExp", "EquipInv",
#'                "SubSahara", "Muslim", "RuleofLaw")
#' print(round(tableFLS[printVars,], 4))
#'
#'
#' ## Replicate third panel of Table I in Amini & Parmeter (2012)
#' data("SDM", package = "BayesVarSel")
#'
#' # rescale response
#' SDM$y <- SDM$y / 100
#'
#' # NOTE: Authors manually scale data, then rescale the resulting coefs and se.
#' X <- model.matrix(y ~ ., data = SDM)
#' Xscaled <- apply(X, MARGIN = 2, function(x) x/max(x))
#' Xscaled <- Xscaled[,-1]
#' scaleVector <- apply(X, MARGIN = 2, function(x) max(x))
#' SDMscaled <- as.data.frame(cbind(y = SDM$y, Xscaled))
#'
#' # NOTE: prescale = FALSE, still used old version of WALS in Magnus et al. (2010).
#' # Not recommended anymore!
#' fitDW <- wals(y ~ 1 | ., data = SDMscaled, prescale = FALSE, eigenSVD = FALSE,
#'               prior = laplace())
#' tableDW <- cbind(coef(fitDW)/scaleVector, sqrt(diag(vcov(fitDW)))/scaleVector)
#' printVars <- c("(Intercept)", "EAST", "P60", "IPRICE1", "GDPCH60L", "TROPICAR")
#' print(round(tableDW[printVars,], 5))
#'
#' @export
wals.formula <- function(formula, data, subset = NULL, na.action = NULL,
                         weights = NULL, offset = NULL, prior = weibull(),
                         model = TRUE, keepY = TRUE,
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
    # Include all as auxiliary regressors and keep only constant as focus
    # regressor, when one-part formula is specified
    formula <- Formula::as.Formula(~ 1, formula(formula))
  } else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
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
  out <- walsFit(X1, X2, Y, sigma, prior, ...)

  # add more elements
  if (keepY) out$y <- Y
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$weights <- weights
  out$offset <- offset
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

#' \code{wals.matrix()} uses prespecified design matrices x (focus) and
#' x2 (auxiliary) and response vector y.
#' @rdname wals
#'
#' @param x Design matrix of focus regressors. Usually includes a constant
#' (column full of 1s) and can be generated using \code{\link[stats]{model.matrix}}.
#' @param x2 Design matrix of auxiliary regressors. Usually does not include
#' a constant column and can also be generated using \code{\link[stats]{model.matrix}}.
#' @param y Response as vector.
#'
#' @returns \code{wals.matrix()} returns an object of class \code{"walsMatrix"},
#' which inherits from \code{"wals"}. This is a list that contains all elements
#' returned from \code{\link[WALS]{walsFit}} and additionally the response \code{y},
#' the list \code{x} with model matrices \code{x1} and \code{x2}, the call
#' \code{cl}, \code{offset} and \code{weights}.
#'
#' @examples
#' ## Example for wals.matrix()
#' X <- model.matrix(mpg ~ disp + hp + wt + vs + am + carb, data = mtcars)
#' X1 <- X[,c("(Intercept)", "disp", "hp", "wt")] # focus
#' X2 <- X[,c("vs", "am", "carb")] # auxiliary
#' y <- mtcars$mpg
#'
#' wals(X1, X2, y, prior = weibull())
#'
#' @export
wals.matrix <- function(x, x2, y, subset = NULL, na.action = NULL,
                        weights = NULL, offset = NULL, prior = weibull(),
                        keepY = TRUE, keepX = FALSE, sigma = NULL,
                        ...) {
  cl <- match.call()
  X1 <- x
  X2 <- x2
  if (!is.null(subset)) {
    X1 <- X1[subset,]; X2 <- X2[subset,]; y <- y[subset]
  }

  out <- walsFit(X1, X2, y, sigma, prior, ...)

  if (keepY) out$y <- y
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$weights <- weights
  out$offset <- offset
  out$call <- cl

  class(out) <- c("walsMatrix", "wals")
  return(out)
}

#' @rdname wals
#'
#' @details
#' \code{wals.default()} raises an error if \code{x} is not an object of class
#' \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise it calls
#' \code{wals.matrix()}. It is a modified version of \code{glmboost.default}
#' from the \code{mboost} package version 2.9-8 (2023-09-06) \insertCite{mboost}{WALS}.
#'
#' @returns \code{wals.default()} raises an error if \code{x} is not an object
#' of class \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise
#' returns an object of class \code{"walsMatrix"}. See above for more details.
#'
#' @export
wals.default <- function(x, ...) {
  # inspired by glmboost.default in mboost.
  if (extends(class(x), "matrix")) {
    return(wals.matrix(x, ...))
  }
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}


#' Fitter function for Weighted Average Least Squares estimation
#'
#' Workhorse function behind \code{\link[WALS]{wals}} and \code{\link[WALS]{walsGLM}}.
#'
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1s) and can be generated using \code{\link[stats]{model.matrix}}.
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using \code{\link[stats]{model.matrix}}.
#' @param y Response as vector.
#' @param sigma if NULL (default), then the variance of the error term is estimated,
#' see p.136 of \insertCite{magnus2016wals;textual}{WALS}. If sigma is specified,
#' then the unrestricted estimator is divided by sigma before performing the
#' Bayesian posterior mean estimation.
#' @param prior Object of class \code{"\link[WALS]{familyPrior}"}. For example
#' \code{\link[WALS]{weibull}} or \code{\link[WALS]{laplace}}.
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
#' @param eigenSVD If \code{TRUE}, then \code{\link[WALS]{semiorthogonalize}}
#' uses \code{\link[base]{svd}} to compute the eigendecomposition of
#' \eqn{\bar{\Xi}} instead of \code{\link[base]{eigen}}. In this case, the
#' tolerances of \code{svdTol} and \code{svdRtol} are used to
#' determine whether \eqn{\bar{\Xi}} is of full rank (need it for \eqn{\bar{\Xi}^{-1/2}}).
#' @param prescale If \code{TRUE} (default), prescales the regressors X1 and X2 with
#' \eqn{\Delta_1} and \eqn{\Delta_2}, respectively, to improve numerical stability
#' and make the coefficients of the auxiliary regressors scale equivariant.
#' See \insertCite{deluca2011stata;textual}{WALS} for more details.
#' \strong{WARNING: It is not recommended to set \code{prescale = FALSE}.}
#' The option \code{prescale = FALSE} only exists for historical reasons.
#' @param postmult If \code{TRUE}, then it computes
#' \deqn{Z_{2} = X_{2} \Delta_{2} T \Lambda^{-1/2} T^{\top},}
#' where \eqn{T} contains the eigenvectors and \eqn{\Lambda} the eigenvalues
#' from the eigenvalue decomposition
#' \deqn{\Xi = \Delta_2 X_{2}^{\top} M_{1} X_{2} \Delta_2 = T \Lambda T^{\top},}
#' instead of
#' \deqn{Z_{2} = X_{2} \Delta_{2} T \Lambda^{-1/2}.}
#' See \insertCite{huynhwals;textual}{WALS} for more details. The latter is used
#' in the original MATLAB code for WALS in the linear regression model
#' \insertCite{magnus2010growth,deluca2011stata,kumar2013normallocation,magnus2016wals}{WALS},
#' see eq. (12) of \insertCite{magnus2016wals;textual}{WALS}.
#' The first form is required in eq. (9) of \insertCite{deluca2018glm;textual}{WALS}.
#' It is not recommended to set \code{postmult = FALSE} when using \code{\link[WALS]{walsGLM}}
#' and \code{\link[WALS]{walsNB}}.
#' @param ... Arguments for internal function \code{\link[WALS]{computePosterior}}.
#'
#'
#' @returns A list containing
#' \item{coef}{Model averaged estimates of all coefficients.}
#' \item{beta1}{Model averaged estimates of the coefficients of the focus regressors.}
#' \item{beta2}{Model averaged estimates of the coefficients of the auxiliary regressors.}
#' \item{gamma1}{Model averaged estimates of the coefficients of the transformed
#' focus regressors.}
#' \item{gamma2}{Model averaged estimates of the coefficients of the transformed
#' auxiliary regressors.}
#' \item{vcovBeta}{Estimated covariance matrix of the regression coefficients.}
#' \item{vcovGamma}{Estimated covariance matrix of the coefficients of the
#' transformed regressors.}
#' \item{sigma}{Estimated or prespecified standard deviation of the error term.}
#' \item{prior}{\code{familyPrior}. The \code{prior} specified in the arguments.}
#' \item{method}{Stores \code{method} used from the arguments.}
#' \item{betaUn1}{If \code{keepUn = TRUE}, contains the unrestricted
#' estimators of the coefficients of the focus regressors.}
#' \item{betaUn2}{If \code{keepUn = TRUE}, contains the unrestricted
#' estimators of the coefficients of the auxiliary regressors.}
#' \item{gammaUn1}{If \code{keepUn = TRUE}, contains the unrestricted
#' estimators of the coefficients of the transformed focus regressors.}
#' \item{gammaUn2}{If \code{keepUn = TRUE}, contains the unrestricted
#' estimators of the coefficients of the transformed auxiliary regressors.}
#' \item{fitted.values}{Estimated conditional means of the data.}
#' \item{residuals}{Residuals, i.e. response - fitted mean.}
#' \item{X1names}{Names of the focus regressors.}
#' \item{X2names}{Names of the auxiliary regressors.}
#' \item{k1}{Number of focus regressors.}
#' \item{k2}{Number of auxiliary regressors.}
#' \item{n}{Number of observations.}
#' \item{condition}{Condition number of the matrix
#' \eqn{\Xi = \Delta_{2} X_{2}^{\top} M_{1} X_{2} \Delta_{2}}.}
#'
#' @seealso [wals], [walsGLM].
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' X <- model.matrix(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth
#'                   + law + tropics + avelf + confucian, data = GrowthMPP)
#' X1 <- X[, c("(Intercept)", "lgdp60", "equipinv", "school60", "life60", "popgrowth")]
#' X2 <- X[, c("law", "tropics", "avelf", "confucian")]
#' y <- GrowthMPP$gdpgrowth
#'
#' walsFit(X1, X2, y, prior = weibull(), method = "svd")
#'
#' @export
walsFit <- function(X1, X2, y, sigma = NULL, prior = weibull(),
                     method = "original", svdTol = .Machine$double.eps,
                     svdRtol = 1e-6, keepUn = FALSE, eigenSVD = TRUE,
                     prescale = TRUE, postmult = FALSE, ...) {
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

  outSemiOrt <- semiorthogonalize(Z2s, X2, Delta2, eigenSVD, postmult)


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
    walsEstimates$betaUn1 <- as.vector(Delta1 * lmUnrestricted$coefficients[1:k1])
    walsEstimates$betaUn2 <- as.vector(outSemiOrt$D2 %*% gammaUnrestricted2)
    names(walsEstimates$betaUn1) <- X1names
    names(walsEstimates$betaUn2) <- X2names

    walsEstimates$gammaUn1 <- lmUnrestricted$coefficients[1:k1]
    walsEstimates$gammaUn2 <- gammaUnrestricted2
    names(walsEstimates$gammaUn2) <- X2names
  }

  walsEstimates$fitted.values <- drop(X1 %*% walsEstimates$beta1 + X2 %*% walsEstimates$beta2)
  walsEstimates$residuals <- drop(y - walsEstimates$fitted.values)

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

  walsEstimates$k1 <- k1
  walsEstimates$k2 <- k2
  walsEstimates$n <- n
  walsEstimates$condition <- outSemiOrt$condition

  return(walsEstimates)
}


## Class methods ---------------------------------------------------------------
#' Methods for wals and walsMatrix Objects
#'
#' Methods for extracting information from fitted model-averaging objects of
#' classes \code{"wals"} and \code{"walsMatrix"}. \code{"walsMatrix"} objects
#' inherit from \code{"wals"}, so the methods for \code{"wals"} also work for
#' objects of class \code{"walsMatrix"}.
#'
#' @param object,x An object of class \code{"wals"}, \code{"walsMatrix"} or
#' \code{"summary.wals"}.
#' @param newdata Optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param na.action Function determining what should be done with missing values
#' in \code{newdata}. The default is to predict \code{NA}.
#' @param type Character specifying the part of the model that should be returned.
#' For details see below.
#' @param transformed Logical specifying whether the coefficients/covariance
#' matrix of original regressors (\code{FALSE}, default) or the transformed
#' regressors (\code{TRUE}) should be returned.
#' @param digits The number of significant digits to display.
#' @param ... Further arguments passed to methods.
#'
#'
#' @details
#' A set of standard extractor functions for fitted model objects is available
#' for objects of class \code{"wals"} and \code{"walsMatrix"}, including methods to
#' the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}
#' which print the model-averaged estimation of the coefficients along with some
#' further information. As usual, the \code{summary} method returns an object of
#' class \code{"summary.wals"} containing the relevant summary statistics which
#' can then be printed using the associated \code{print} method.
#' Inspired by \insertCite{deluca2011stata;textual}{WALS}, the summary statistics
#' also show \code{Kappa} which is an indicator for the numerical stability of
#' the method, i.e. it shows the square root of the condition number of the
#' matrix \eqn{\Xi = \Delta_{2} X_{2}^{\top} M_{1} X_{2} \Delta_{2}}.
#' The summary further provides information on the prior used along with its
#' parameters. The \code{summary()}, \code{print.summary()},
#' \code{print()} and \code{logLik()} methods are also inspired by the corresponding
#' methods for objects of class \code{"lm"} in \code{\link[stats]{stats}} version
#' 4.3.1 (2023-06-16) \insertCite{R2023}{WALS}, see e.g. \code{\link[stats]{print.summary.lm}}.
#'
#' The \code{\link[stats]{residuals}} method computes raw residuals
#' (observed - fitted).
#'
#' For \code{\link[stats]{coef}} and \code{\link[stats]{vcov}}, the \code{type}
#' argument, either \code{"all"}, \code{"focus"} or \code{"aux"}, specifies which
#' part of the coefficient vector/covariance matrix of the estimates should be
#' returned. Additionally, the \code{transformed} argument specifies whether to
#' return the estimated  coefficients/covariance matrix for the original
#' regressors \eqn{X} or of the transformed regressors \eqn{Z}.
#'
#' The extractors \code{\link[stats]{terms}} and \code{\link[stats]{model.matrix}}
#' behave similarly to \code{coef}, but they only allow \code{type = "focus"}
#' and \code{type = "aux"}. They extract the corresponding component of the model.
#' This is similar to the implementation of these extractors in \code{countreg}
#' version 0.2-1 (2023-06-13) \insertCite{countreg,countreghurdle}{WALS}, see e.g.
#' \code{terms.hurdle()}.
#'
#' @returns \code{predict.wals()} and \code{predict.walsMatrix()} return a vector
#' containing the predicted means.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [wals]
#'
#' @examples
#' ## Example for wals objects
#' fitGrowth <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
#'                 law + tropics + avelf + confucian, data = GrowthMPP,
#'                 prior = laplace())
#' summary(fitGrowth)
#' fitted(fitGrowth)
#' vcov(fitGrowth, type = "aux")
#' familyPrior(fitGrowth)
#' nobs(fitGrowth)
#'
#' ## Example for walsMatrix objects
#' X1 <- model.matrix(fitGrowth, type = "focus")
#' X2 <- model.matrix(fitGrowth, type = "aux")
#' y <- GrowthMPP$gdpgrowth
#' fitGrowthMatrix <- wals(X1, X2, y, prior = laplace())
#' coef(fitGrowthMatrix)
#'
#' @export
predict.wals <- function(object, newdata, na.action = na.pass, ...) {
  # TODO: include offsets
  if (missing(newdata)) {
    return(object$fitted.values)
  } else {
    newMatrices <- genNewdata(object$terms, object$contrasts, newdata,
                              na.action = na.action, xlev = object$levels)
    return(drop(newMatrices$X1 %*% object$beta1 + newMatrices$X2 %*% object$beta2))
  }
}

#' @rdname predict.wals
#' @param newX1 Focus regressors matrix to be used for the prediction.
#' @param newX2 Auxiliary regressors matrix to be used for the prediction.
#' @export
predict.walsMatrix <- function(object, newX1, newX2, ...) {
  # TODO: include offsets
  # Sanity checks
  if (missing(newX1) || missing(newX2)) {
    if (missing(newX1) && missing(newX2)) {
      stop("Missing newX1 and newX2.")
    } else {
      missingArgs <- c("newX1", "newX2")
      stop(paste0("Missing ", missingArgs[c(missing(newX1), missing(newX2))], "."))
    }
  }

  mismatchX1 <- any(colnames(newX1) != object$X1names)
  mismatchX2 <- any(colnames(newX2) != object$X2names)
  if (mismatchX1 || mismatchX2) {
    if (mismatchX1 && mismatchX2) {
      stop("newX1 and newX2 do not contain the same variables as fitted model.")
    } else {
      mismatchArgs <- c("newX1", "newX2")
      stop(paste(mismatchArgs[c(mismatchX1, mismatchX2)],
                 "does not contain the same variables as fitted model."))
    }
  }

  stopifnot(is.matrix(newX1) && is.matrix(newX2))
  stopifnot(nrow(newX1) == nrow(newX2))

  return(drop(newX1 %*% object$beta1 + newX2 %*% object$beta2))
}

#' @rdname predict.wals
#'
#' @returns \code{fitted.wals()} returns a vector containing the fitted means
#' for the data used in fitting.
#'
#' @export
fitted.wals <- function(object, ...) return(object$fitted.values)

#' @rdname predict.wals
#'
#' @returns \code{residuals.wals()} returns the raw residuals of the fitted
#' model, i.e. response - fitted mean.
#'
#' @export
residuals.wals <- function(object, ...) return(object$residuals)

#' @rdname predict.wals
#'
#' @returns \code{print.wals()} invisibly returns its input argument \code{x},
#' i.e. an object of object of class \code{"wals"}.
#'
#' @export
print.wals <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  cat("Focus coefficients: \n")
  print.default(format(x$beta1, digits = digits), print.gap = 2,
                quote = FALSE)

  cat("\nAuxiliary coefficients: \n")
  print.default(format(x$beta2, digits = digits), print.gap = 2,
                quote = FALSE)

  priorPars <- paste(names(x$prior$printPars), signif(x$prior$printPars, digits),
                     sep = " = ", collapse = ", ")
  cat(paste0("\nPrior: ", x$prior$prior, "(", priorPars, ")\n"))

  invisible(x)
}

#' @rdname predict.wals
#'
#' @returns \code{summary.wals} returns an object of class \code{"summary.wals"}
#' which contains the necessary fields for printing the summary in
#' \code{print.summary.wals()}.
#'
#' @export
summary.wals <- function(object, ...) {
  k1 <- object$k1
  k2 <- object$k2
  se <- sqrt(diag(object$vcovBeta))

  object$focusCoefs <- cbind(object$coef[1:k1], se[1:k1])
  object$auxCoefs <- cbind(object$coef[(k1 + 1):(k1 + k2)], se[(k1 + 1):(k1 + k2)])
  colnames(object$focusCoefs) <- colnames(object$auxCoefs) <- c("Estimate", "Std. Error")

  class(object) <- "summary.wals"
  return(object)
}

#' @rdname predict.wals
#'
#' @returns \code{print.summary.wals()} invisibly returns its input argument
#' \code{x}, i.e. an object of object of class \code{"summary.wals"}.
#'
#' @export
print.summary.wals <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")
  printCallCoefs(x, digits, ...)
  printPriorNKappa(x, digits)
  invisible(x)
}

#' @rdname predict.wals
#'
#' @returns \code{coef.wals()} returns a vector containing the fitted coefficients.
#' If \code{type = "focus"}, only the coefficients of the focus regressors are
#' returned and if \code{type = "aux"}, only the coefficients of auxiliary
#' regressors are returned. Else if \code{type = "all"}, the coefficients
#' of both focus and auxiliary regressors are returned. Additionally if
#' \code{transformed = FALSE}, \code{coef.wals()} returns the estimated
#' coefficients for the original regressors \eqn{X} (\eqn{\beta} coefficients)
#' and else if \code{transformed = TRUE} the coefficients of the transformed
#' regressors \eqn{Z} (\eqn{\gamma} coefficients).
#'
#' @export
coef.wals <- function(object, type = c("all", "focus", "aux"),
                      transformed = FALSE, ...) {
  type <- match.arg(type)

  if (transformed) {
    out <- switch(type,
                  "all" = c(object$gamma1, object$gamma2),
                  "focus" = object$gamma1,
                  "aux" = object$gamma2)
    return(out)
  } else {
    out <- switch(type,
                  "all" = object$coef,
                  "focus" = object$beta1,
                  "aux" = object$beta2)
    return(out)
  }
}

#' @rdname predict.wals
#'
#' @returns \code{vcov.wals()} returns a matrix containing the estimated
#' (co-)variances of the fitted regression coefficients. If \code{type = "focus"},
#' only the submatrix belonging to the focus regressors is returned and if
#' \code{type = "aux"}, only the submatrix corresponding to the auxiliary
#' regressors is returned. Else if \code{type = "all"}, the complete covariance
#' matrix is returned. Additionally if \code{transformed = FALSE},
#' \code{vcov.wals()} returns the estimated covariance matrix for the original
#' regressors \eqn{X} (\eqn{\beta} coefficients) and else if
#' \code{transformed = TRUE} the covariance matrix of the transformed regressors
#' \eqn{Z} (\eqn{\gamma} coefficients).
#'
#' @export
vcov.wals <- function(object, type = c("all", "focus", "aux"),
                      transformed = FALSE, ...) {
  type <- match.arg(type)
  k1 <- object$k1; k2 <- object$k2
  if (transformed) {
    out <- switch(type,
                  "all" = object$vcovGamma,
                  "focus" = object$vcovGamma[1:k1, 1:k1],
                  "aux" = object$vcovGamma[(k1 + 1):(k1 + k2), (k1 + 1):(k1 + k2)])
    return(out)
  } else {
    out <- switch(type,
                  "all" = object$vcovBeta,
                  "focus" = object$vcovBeta[1:k1, 1:k1],
                  "aux" = object$vcovBeta[(k1 + 1):(k1 + k2), (k1 + 1):(k1 + k2)])
    return(out)
  }
}

#' @rdname predict.wals
#'
#' @returns \code{nobs.wals()} returns the number of observations used for
#' fitting the model.
#'
#' @export
nobs.wals <- function(object, ...) return(object$n)

#' @rdname predict.wals
#'
#' @returns \code{terms.wals()} returns the *terms* representation of the fitted
#' model. It is of class \code{c("terms", "formula")}, see \code{\link[stats]{terms}}
#' and \code{\link[stats]{terms.object}} for more details. If \code{type = "focus"},
#' then returns the terms for the focus regressors, else if \code{type = "aux"}
#' returns the terms for the auxiliary regressors.
#'
#' @export
terms.wals <- function(x, type = c("focus", "aux"), ...) {
  return(x$terms[[match.arg(type)]])
}

#' @rdname predict.wals
#'
#' @returns \code{model.matrix.wals()} either returns the design matrix of the
#' focus regressors (\code{type = "focus"}) or of the auxiliary regressors
#' (\code{type = "aux"}). See \code{\link[stats]{model.matrix}} for more details.
#'
#' @export
model.matrix.wals <- function(object, type = c("focus", "aux"), ...) {
  type <- match.arg(type)
  if (!is.null(object$x)) {
    return(object$x[[type]])
  } else if (!is.null(object$model)) {
    out <- model.matrix(object$terms[[type]], object$model, contrasts = object$contrasts[[type]])

    # HACK: remove intercept from auxiliary model matrix
    if (type == "aux") out <- out[,-1]
    return(out)
  } else stop("not enough information in fitted model to return model.matrix")
}

#' @rdname familyPrior
#' @export
familyPrior.wals <- function(object, ...) return(object$prior)
