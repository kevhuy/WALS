#' Weighted Average Least Squares for Generalized Linear Models
#'
#' Fits an NB2 regression model using the Weighted-Average Least Squares method
#' described in \insertCite{deluca2018glm;textual}{WALS}.
#'
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @export
walsGLM <- function(x, ...) UseMethod("walsGLM", x)

#' \code{walsGLM.formula} uses formulas to specify the design matrix.
#' @rdname walsGLM
#'
#' @param formula an object of class "\link{Formula}"
#' (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' The details of model specification are given under ‘Details’.
#' @param family Object of class \link{familyWALS}.
#' @param data an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment(formula), typically the environment from which walsNB is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action **not implemented yet.**
#' @param weights **not implemented yet.**
#' @param offset **not implemented yet.**
#' @param prior Object of class \code{familyPrior}. For example \link[WALS]{weibull}
#' or \link[WALS]{laplace}.
#' @param controlGLMfit Controls estimation of starting values for one-step ML,
#' passed to \link[stats]{glm.fit}. See also \link[stats]{glm.control}.
#' @param model if \code{TRUE} (default), then the model.frame is stored in
#' the return.
#' @param keepY if \code{TRUE} (default), then the response is stored in
#' the return.
#' @param keepX if \code{TRUE}, then the model matrix is stored in the return.
#' the return.
#' @param iterate if TRUE then the WALS algorithm is iterated using the previous
#' estimates as starting values
#' @param tol Only used if iterate = TRUE and nIt = NULL. If the Euclidean distance
#' between the previous beta and current beta falls below tol, then the algorithm stops.
#' @param maxIt Only used it iterate = TRUE and nIt = NULL. Aborts iterative fitting
#' when number of iterations exceed maxIt
#' @param nIt Only used if iterate = TRUE. If this is specified, then tol is ignored
#' and the algorithm iterates nIt times. This option should not be used unless
#' the user has a specific reason to run the algorithm nIt times, e.g. for
#' replication purposes.
#' @param verbose If verbose = TRUE, then it prints the iteration process
#' (only relevant if iterate = TRUE).
#' @param ... Arguments for workhorse \link[WALS]{walsGLMfit}.
#'
#' @details
#' Formulas should always contain two parts, i.e. they should be of the form
#' "y ~ X11 + X12 | X21 + X22", where the variables before "|" are the focus
#' regressors (includes a constant by default) and the ones after "|" are the
#' auxiliary regressors.
#'
#' **WARNING:** Interactions in formula do work work properly yet.
#' It is recommended to manually create the interactions beforehand and then
#' to insert them as 'linear terms' in the formula.
#'
#' @returns For \code{walsGLM.formula}, it returns an object of class
#' \code{walsGLM} which inherits from \code{walsGLM}. This is a list that contains
#' all elements returned from \link[WALS]{walsGLMfitIterate} and additionally
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
#' See returns of \link[WALS]{walsGLMfit} and \link[WALS]{walsGLMfitIterate}
#' for more details.
#'
#' @examples
#' data("HMDA", package = "AER")
#' fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
#'                        selfemp + afam, data = HMDA, family = binomialWALS(),
#'                        prior = weibull())
#' summary(fitBinomial)
#'
#' data("NMES1988", package = "AER")
#' fitPoisson <- walsGLM(emergency ~ health + chronic + age + gender |
#'                       I((age^2)/10) + married + region, data = NMES1988,
#'                       family = poissonWALS(), prior = laplace())
#' summary(fitPoisson)
#'
#' @export
walsGLM.formula <- function(formula, family, data, subset = NULL,
                            na.action = NULL, weights = NULL, offset = NULL,
                            prior = weibull(), controlGLMfit = list(),
                            model = TRUE, keepY = TRUE, keepX = FALSE,
                            iterate = FALSE, tol = 1e-6, maxIt = 50, nIt = NULL,
                            verbose = FALSE, ...) {

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

  # remove intercept in 2nd part by default
  # if 2nd part already includes "-1" it will remain unchanged
  # formula <- update(formula, . ~ . | . + -1)

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

  ## extract terms, model matrix, response
  mm <- extractModel(formula, mf, data)

  # extract objects from mm
  Y <- mm$Y; X1 <- mm$X1; X2 <- mm$X2; mt <- mm$mt; mtX1 <- mm$mtX1; mtX2 <- mm$mtX2
  cont <- mm$cont
  n <- length(Y)

  ## weights (not used yet)
  weights <- processWeights(weights, mf, n)

  ## offsets (not used yet)
  offset <- getOffset(formula, mf, cl, n)

  ## Fit model
  out <- walsGLMfitIterate(Y, X1, X2, family, na.action, weights, offset, prior,
                           controlGLMfit, keepY, keepX, iterate, tol, maxIt,
                           nIt, verbose, ...)

  # add more elements
  out$call <- cl
  out$formula <- oformula
  out$terms <- list(focus = mtX1, aux = mtX2, full = mt)
  out$levels <- list(focus = .getXlevels(mtX1, mf), aux = .getXlevels(mtX2, mf),
                     full = .getXlevels(mt, mf))
  out$contrasts <- cont
  if (model) out$model <- mf


  class(out) <- c("walsGLM", "wals")
  return(out)
}

#' \code{walsGLM.matrix()} uses prespecified design matrices x (focus) and
#' x2 (auxiliary) and response vector y.
#' @rdname walsGLM
#'
#' @param x Design matrix for focus regressors. Usually includes a constant
#' (column full of 1s) and can be generated using model.matrix().
#' @param x2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param y Count response as vector
#'
#' @returns For \code{walsGLM.matrix}, it returns an object of class
#' \code{walsGLMmatrix}, which inherits from \code{walsGLM}, \code{walsMatrix}
#' and \code{wals}. This is a list that contains all elements returned from
#' \link[WALS]{walsGLMfitIterate} and additionally the call in \code{cl}.
#'
#' @examples
#' ## Example for walsGLM.matrix()
#' data("HMDA", package = "AER")
#' X <- model.matrix(deny ~ pirat + hirat + lvrat + chist + mhist + phist + selfemp + afam,
#'                   data = HMDA)
#' X1 <- X[,c("(Intercept)", "pirat", "hirat", "lvrat", "chist2", "chist3",
#'         "chist4", "chist5", "chist6", "mhist2", "mhist3", "mhist4", "phistyes")]
#' X2 <- X[,c("selfempyes", "afamyes")]
#' y <- HMDA$deny
#' fit <- walsGLM(X1, X2, y, family = binomialWALS(), prior = weibull())
#' summary(fit)
#'
#' @export
walsGLM.matrix <- function(x, x2, y, family, subset = NULL, na.action = NULL,
                           weights = NULL, offset = NULL,
                           prior = weibull(), controlGLMfit = list(),
                           keepY = TRUE, keepX = FALSE,
                           iterate = FALSE, tol = 1e-6, maxIt = 50, nIt = NULL,
                           verbose = FALSE, ...) {
  cl <- match.call()
  X1 <- x
  X2 <- x2
  if (!is.null(subset)) {
    X1 <- X1[subset,]; X2 <- X2[subset,]; y <- y[subset]
  }

  out <- walsGLMfitIterate(y, X1, X2, family, na.action, weights, offset, prior,
                           controlGLMfit, keepY, keepX, iterate, tol, maxIt,
                           nIt, verbose, ...)

  out$call <- cl

  class(out) <- c("walsGLMmatrix", "walsGLM", "walsMatrix", "wals")
  return(out)
}

#' @rdname walsGLM
#' @export
walsGLM.default <- function(x, ...) {
  # inspired by glmboost.default in mboost.
  if (extends(class(x), "matrix")) {
    return(walsGLM.matrix(x, ...))
  }
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}

#' Fitter function for Weighted Average Least Squares estimation of GLMs
#'
#' Workhorse function behind \link[WALS]{walsGLM} and used internally in
#' \link[WALS]{walsGLMfitIterate}.
#'
#' @usage walsGLMfit(
#'  X1,
#'  X2,
#'  y,
#'  betaStart1,
#'  betaStart2,
#'  family,
#'  prior = weibull(),
#'  ...)
#'
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1s) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param y response as vector.
#' @param betaStart1 Starting values for coefficients of focus regressors X1.
#' @param betaStart2 Starting values for coefficients of auxiliary regressors X2.
#' @param family Object of class \link{familyWALS}.
#' @param prior Object of class \code{familyPrior}. For example \link[WALS]{weibull}
#' or \link[WALS]{laplace}.
#' @param ... Further arguments passed to \link[WALS]{walsFit}.
#'
#' @details
#' Uses \link[WALS]{walsFit} under the hood after transforming the regressors
#' \code{X1} and \code{X2} and the response \code{y}. For more details, see
#' \insertCite{huynhwals}{WALS} and \insertCite{deluca2018glm;textual}{WALS}.
#'
#' @returns A list containing all elements returned by \link[WALS]{walsFit},
#' except for \code{residuals}, and additionally (some fields are replaced)
#' \item{condition}{Condition number of the matrix
#' \eqn{\bar{\Xi} = \bar{\Delta}_{2} \bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2} \bar{\Delta}_{2}}.}
#' \item{family}{Object of class \link{familyWALS}. The family used.}
#' \item{betaStart}{Starting values of the regression coefficients for the
#' one-step ML estimators.}
#' \item{fitted.link}{Linear link fitted to the data.}
#' \item{fitted.values}{Estimated conditional mean for the data. Lives on the
#' scale of the response.}
#'
#'
#' @seealso [walsGLM], [walsGLMfitIterate], [walsFit].
#'
#' @references
#' \insertAllCited{}
#'
#' @export
walsGLMfit <- function(X1, X2, y, betaStart1, betaStart2,
                       family, prior = weibull(), ...) {
  X1names <- colnames(X1)
  X2names <- colnames(X2)
  Xnames <- c(X1names, X2names)

  etaStart <- X1 %*% betaStart1 + X2 %*% betaStart2
  X1start <- family$transformX(X1, etaStart, y)
  X2start <- family$transformX(X2, etaStart, y)
  yStart <- family$transformY(y, X1start, X2start, betaStart1, betaStart2, etaStart)


  # use generic WALS algo for linear models
  fit <- walsFit(X1 = X1start, X2 = X2start, y = yStart, sigma = 1,
                  prior = prior, prescale = TRUE, ...)

  fit$family <- family
  fit$betaStart <- c(betaStart1, betaStart2)
  fit$fitted.link <- drop(X1 %*% fit$beta1 + X2 %*% fit$beta2)
  fit$fitted.values <- family$linkinv(fit$fitted.link)
  fit$residuals <- NULL

  return(fit)
}

#' Iteratively fitting walsGLM, internal function for walsGLM.formula and
#' walsGLM.matrix.
#'
#' See description of \link[WALS]{walsGLM}.
#'
#' @param y Response as vector.
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1s) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param family Object of class \link{familyWALS}.
#' @param na.action Not implemented yet.
#' @param weights Not implemented yet.
#' @param offset Not implemented yet.
#' @param prior Object of class \code{familyPrior}, e.g. \link[WALS]{weibull}.
#' @param controlGLMfit Controls estimation of starting values for one-step ML,
#' passed to \link[stats]{glm.fit}. See also \link[stats]{glm.control}.
#' @param keepY If \code{TRUE}, then output keeps response.
#' @param keepX If \code{TRUE}, then output keeps the design matrices.
#' @param iterate if TRUE then the WALS algorithm is iterated using the previous
#' estimates as starting values.
#' @param tol Only used if iterate = TRUE and nIt = NULL. If the Euclidean distance
#' between the previous beta and current beta falls below tol and the absolute
#' difference between the previous and current rho falls below tol, then
#' the algorithm stops.
#' @param maxIt Only used if iterate = TRUE and nIt = NULL. Aborts iterative fitting
#' when number of iterations exceed maxIt.
#' @param nIt Only used if iterate = TRUE. If this is specified, then tol is ignored
#' and the algorithm iterates nIt times.
#' @param verbose If verbose = TRUE, then it prints the iteration process
#' (only relevant if iterate = TRUE).
#' @param ... Arguments to be passed to the workhorse function walsGLMfit.
#'
#' @returns A list containing all elements returned from \link[WALS]{walsGLMfit}
#' and additionally the following elements:
#' \item{y}{If \code{keepY = TRUE}, contains the response vector.}
#' \item{x}{list. If \code{keepX} is true, then it is a list with elements
#' \code{x1} and \code{x2} containing the design matrices of the focus and
#' auxiliary regressors, respectively.}
#' \item{weights}{returns the argument \code{weights}.}
#' \item{offset}{returns the argument \code{offset}.}
#' \item{converged}{Logical. Only relevant if \code{iterate = TRUE}. Equals
#' \code{TRUE} if iterative fitting converged, else \code{FALSE}. Is \code{NULL}
#' if \code{iterate = FALSE}.}
#' \item{it}{Number of iterations run in the iterative fitting algorithm.
#' \code{NULL} if \code{iterate = FALSE}.}
#' \item{deviance}{Deviance of the fitted regression model.}
#' \item{residuals}{Raw residuals, i.e. response - fitted mean.}
#'
#' @seealso [walsGLM], [walsGLMfit].
#'
#' @export
walsGLMfitIterate <- function(y, X1, X2, family, na.action = NULL,
                              weights = NULL, offset = NULL,
                              prior = weibull(), controlGLMfit = list(),
                              keepY = TRUE, keepX = FALSE, iterate = FALSE,
                              tol = 1e-6, maxIt = 50, nIt = NULL,
                              verbose = FALSE, ...) {
  # Useful quantities
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  # check if X1 and X2 contain the same variables
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")

  # generate starting values
  betaStart <- glm.fit(cbind(X1, X2), y, family = family,
                       control = controlGLMfit)$coefficients
  betaStart1 <- betaStart[1L:k1]
  betaStart2 <- betaStart[(k1 + 1L):(k1 + k2)]

  # simply reuse iterative code by setting nIt = 1
  if (!iterate) nIt <- 1

  betaOld <- rep(NA, length(betaStart)) # make sure loop starts
  betaCurrent <- betaStart

  it <- 0
  converged <- FALSE
  if (!is.null(nIt)) {
    maxIt <- nIt
  }

  it <- 0
  for (i in 1:maxIt) {

    ## update values
    betaOld <- betaCurrent

    ## call workhorse
    out <- walsGLMfit(X1 = X1, X2 = X2, y = y, betaStart1 = betaCurrent[1L:k1],
                       betaStart2 = betaCurrent[(k1 + 1L):(k1 + k2)],
                       family = family, prior = prior,
                       ...)

    betaCurrent <- out$coef
    it <- it + 1

    if (verbose) cat(paste("\r finished iteration", it))

    if (is.null(nIt) && (norm(betaOld - betaCurrent, type = "2") < tol)) {
      converged <- TRUE
      cat("\n algorithm converged")
      break
    }

  }

  if (!is.null(nIt)) {
    converged <- NULL
  } else if (!converged) cat("\n algorithm failed to converge")

  # replace starting values with original starting values
  out$betaStart <- betaStart

  # add more elements
  if (keepY) out$y <- family$initializeY(y) # e.g. convert logical to 0s and 1s.
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$weights <- weights
  out$offset <- offset
  out$converged <- converged
  out$it <- if (iterate) it else NULL

  # deviance & residuals
  wt <- if (is.null(weights)) rep(1, nrow(X1)) else weights
  mu <- out$fitted.values
  out$deviance <- sum(family$dev.resids(out$y, mu, wt))
  out$residuals <- out$y - mu

  return(out)
}


# Class methods ----------------------------------------------------------------

#' @export
print.walsGLM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print.wals(x, digits, ...)
  cat(paste0("\nResidual Deviance: ", signif(x$deviance, digits), "\n"))
  invisible(x)
}

#' @export
summary.walsGLM <- function(object, ...) {
  object <- summary.wals(object, ...)

  # inspired by summary.glm() in stats
  if (!is.na(object$family$dispersion)) {
    object$dispersion <- object$family$dispersion
  } else object$dispersion <- NA


  class(object) <- "summary.walsGLM"
  return(object)
}

#' @export
print.summary.walsGLM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  printCallCoefs(x, digits, ...)

  cat(paste0("\nStarting values for estimation: \n"))
  print.default(format(x$betaStart, digits = digits), print.gap = 2, quote = FALSE)

  # inspired by print.summary.glm() from stats
  if (!is.na(x$dispersion)) {
    cat(paste0("\n(Dispersion parameter for ", x$family$family,
        " family taken to be ", format(x$dispersion),")\n"))
  } else {
    cat("\n(", x$family$family, " family used)\n")
  }

  cat(paste0("\nResidual deviance: ", signif(x$deviance, max(5L, digits + 1L)),
             "\n"))

  printPriorNKappa(x, digits)

  if (!is.null(x$it)) {
    if (!is.null(x$converged)) {
      convStatus <- if (x$converged) "converged " else "did not converge "
      cat(paste0("\nFitting algorithm ", convStatus, "in ", x$it, " iterations.\n"))
    } else if (is.null(x$converged)) {
      # manually ran nIt iterations of fitting algo without checking for conv.
      cat(paste0("\nFitting algorithm run for ", x$it, " iterations.\n"))
    }
  }

  invisible(x)
}


#' @export
predict.walsGLM <- function(object, newdata,
                            type = c("response", "link", "variance", "prob",
                                     "density", "logDens"),
                            at = NULL,
                            na.action = na.pass, log = FALSE, ...) {
  # TODO: include offsets
  type <- match.arg(type)

  # convenience type for scores.R
  if (type == "logDens") {
    type <- "density"
    log <- TRUE
  }

  if (missing(newdata)) {
    link <- object$fitted.link

    if (type == "density") {
      y <- residuals(object, type = "response") + fitted(object)
    } else y <- NULL

  } else {
    # compute link
    newMatrices <- genNewdata(object$terms, object$contrasts, newdata,
                              na.action = na.action, xlev = object$levels)
    link <- drop(newMatrices$X1 %*% object$beta1 + newMatrices$X2 %*% object$beta2)

    if (type == "density") {
      y <- getY(terms(object, "focus"), newdata, na.action = na.action)
    } else y <- NULL

  }
  return(.predictGLM(object, link, y, type, at, log))
}

#' @export
predict.walsGLMmatrix <- function(object, newX1, newX2, newY = NULL,
                                  type = c("response", "link", "variance", "prob",
                                           "density", "logDens"),
                                  at = NULL, log = FALSE, ...) {
  # TODO: include offsets
  type <- match.arg(type)

  # convenience type
  if (type == "logDens") {
    type <- "density"
    log <- TRUE
  }

  if (missing(newX1) || missing(newX2)) {
    return(predict.walsGLM(object, type = type, at = at, log = log, ...))
  } else {
    link <- newX1 %*% object$beta1 + newX2 %*% object$beta2
    return(.predictGLM(object, link, newY, type, at, log))
  }
}

.predictGLM <- function(object, link, y, type, at, log, ...) {
  switch(type,
         "response" = {
           return(object$family$linkinv(link))
         },
         "link" = {
           return(link)
         },
         "variance" = {
           mu <- object$family$linkinv(link)
           return(object$family$variance(mu))
         },
         "prob" = {
           if (!is.null(object$y)) y <- object$y
           else if (!is.null(object$model)) y <- model.response(object$model)
           else if (!is.null(at)) y <- at
           else stop(c("predicted probabilities cannot be
                         computed for fits with y = FALSE, model = FALSE
                         and at = NULL"))

           if (any(at < 0)) stop("prediction at count < 0 not allowed")

           yUnique <- if (is.null(at)) 0:max(y) else at
           return(predictCounts(object$family, yUnique = yUnique,
                                rowNames = names(link),
                                eta = link, ...))
         },
         "density" = {
           return(object$family$density(y, link, log = log))
         })
}

#' @export
residuals.walsGLM <- function(object, type = c("deviance", "pearson", "response"),
                              ...) {
  type <- match.arg(type)
  y <- object$residuals + fitted(object)
  mu <- fitted(object)
  wt <- if (is.null(object$weights)) rep(1, length(y)) else object$weights

  switch(type,
         deviance = { # inspired by stats:::residuals.glm()
           dres <- sqrt(pmax((object$family$dev.resids)(y, mu, wt), 0))
           return(ifelse(y > mu, dres, -dres))
         },
         pearson = {
           return((y - mu) * sqrt(wt)/sqrt(object$family$variance(mu)))
         },
         response = {return(object$residuals)}
  )
}

#' @export
logLik.walsGLM <- function(object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  y <- residuals(object, type = "response") + fitted(object)
  return(sum(object$family$density(y, object$fitted.link, log = TRUE)))
}

#' @export
familyWALS.walsGLM <- function(object, ...) return(object$family)
