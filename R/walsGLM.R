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


#' Fitter function for Weighted Average Least Squares estimation of GLMs
#'
#' @usage walsGLM.fit(
#'  X1,
#'  X2,
#'  y,
#'  betaStart1,
#'  betaStart2,
#'  family,
#'  prior = weibull(),
#'  ...)
#'
#' @export walsGLM.fit
walsGLM.fit <- function(X1, X2, y, betaStart1, betaStart2,
                        family, prior = weibull(), ...) {

  X1names <- colnames(X1)
  X2names <- colnames(X2)
  Xnames <- c(X1names, X2names)

  etaStart <- X1 %*% betaStart1 + X2 %*% betaStart2

  X1start <- family$transformX(X1, etaStart, y)
  X2start <- family$transformX(X2, etaStart, y)
  yStart <- family$transformY(y, X1start, X2start, betaStart1, betaStart2, etaStart)


  # use generic WALS algo for linear models
  fit <- wals.fit(X1 = X1start, X2 = X2start, y = yStart, sigma = 1,
                  prior = prior, prescale = TRUE, ...)

  fit$family <- family
  fit$betaStart <- c(betaStart1, betaStart2)
  fit$fitted.link <- drop(X1 %*% fit$beta1 + X2 %*% fit$beta2)
  fit$fitted.values <- family$linkinv(fit$fitted.link)

  # class(fit) <- c("walsGLM", class(fit))
  return(fit)
}

#' @export
walsGLMfitIterate <- function(y, X1, X2, family, na.action = NULL,
                              weights = NULL, offset = NULL,
                              prior = weibull(), controlGlmFit = list(),
                              keepY = TRUE, keepX = FALSE, iterate = FALSE,
                              tol = 1e-6, maxIt = 10000, nIt = NULL,
                              verbose = FALSE, ...) {
  # Useful quantities
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  # check if X1 and X2 contain the same variables
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")

  # generate starting values
  betaStart <- glm.fit(cbind(X1, X2), y, family = family,
                       control = controlGlmFit)$coefficients
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
    out <- walsGLM.fit(X1 = X1, X2 = X2, y = y, betaStart1 = betaCurrent[1L:k1],
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
  out$converged <- converged

  # add more elements
  if (keepY) out$y <- y
  if (keepX) out$x <- list(focus = X1, aux = X2)

  return(out)
}


#' @rdname walsGLM
#' @export
walsGLM.default <- function(x, ...) {
  # inspired by glmboost.default in mboost.
  if (extends(class(x), "matrix")) {
    return(walsGLM.matrix(X1 = x, ...))
  }
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}

#' @rdname walsGLM
#' @export
walsGLM.matrix <- function(X1, X2, y, family, subset = NULL, na.action = NULL,
                           weights = NULL, offset = NULL,
                           prior = weibull(), controlGlmFit = list(),
                           keepY = TRUE, keepX = FALSE,
                           iterate = FALSE, tol = 1e-6, maxIt = 10000, nIt = NULL,
                           verbose = FALSE, ...) {
  cl <- match.call()

  if (!is.null(subset)) {
    X1[subset,] <- X1; X2[subset,]; y <- y[subset]
  }

  out <- walsGLMfitIterate(y, X1, X2, family, na.action, weights, offset, prior,
                           controlGlmFit, keepY, keepX, iterate, tol, maxIt,
                           nIt, verbose, ...)

  out$call <- cl

  class(out) <- c("walsMatrix", "wals")
  return(out)
}

#' Fits a GLM with the Weighted-Average Least Squares method
#'
#' **WARNING:** Interactions in formula do work work properly yet.
#' It is recommended to manually create the interactions beforehand and then
#' to insert them as 'linear terms' in the formula.
#'
#' @rdname walsGLM
#'
#' @param formula an object of class "\link{Formula}"
#' (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' The details of model specification are given under ‘Details’.
#' @param data an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment(formula), typically the environment from which walsNB is called.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights **not implemented yet.**
#' @param offset **not implemented yet.**
#' @param na.action **not implemented yet.**
#' @param iterate if TRUE then the WALS algorithm is iterated using the previous
#' estimates as starting values
#' @param tol Only used if iterate = TRUE and nIt = NULL. If the Euclidean distance
#' between the previous beta and current beta falls below tol, then the algorithm stops.
#' @param maxIt Only used it iterate = TRUE and nIt = NULL. Aborts iterative fitting
#' when number of iterations exceed maxIt
#' @param nIt Only used if iterate = TRUE. If this is specified, then tol is ignored
#' and the algorithm iterates nIt times.
#' @param verbose If verbose = TRUE, then it prints the iteration process
#' (only relevant if iterate = TRUE).
#' @param ... Arguments for workhorse \link[WALS]{walsGLM.fit}.
#'
#' @details
#' Formulas should always contain two parts, i.e. they should be of the form
#' "y ~ X11 + X12 | X21 + X22", where the variables before "|" are the focus
#' regressors (includes a constant by default) and the ones after "|" are the
#' auxiliary regressors.
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
                            prior = weibull(), controlGlmFit = list(),
                            model = TRUE, keepY = TRUE, keepX = FALSE,
                            iterate = FALSE, tol = 1e-6, maxIt = 10000, nIt = NULL,
                            verbose = FALSE, ...) {

  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
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
                           controlGlmFit, keepY, keepX, iterate, tol, maxIt,
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


# Class methods ----------------------------------------------------------------

#' @export
summary.walsGLM <- function(object, ...) {

  se <- sqrt(diag(object$vcovBeta))
  object$se <- se
  object$coefficients <- cbind(object$coef, se)
  colnames(object$coefficients) <- c("Estimate", "Std. Error")

  class(object) <- "summary.walsGLM"
  return(object)
}

#' @export
print.summary.walsGLM <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  dist <- x$family$family
  link <- x$family$link
  prior <- x$prior$prior

  cat(paste0("\nCoefficients (", dist," with ", link, " link", " and ",
             prior, " prior): \n"))
  printCoefmat(x$coefficients, digits = digits,...)

  cat(paste0("\nStarting values for estimation: \n"))
  print(signif(x$betaStart, digits))

  invisible(x)
}

#' @export
coef.walsGLM <- function(object) return(object$coef)

#' @export
vcov.walsGLM <- function(object) return(object$vcovBeta)

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
    switch(type,
           "response" = {
             return(object$fitted.values)
           },
           "link" = {
             return(object$family$linkfun(object$fitted.values))
           },
           "variance" = {
             return(object$family$variance(object$fitted.values))
           },
           "prob" = {
             stop("predicted probabilities cannot be computed with missing newdata")
           }
           # TODO: add density prediction
           )

  } else {

    # compute link
    newMatrices <- genNewdata(object$terms, object$contrasts, newdata,
                              na.action = na.action, xlev = object$levels)

    link <- drop(newMatrices$X1 %*% object$beta1 + newMatrices$X2 %*% object$beta2)

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
                                        rowNames = rownames(newMatrices$X1),
                                        eta = link, ...))
           },
           "density" = {
             y <- getY(terms(object, "focus"), newdata, na.action = na.action)
             return(object$family$density(y, link, log = log))
           })

  }

}




