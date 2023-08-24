#' WALS for zerotruncated count data models
#'
#' Implements the WALS method for count data models (only for GLMs so far)
#' WARNING: Interactions in formula do work work properly yet, see bugs.txt
#'
#' @export
walsZerotrunc <- function(x, ...) UseMethod("walsZerotrunc", x)

#' @export
walsZerotrunc.formula <- function(formula, data, subset, na.action, weights, offset,
                                  family, prior = weibull(), methodStart = "BFGS",
                                  controlStart = list(),
                                  model = TRUE, y = TRUE, x = FALSE, ...) {
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
  mt <- terms(formula, data = data)
  mtX1 <- terms(formula, data = data, rhs = 1L)
  mtX2 <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")

  if (any(Y <= 0)) stop("response Y <= 0 not allowed")

  n <- length(Y)
  X1 <- model.matrix(mtX1, mf)
  X2 <- model.matrix(mtX2, mf)
  
  # save contrasts before dropping constant in X2
  # lose attributes when we apply HACK
  cont <- list(focus = attr(X1, "contrasts"), aux = attr(X2, "contrasts"))
  
  ## HACK ##
  # remove intercept in X2
  # TODO: can we do this more elegantly via formula or terms or contrasts?
  # issue if add "-1" in formula --> estimate all levels of a factor because
  # without intercept, can estimate one additional level. But we do not want this.
  # want to keep the same levels as before and only remove constant column.
  X2 <- X2[,-1L] # intercept is always first column

  k1 <- ncol(X1)
  k2 <- ncol(X2)

  # check if X1 and X2 contain the same variables
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")

  ## weights (not used yet)
  weights <- processWeights(weights, mf, n)

  ## offsets (not used yet)
  offset <- getOffset(formula, mf, cl, n)

  # generate starting values
  # TODO: weights and offset not implemented yet
  betaStart <- zerotruncBetaStart(cbind(X1, X2), Y, family = family,
                                  weights = rep(1.0, times = length(Y)),
                                  offset = rep(0.0, times = length(Y)),
                                  method = methodStart,
                                  control = controlStart)$par
  betaStart1 <- betaStart[1L:k1]
  betaStart2 <- betaStart[(k1 + 1L):(k1 + k2)]


  ## call workhorse
  out <- walsZerotrunc.fit(X1 = X1, X2 = X2, y = Y, betaStart1 = betaStart1,
                     betaStart2 = betaStart2, family = family, prior = prior,
                     ...)

  # add more elements
  out$call <- cl
  out$formula <- oformula
  out$terms <- list(focus = mtX1, aux = mtX2, full = mt)
  out$levels <- list(focus = .getXlevels(mtX1, mf), aux = .getXlevels(mtX2, mf),
                     full = .getXlevels(mt, mf))
  out$contrasts <- cont
  # overwrite family from walsGLM.fit which returns only family$zerotrunc
  out$family <- family
  # out$fitted.link <- X1 %*% out$beta1 + X2 %*% out$beta2
  if (model) out$model <- mf
  if (y) out$y <- Y
  if (x) out$x <- list(focus = X1, aux = X2)


  class(out) <- c("walsZerotrunc", "wals")
  return(out)
}


#' @export walsZerotrunc.fit
walsZerotrunc.fit <- function(X1, X2, y, betaStart1, betaStart2,
                              family, prior = weibull(), ...) {
  # For GLMs it uses default walsGLM.fit, for NegBin I need to modify

  walsGLM.fit(X1, X2, y, betaStart1, betaStart2, family = family$zerotrunc,
              prior = prior, ...)
}


# Class methods ----------------------------------------------------------------

#' @export
coef.walsGLM <- function(object) return(object$coef)

#' @export
vcov.walsGLM <- function(object) return(object$vcovBeta)

#' @export
summary.walsZerotrunc <- function(object, ...) {
  out <- summary.walsGLM(object, ...)
  class(out) <- "summary.walsZerotrunc"
  return(out)
}

#' @export
print.summary.walsZerotrunc <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  dist <- x$family$untrunc$family
  link <- x$family$untrunc$link
  prior <- x$prior$prior

  cat(paste0("\nCoefficients (zerotruncated ", dist," with ", link, " link", " and ",
             prior, " prior): \n"))
  printCoefmat(x$coefficients, digits = digits,...)

  cat(paste0("\nStarting values for estimation: \n"))
  print(signif(x$betaStart, digits))

  invisible(x)
}



#' Predict method for walsZerotrunc Fits
#'
#' Obtains predictions from a fitted walsZerotrunc object
#'
#' @param type the type of prediction required. The default (\code{type = "response"})
#' is the mean of the zerotruncated ddistribution, i.e. it is on the scale of the response variable.
#' \code{type = "link"} gives the linear predictor, \code{type = "variance"}
#' the variance of the zerotruncated distribution. Further, \code{type = "zero"}
#' returns the probability for observing counts larger than zero (i.e. P(Y > 0))
#' and \code{type = "count"} returns the predicted mean of the untruncated distribution.
#' Finally, \code{type = "prob"} returns the probability for counts.
#'
#' @export
predict.walsZerotrunc <- function(object, newdata,
                            type = c("response", "link", "variance",
                                     "zero", "count", "prob"),
                            na.action = na.pass, at = NULL, ...) {
  # TODO: include offsets
  type <- match.arg(type)
  fam <- object$family$untrunc$family

  if (missing(newdata)) {
    switch(type,
           "response" = {
             return(object$fitted.values)
           },
           "link" = {
             return(object$fitted.link)
           },
           "variance" = {
             return(object$family$zerotrunc$variance(object$fitted.values))
           },
           "zero" = {
             # TODO: Currently only guaranteed to work for poisson...
             # Implement for others as well

             # if (fam == "poisson") {
             #   mu <- object$family$untrunc$linkinv(object$fitted.link)
             #   return(ppois(0, lambda = mu))
             # } else stop(paste(fam, "not supported yet"))

             return(object$family$untrunc$distFun(0, eta = object$fitted.link,
                                        lower.tail = FALSE, ...))

           },
           "count" = {
             return(object$family$untrunc$linkinv(object$fitted.link))
           },
           "prob" = {
             stop("predicted probabilities cannot be computed with missing newdata")
           })

  } else {

    # compute link
    newMatrices <- genNewdata(terms = object$terms,
                              contrasts = object$contrasts,
                              newdata = newdata, na.action = na.action,
                              xlev = object$levels)

    link <- newMatrices$X1 %*% object$beta1 + newMatrices$X2 %*% object$beta2

    switch(type,
           "response" = {
             return(object$family$zerotrunc$linkinv(link))
           },
           "link" = {
             return(link)
           },
           "variance" = {
             mu <- object$family$untrunc$linkinv(link)
             return(object$family$zerotrunc$variance(mu))
           },
           "zero" = {
             return(object$family$untrunc$distFun(0, eta = link,
                                        lower.tail = FALSE, ...))
           },
           "count" = {
             return(object$family$untrunc$linkinv(link))
           },
           "prob" = {
             if (!is.null(object$y)) y <- object$y
             else if (!is.null(object$model)) y <- model.response(object$model)
             else if (!is.null(at)) y <- at
             else stop(c("predicted probabilities cannot be
                         computed for fits with y = FALSE, model = FALSE
                         and at = NULL"))
             
             yUnique <- if (is.null(at)) 1:max(y) else at
             return(predictCounts(object$family, yUnique = yUnique,
                                  rowNames = rownames(newMatrices$X1),
                                  eta = link, ...))
           })

  }

}


# Helper functions -------------------------------------------------------------

#' Get starting values for WALS estimation
#'
#' Get starting values for WALS estimation with walsZerotrunc.fit via
#' maximizing likelihood with optim. Initial values are obtained via
#' maximizing the likelihood of a Poisson model via glm.fit.
zerotruncBetaStart <- function(X, Y, family, weights, offset, method, control) {
  start <- glm.fit(X, Y, family = poisson(), weights = weights,
                   offset = offset)$coefficients
  # optim minimizes, so we need to multiply loglik and gr by -1 to maximize
  optim(par = start, fn = family$zerotrunc$nllfun, gr = family$zerotrunc$ngrad,
        X = X, Y = Y, linkinv=family$untrunc$linkinv,
        weights = weights, offset = offset, method = method)
}

