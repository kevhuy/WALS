
#' @export
walsHurdle <- function(x, ...) UseMethod("walsHurdle", x)

#' Fits Hurdle count models with the Weighted-Average Least Squares method
#'
#' WARNING: Interactions in formula do work work properly yet, see bugs.txt
#' @export
walsHurdle.formula <- function(formula, data, subset, na.action, weights, offset,
                               family, prior = weibull(),
                               priorZero = weibull(),
                               methodStart = "BFGS",
                               controlStart = list(count = list(), zero = list()),
                               model = TRUE, y = TRUE, x = FALSE, method = "original",
                               ...) {
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

  if (length(formula)[2L] < 4L) {
    # TODO: Implement what happens when only one part is specified in formula
    stop("One part formula not implemented yet")
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  # mt <- terms(formula, data = data)
  Y <- model.response(mf, "numeric")
  Ybinary <- Y > 0
  Ytrunc <- Y[Ybinary]
  mtX <- terms(formula, data = data, rhs = 1L:2L)
  mtX1 <- terms(formula, data = data, rhs = 1L)
  mtX2 <- delete.response(terms(formula, data = data, rhs = 2L))
  mtZ <- terms(formula, data = data, rhs = 3L:4L)
  mtZ1 <- delete.response(terms(formula, data = data, rhs = 3L))
  mtZ2 <- delete.response(terms(formula, data = data, rhs = 4L))

  if (any(Y < 0)) stop("response Y < 0 not allowed")

  n <- length(Y)
  X1 <- model.matrix(mtX1, mf)
  X2 <- model.matrix(mtX2, mf)
  Z1 <- model.matrix(mtZ1, mf)
  Z2 <- model.matrix(mtZ2, mf)


  # save contrasts before dropping constants in X2 and Z2
  # lose attributes when we apply the HACK below
  cont <- list(count = NULL, zero = NULL)
  cont$count <- list(focus = attr(X1, "contrasts"), aux = attr(X2, "contrasts"))
  cont$zero <- list(focus = attr(Z1, "contrasts"), aux = attr(Z2, "contrasts"))

  ## HACK ##
  # remove intercept in X2 and Z2
  # TODO: can we do this more elegantly via formula or terms or contrasts?
  # issue if add "-1" in formula --> estimate all levels of a factor because
  # without intercept, can estimate one additional level. But we do not want this.
  # want to keep the same levels as before and only remove constant column.
  X2 <- X2[,-1L] # intercept is always first column
  Z2 <- Z2[,-1L]

  k1x <- ncol(X1)
  k2x <- ncol(X2)
  k1z <- ncol(Z1)
  k2z <- ncol(Z2)

  # check if X1 and X2 contain the same variables, same for Z1 and Z2
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")
  if (any(colnames(Z1) %in% colnames(Z2))) stop("Z1 and Z2 contain the same variables")


  ## weights (not used yet)
  weights <- processWeights(weights, mf, n)

  ## offsets (not used yet)
  offset <- getOffset(formula, mf, cl, n)

  # generate starting values for count part
  # TODO: weights and offset not implemented yet
  betaStart <- zerotruncBetaStart(cbind(X1[Ybinary,], X2[Ybinary,]), Ytrunc,
                                  family = family$count,
                                  weights = rep(1.0, times = length(Ytrunc)),
                                  offset = rep(0.0, times = length(Ytrunc)),
                                  method = methodStart,
                                  control = controlStart$count)$par
  betaStart1 <- betaStart[1L:k1x]
  betaStart2 <- betaStart[(k1x + 1L):(k1x + k2x)]

  # generate starting values for zero part
  gammaStart <- glm.fit(cbind(Z1, Z2), factor(Ybinary), family = family$zero,
                        control = controlStart$zero)$coefficients
  gammaStart1 <- gammaStart[1L:k1z]
  gammaStart2 <- gammaStart[(k1z + 1L):(k1z + k2z)]



  ## call workhorse
  out <- walsHurdle.fit(X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, y = Ytrunc,
                        yBinary = Ybinary,
                        betaStart1 = betaStart1, betaStart2 = betaStart2,
                        gammaStart1 = gammaStart1, gammaStart2 = gammaStart2,
                        family = family$count, familyZero = family$zero, prior = prior,
                        priorZero = priorZero, method = method, ...)

  # add more elements
  out$call <- cl
  out$formula <- oformula
  out$count$terms <- list(focus = mtX1, aux = mtX2, full = mtX)
  out$zero$terms <- list(focus = mtZ1, aux = mtZ2, full = mtZ)
  out$count$levels <- list(focus = .getXlevels(mtX1, mf), aux = .getXlevels(mtX2, mf),
                     full = .getXlevels(mtX, mf))
  out$zero$levels <- list(focus = .getXlevels(mtZ1, mf), aux = .getXlevels(mtZ2, mf),
                          full = .getXlevels(mtZ, mf))
  out$count$contrasts <- cont$count
  out$zero$contrasts <- cont$zero

  # overwrite count$family because walsHurdle only returns zerotrunc family
  # --> also need untruncated family for easy printing in summary
  out$count$family <- family$count
  # out$zero$family <- familyZero
  out$family <- family
  if (model) out$model <- mf
  if (y) out$y <- Y
  if (x) out$x <- list(count = list(focus = X1, aux = X2),
                       zero = list(focus = Z1, aux = Z2))

  ## Add fitted values
  ## TODO: Currently only works for poisson, how do we solve this elegantly
  ## for negbin where we need to input 2nd parameter apart from eta?

  # overwrite fitted.link and fitted.values
  # from count to use all observations and not only those with Y > 0
  out$count$fitted.link <- X1 %*% out$count$beta1 + X2 %*% out$count$beta2
  out$count$fitted.values <- family$count$zerotrunc$linkinv(out$count$fitted.link)

  p0Zero <- out$zero$fitted.values
  out$fitted.values <- p0Zero * out$count$fitted.values


  class(out) <- c("walsHurdle", "wals")
  return(out)

}

#' Fitting function for WALS Hurdle models
#'
#' Uses walsGLM.fit for the binary response (y > 0 vs y = 0)
#' and walsZerotrunc.fit for fitting a count data model on all y > 0
#' @export walsHurdle.fit
walsHurdle.fit <- function(X1, X2, Z1, Z2, y, yBinary, betaStart1, betaStart2,
                           gammaStart1, gammaStart2, family, familyZero,
                           prior, priorZero, method, ...) {
  out <- list(count = NULL, zero = NULL)
  # drop = FALSE so we keep it as a matrix if X1 is column vector (n x 1)
  out$count <- walsZerotrunc.fit(X1 = X1[yBinary,,drop = FALSE],
                                 X2 = X2[yBinary,,drop = FALSE],
                                 y = y, betaStart1 = betaStart1,
                                 betaStart2 = betaStart2, family = family,
                                 prior = prior, method, ...)
  out$zero <- walsGLM.fit(X1 = Z1, X2 = Z2, y = yBinary, betaStart1 = gammaStart1,
                         betaStart2 = gammaStart2, family = familyZero,
                         prior = priorZero, method, ...)

  return(out)
}


# Class methods ----------------------------------------------------------------

#' @export
coef.walsHurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  switch(model,
         "full" = {
          return(structure(c(object$count$coef, object$zero$coef),
                          .Names = c(paste("count", names(object$count$coef),
                                          sep = "_"),
                                    paste("zero", names(object$zero$coef),
                                          sep = "_")))
                         )},
        "count" = {
          return(object$count$coef)},
        "zero" = {
          return(object$zero$coef)}
  )
}

#' @export
vcov.walsHurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)

  kc <- length(object$count$coef)
  kz <- length(object$zero$coef)
  varnames <- names(coef(object, model))
  vcov <- switch(model,
         "full" = {
           rbind(cbind(object$count$vcovBeta, matrix(0, kc, kz)),
                 cbind(matrix(0, kz, kc), object$zero$vcovBeta))
         },
         "count" = object$count$vcovBeta,
         "zero" = object$zero$vcovBeta
         )

  colnames(vcov) <- rownames(vcov) <- varnames
  return(vcov)
}

#' @export
summary.walsHurdle <- function(object, ...) {
  out <- list(count = NULL, zero = NULL)
  out$count <- summary.walsZerotrunc(object$count)
  out$zero <- summary.walsGLM(object$zero)
  class(out) <- "summary.walsHurdle"
  return(out)
}

#' @export
print.summary.walsHurdle <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCount model: ------------------ \n")
  print.summary.walsZerotrunc(x$count, digits = digits, ...)
  cat("\nZero model: ------------------- \n")
  print.summary.walsGLM(x$zero, digits = digits, ...)
}

#' @export
terms.walsHurdle <- function(object, model = c("count", "zero"), ...) {
  model <- match.arg(model)
  switch(model,
         "count" = return(object$count$terms),
         "zero" = return(object$zero$terms),
  )

}

#' @export
predict.walsHurdle <- function(object, newdata,
                               type = c("response", "prob", "count", "zero",
                                        "density"),
                               na.action = na.pass, at = NULL,
                               log.p = FALSE, # log.p is for type = "prob"
                               log = FALSE, # log is for type = "density"
                               ...) {
  # TODO: include offsets
  type <- match.arg(type)
  famCount <- object$family$count
  famZero <- object$family$zero

  if (missing(newdata)) {
    switch(type,
           "response" = {
             return(object$fitted.values)
           },
           "zero" = {
             # TODO: Currently only guaranteed to work for poisson...
             # Implement for others as well
             p0Count <- famCount$untrunc$distFun(0,eta = object$count$fitted.link,
                                                 lower.tail = FALSE,
                                                 log.p = TRUE, log = FALSE, ...)
             p0Zero <- log(object$zero$fitted.values)

             if (log.p) return(p0Zero - p0Count) else return(exp(p0Zero - p0Count))

           },
           "count" = {
             return(famCount$untrunc$linkinv(object$count$fitted.link))
           },
           "prob" = {
             stop("predicted probabilities cannot be computed with missing newdata")
           },
           "density" = {
             if (is.null(object$y)) {
               stop("predicted density cannot be computed without y")
             } else {
               return(object$family$density(y, object$count$fitted.link,
                                  object$zero$fitted.link, log = log, ...))
             }
           })

  } else {

    # compute link
    newMatricesCount <- genNewdata(terms = object$count$terms,
                                   contrasts = object$count$contrasts,
                                   newdata = newdata, na.action = na.action,
                                   xlev = object$count$levels)
    newMatricesZero <- genNewdata(terms = object$zero$terms,
                                  contrasts = object$zero$contrasts,
                                  newdata = newdata, na.action = na.action,
                                  xlev = object$zero$levels)


    linkCount <- (newMatricesCount$X1 %*% object$count$beta1 +
                  newMatricesCount$X2 %*% object$count$beta2)
    linkZero <- (newMatricesZero$X1 %*% object$zero$beta1 +
                   newMatricesZero$X2 %*% object$zero$beta2)

    # TODO: only works if famZero = binomialWALS() --> mean = probability
    # for bernoulli! Reimplement if allow other distributions.
    p0Zero <- famZero$linkinv(linkZero)

    switch(type,
           "response" = {
             return(p0Zero * famCount$zerotrunc$linkinv(linkCount))
           },
           "zero" = {
             # TODO: Currently only guaranteed to work for poisson...
             # Implement for others as well
             p0Count <- famCount$untrunc$distFun(0,eta = linkCount,
                                                 lower.tail = FALSE,
                                                 log.p = TRUE, ...)
             if (log.p) return(p0Zero - p0Count) else return(exp(p0Zero - p0Count))
           },
           "count" = {
             return(famCount$untrunc$linkinv(linkCount))
           },
           "prob" = {
             # TODO: Currently only guaranteed to work for poisson
             # because of definition of p0Count
             # Implement for others as well
             if (!is.null(object$y)) y <- object$y
             else if (!is.null(object$model)) y <- model.response(object$model)
             else if (!is.null(at)) y <- at
             else stop(c("predicted probabilities cannot be
                         computed for fits with y = FALSE, model = FALSE
                         and at = NULL"))

             if (any(at < 0)) stop("prediction at count < 0 not allowed")
             p0Count <- famCount$untrunc$distFun(0,eta = linkCount,
                                                 lower.tail = FALSE,
                                                 log.p = TRUE, ...)
             yUnique <- if (is.null(at)) 0:max(y) else at
             return(predictCountsHurdle(famCount, famZero,
                                        yUnique = yUnique,
                                        rowNames = rownames(newMatricesCount$X1),
                                        etaCount = linkCount,
                                        p0Zero = p0Zero,
                                        etaZero = linkZero,...))
           },
           "density" = {
             y <- getY(terms(object, model = "count")$focus, newdata, na.action = na.action)
             return(object$family$density(y, linkCount, linkZero, log = log))
           })

  }
}
