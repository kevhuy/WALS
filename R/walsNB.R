#' Weighted Average Least Squares Negative Binomial
#'
#' Fits an NB2 regression model using the Weighted-Average Least Squares method
#' of \insertCite{huynhwalsnb;textual}{WALS}.
#'
#'
#' @return For \code{walsNB.formula}, it returns an object of class
#' \code{walsNB} which inherits from \code{walsGLM} and \code{wals}.
#' See return of \link[WALS]{walsNB.fit} for more details.
#'
#' For \code{walsNB.matrix}, it returns an object of class \code{walsNBmatrix},
#' which inherits from \code{walsNB}, \code{walsGLM} and \code{wals}.
#' See return of \link[WALS]{walsNB.fit} for more details.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
walsNB <- function(x, ...) UseMethod("walsNB", x)


#' Fits a NB2 regression with the Weighted-Average Least Squares method
#'
#' **WARNING:** Interactions in formula do work work properly yet.
#' It is recommended to manually create the interactions beforehand and then
#' to insert them as 'linear terms' in the formula.
#'
#' \code{walsNB.formula} uses formulas to specify the design matrix.
#' @rdname walsNB
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
#' @param link specifies the link function, currently only "log" is supported.
#' @param prior Object of class \code{familyPrior}. For example \link{weibull}
#' or \link{laplace}. Not tested with other priors.
#' @param model if \code{TRUE} (default), then the model.frame is stored in
#' the return.
#' @param keepY if \code{TRUE} (default), then the response is stored in
#' the return.
#' @param keepX if \code{TRUE}, then the model matrix is stored in the return.
#' the return.
#' @param controlInitNB Controls estimation of starting values for one-step ML,
#' see \link[WALS]{controlNB}.
#' @param iterate if TRUE then the WALS algorithm is iterated using the previous
#' estimates as starting values
#' @param tol Only used if iterate = TRUE and nIt = NULL. If the Euclidean distance
#' between the previous beta and current beta falls below tol and the absolute difference between
#' the previous and current rho falls below tol, then the algorithm stops.
#' @param maxIt Only used it iterate = TRUE and nIt = NULL. Aborts iterative fitting
#' when number of iterations exceed maxIt
#' @param nIt Only used if iterate = TRUE. If this is specified, then tol is ignored
#' and the algorithm iterates nIt times. This option should not be used unless
#' the user has a specific reason to run the algorithm nIt times, e.g. for
#' replication purposes.
#' @param verbose If verbose = TRUE, then it prints the iteration process of
#' internal function \link[WALS]{walsNB.fitIterate} (only relevant if iterate = TRUE).
#' @param ... Arguments for workhorse \link[WALS]{walsNB.fit}.
#'
#'
#' @details
#' Formulas should always contain two parts, i.e. they should be of the form
#' "y ~ X11 + X12 | X21 + X22", where the variables before "|" are the focus
#' regressors (includes a constant by default) and the ones after "|" are the
#' auxiliary regressors.
#'
#' @examples
#' data("NMES1988", package = "AER")
#'
#' fitWeibull <- walsNB(visits ~ health + chronic + age + gender | I((age^2)/10) +
#'                      married + region, data = NMES1988, prior = weibull())
#' summary(fitWeibull)
#'
#' fitLaplace <- walsNB(visits ~ health + chronic + age + gender | I((age^2)/10) +
#'                      married + region, data = NMES1988, prior = laplace())
#' summary(fitLaplace)
#'
#' @export
walsNB.formula <- function(formula, data, subset = NULL, na.action = NULL,
                           weights = NULL, offset = NULL,
                           link = "log", prior = weibull(),
                           controlInitNB = controlNB(),
                           model = TRUE, keepY = TRUE, keepX = FALSE,
                           iterate = FALSE, tol = 1e-6, maxIt = 10000, nIt = NULL,
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

  # extract terms, model matrix, response
  mm <- extractModel(formula = formula, mf = mf, data = data)

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
  out <- walsNB.fitIterate(Y, X1, X2, link, na.action, weights, offset,
                           prior, controlInitNB, keepY, keepX,
                           iterate, tol, maxIt, nIt, verbose, ...)


  # add more elements
  out$call <- cl
  out$formula <- oformula
  out$terms <- list(focus = mtX1, aux = mtX2, full = mt)
  out$levels <- list(focus = .getXlevels(mtX1, mf), aux = .getXlevels(mtX2, mf),
                     full = .getXlevels(mt, mf))
  out$contrasts <- cont
  if (model) out$model <- mf


  class(out) <- c("walsNB", "walsGLM", "wals")
  return(out)
}



#' Fits a NB2 regression with the Weighted-Average Least Squares method
#'
#' \code{walsNB.matrix()} uses prespecified design matrices X1 and X2 and
#' response vector y.
#' @rdname walsNB
#'
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1's) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param y Count response as vector
#' @export
walsNB.matrix <- function(X1, X2, y, link = "log", subset = NULL,
                          na.action = NULL, weights = NULL, offset = NULL,
                          prior = weibull(), controlInitNB = controlNB(),
                          model = TRUE, keepY = TRUE, keepX = FALSE,
                          iterate = FALSE, tol = 1e-6, maxIt = 10000, nIt = NULL,
                          verbose = FALSE, ...) {
  cl <- match.call()

  if (!is.null(subset)) {
    X1[subset,] <- X1; X2[subset,]; y <- y[subset]
  }

  out <- walsNB.fitIterate(y, X1, X2, link, na.action, weights, offset,
                           prior, controlInitNB, keepY, keepX,
                           iterate, tol, maxIt, nIt,
                           verbose, ...)

  out$call <- cl
  class(out) <- c("walsNBmatrix", "walsNB", "walsGLMmatrix", "walsGLM", "wals")
  return(out)
}

#' @rdname walsNB
#' @export
walsNB.default <- function(x, ...) {
  # inspired by glmboost.default in mboost.
  if (extends(class(x), "matrix")) {
    return(walsNB.matrix(X1 = x, ...))
  }
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}



#' Fitter function for Weighted Average Least Squares estimation of NB2 regression model
#'
#'
#' Workhorse function behind \link{walsNB} and used internally in
#' \link{walsNB.fitIterate}.
#'
#' @usage walsNB.fit(
#'  X1,
#'  X2,
#'  y,
#'  betaStart1,
#'  betaStart2,
#'  rhoStart,
#'  family,
#'  prior,
#'  method = c("fullSVD", "original"),
#'  svdTol = .Machine$double.eps,
#'  svdRtol = 1e-6,
#'  keepUn = FALSE,
#'  keepR = FALSE,
#'  eigenSVD = TRUE,
#'  postmult = TRUE,
#'  ...
#'  )
#'
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1's) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param y Count response as vector
#' @param betaStart1 Starting values for coefficients of focus regressors X1.
#' @param betaStart2 Starting values for coefficients of auxiliary regressors X2.
#' @param rhoStart Starting value for log-dispersion parameter of NB2
#' @param family Object of class \code{family}. Currently only supports
#' \code{negativeBinomial()}.
#' @param prior Object of class \code{familyPrior}. For example \link[WALS]{weibull}
#' or \link[WALS]{laplace}. Not tested with other priors.
#' @param method Specifies method used. Available methods are \code{"fullSVD"}
#' (default) or \code{"original"}. See details.
#' @param svdTol Tolerance for rank of matrix \eqn{\bar{Z}_{1}} and \eqn{\bar{Z}}.
#' Only used if \code{method = "fullSVD"}.
#' Checks if smallest eigenvalue in SVD of \eqn{\bar{Z}_1} and \eqn{\bar{Z}}
#' is larger than \code{svdTol}, otherwise reports a rank deficiency.
#' @param svdRtol Relative tolerance for rank of matrix \eqn{\bar{Z}_{1}} and \eqn{\bar{Z}}.
#' Only used if \code{method = "fullSVD"}. Checks if ratio of largest to smallest
#' eigenvalue in SVD of \eqn{\bar{Z}_1} and \eqn{\bar{Z}} is larger than
#' \code{svdRtol}, otherwise reports a rank deficiency.
#' @param keepUn If \code{TRUE}, keeps the one-step ML estimators of the
#' unrestricted model, i.e. \eqn{\tilde{\gamma}_{u}} and \eqn{\tilde{\beta}_{u}}.
#' @param keepR If \code{TRUE}, keeps the one-step ML estimators of the fully
#' restricted model, i.e. \eqn{\tilde{\gamma}_{r}} and \eqn{\tilde{\beta}_{r}}.
#' @param eigenSVD If \code{TRUE}, then \code{semiorthogonalize()} uses \code{svd()}
#' to compute the eigendecomposition of \eqn{\bar{Xi}} instead of \code{eigen()}.
#' In this case, the tolerances of \code{svdTol} and \code{svdRtol} are used to
#' determine whether \eqn{\bar{Xi}} is of full rank (need it for \eqn{\bar{Xi}^{-1/2}}).
#' @param postmult If \code{TRUE} (default), then it computes
#' \deqn{\bar{\Xi}^{-1/2} = H \Lambda^{-1/2} H^{\top}} instead of
#' \deqn{\bar{\Xi}^{-1/2} = H \Lambda^{-1/2}.}
#' The latter is used in the original MATLAB code for WALS in the linear regression model
#' \insertCite{magnus2010growth,deluca2011stata,kumar2013normallocation,magnus2016wals}{WALS},
#' see eq. (12) of \insertCite{magnus2016wals;textual}{WALS} for more details.
#' @param ... Arguments for internal function \code{computePosterior()}.
#'
#'
#'
#' @return A list containing the coefficients in \code{coefs}.
#'
#' \code{beta1} and \code{beta2} are the model-averaged estimates of coefficients for the focus
#' and auxiliary regressors, respectively.
#'
#' \code{gamma1} and \code{gamma2} are the model-averaged estimates of the
#' coefficients of the transformed focus and auxiliary regressors.
#'
#' \code{rho} is the model-averaged estimate of the
#' log-dispersion parameter of the NB2 distribution.
#'
#' \code{prior} contains the \code{prior} specified in the arguments.
#'
#' If \code{keepUn = TRUE}, then \code{betaUn1} and \code{betaUn2} contain the
#' unrestricted one-step ML estimators of the coefficients of the focus and
#' auxiliary regressors. \code{gammaUn1} and \code{gammaUn2} are the fully
#' unrestricted estimates for the transformed focus and auxiliary regressors.
#'
#' If \code{keepR}, then \code{gamma1r} contains the fully restricted one-step
#' ML estimator for the transformed regressors (only focus regressors).
#'
#' \code{family} contains the \link[stats]{family} object from
#' \link[WALS]{negbinWALS} used later for predictions.
#'
#' \code{betaStart} and \code{rhoStart} contain the starting values for the
#' one-step ML estimators.
#'
#' \code{fitted.link} contains the linear link fitted to the data.
#'
#' \code{fitted.values} contains the estimated conditional mean for the data.
#' Lives on the scale of the response.
#'
#'
#'
#' @details The method to be specified in \code{method} mainly differ in the way
#' they compute the fully restricted and unrestricted estimators for the
#' transformed regressors \eqn{Z}, i.e. \eqn{\tilde{\gamma_{1r}}},
#' \eqn{\gamma_{1u}} and \eqn{\tilde{\gamma_{2u}}}.
#'
#' \itemize{
#' \item{"fullSVD"}{Recommended approach. First applies an SVD to \eqn{\bar{Z}_{1}}
#' to compute \eqn{\bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2}}:
#' It is used for computing the inverse of
#'
#' \deqn{\bar{X}_{1}^{\top}\bar{X}_{1}
#'  + \bar{g} \bar{\epsilon} X_{1}^{\top}\bar{q} \bar{q}^{\top} X_{1}},
#'
#' when using the Sherman-Morrison-Woodbury formula. We further leverage the
#' SVD of \eqn{\bar{Z}_1} and additionally \eqn{\bar{Z}} to compute the fully
#' unrestricted estimator \eqn{\tilde{\gamma}_{u}} and the fully restricted
#' estimator \eqn{\tilde{\gamma}_{r}}. For \eqn{\tilde{\gamma}_{u}}, we simply
#' use the SVD of \eqn{\bar{Z}} to solve the full equation system derived from
#' the one-step ML problem
#' (see section "Simplifications for \eqn{\tilde{\gamma}_{u}}" in
#' the appendix of \insertCite{huynhwals;textual}{WALS} for more details.
#' The SVD of \eqn{\bar{Z}_1} is further used in computing the model averaged
#' estimator for the focus regressors \eqn{\hat{\gamma}_1} (see section
#' "Simplifications for \eqn{\hat{\gamma}_1}" in
#' \insertCite{huynhwals;textual}{WALS} for more details).
#'
#' Described in more detail in the appendix
#' of \insertCite{huynhwals;textual}{WALS}.}
#'
#'
#' \item{"original"} {Computes all inverses directly using \code{solve()} and
#' does not make use of the Sherman-Morrison-Woodbury formula for certain
#' inverses. Specifically: directly inverts the matrix
#' \eqn{\bar{Z}_{1}^{\top} \bar{Z}_{1}} using \code{solve()}
#' in order to compute \eqn{\bar{M}_1}. Moreover, it computes the fully
#' unrestricted estimators of the focus regressors
#' \eqn{\tilde{\gamma}_{1u}} and of the auxiliary regressors
#' \eqn{\tilde{\gamma}_{2u}} and the fully restricted estimator
#' \eqn{\tilde{\gamma}_{1r}} by directly implementing the formulas derived
#' in \insertCite{huynhwalsnb;textual}{WALS}.
#' This method should only be used as reference and for easier
#' debugging.}
#' }
#'
#' All variables in the code that are contain "start" in their name feature the
#' starting values for the one-step ML estimation of submodels. See section
#' "One-step ML estimator" of \insertCite{huynhwalsnb}{WALS} for details.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @export walsNB.fit
walsNB.fit <- function(X1, X2, y, betaStart1, betaStart2, rhoStart, family,
                       prior, method = c("fullSVD", "original"),
                       svdTol = .Machine$double.eps,
                       svdRtol = 1e-6, keepUn = FALSE, keepR = FALSE,
                       eigenSVD = TRUE, postmult = TRUE, ...) {
  ## Sanity checks
  method <- match.arg(method)

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


  X1names <- colnames(X1)
  X2names <- colnames(X2)
  Xnames <- c(X1names, X2names)

  # transform X and y
  etaStart <- X1 %*% betaStart1 + X2 %*% betaStart2
  X1start <- family$transformX(X1, etaStart, y)
  X2start <- family$transformX(X2, etaStart, y)
  y0Start <- family$transformY0(y, X1start, X2start, betaStart1, betaStart2,
                                etaStart)


  ## Step 2.a: Scaling X1 so all diag elements of (X1*Delta1)'X1*Delta1 = 1 ----
  ## Corresponds to equation (8) of De Luca et al. 2018 except division by n
  ## missing
  Delta1 <- 1.0 / sqrt(colSums(X1start^2.0))

  # multiply each row with Delta1
  Z1 <- multAllRows(X1, Delta1)
  Z1start <- multAllRows(X1start, Delta1)

  # generate useful vectors
  qStart <- as.vector(family$q(etaStart, y))
  gStart <- family$g()
  epsilonStart <- family$epsilon(etaStart, y)
  psiStart <- family$psi(etaStart, y)

  ## Step 2.b: Scaling X2 so diag els. of (X2*Delta2)'M1*X2*Delta2 = 1 ---------
  if (method == "fullSVD") {
    svdZ1start <- svd(Z1start) # does economy SVD
    singularVals <- svdZ1start$d

    checkSingularitySVD(singularVals, tol = svdTol, rtol = svdRtol)

    ellStart <- crossprod(svdZ1start$u, qStart/sqrt(psiStart))
    UellStart <- svdZ1start$u %*% ellStart


    ## Step 2.c Constructing Xi ------------------------------
    B <- gStart * epsilonStart * sum(ellStart^2)
    geB <- (gStart*epsilonStart / (1 + B))

    X2M1X2 <- computeX2M1X2(X2, X2start, qStart, svdZ1start$u, UellStart,
                            ellStart, psiStart,
                            gStart, epsilonStart, geB)

  } else if (method == "original") {
    # empty entries for svd
    svdZ1start <- NULL

    ## Sherman-Morrison-Woodbury for (X1bar'X1bar + g*epsilon  * X1'q*q'X1)^-1
    # Z1startinv <- solve(crossprod(Z1start, Z1start))
    #
    # Delta1 * (Z1bar'Z1bar)^-1 Delta1 = (X1bar'X1bar)^-1
    # first multiply all columns with Delta1, then all rows!
    # X1startinv <- multAllRows(Delta1 * Z1startinv, Delta1)
    # X1X1plusinv <- ( X1startinv +
    #                    ( (X1startinv %*% tcrossprod(X1tq, X1tq) %*% X1startinv) /
    #                        as.numeric((1 + gStart*epsilonStart*crossprod(X1tq, X1startinv %*% X1tq)))
    #                    )
    # )

    X1tq <- colSums(X1 * qStart)
    X1X1plusinv <-  solve(crossprod(X1start, X1start) + gStart*epsilonStart*tcrossprod(X1tq, X1tq))

    X2tq <- colSums(X2 * qStart)
    A <- crossprod(X2start, X1start) +  gStart*epsilonStart*tcrossprod(X2tq, X1tq)
    X2M1X2 <- (crossprod(X2start, X2start) + gStart*epsilonStart*tcrossprod(X2tq, X2tq) -
                 tcrossprod(A %*% X1X1plusinv, A))

  } else stop(paste("method ==", method, "is not implemented."))

  Delta2 <- 1.0 / sqrt(diag(X2M1X2))
  Xibar <- multAllRows(Delta2 * X2M1X2, Delta2) # Xibar = Delta2 * X2bar'M1 X2bar * Delta2
  # note we ignore normalization by n! --> can ignore all 1/n in the derived formulas!
  # They will all cancel anyways!


  ## Step 3: Transformation of X2 to Z2 --------------------------------------
  outSemiOrt <- semiorthogonalize(Xibar, X2, Delta2, eigenSVD, postmult)


  ## Step 4: Unrestricted  estimators ----------------------------------------
  Z2 <- outSemiOrt$Z2
  Z2start <- family$transformX(Z2, etaStart, y) # \bar{Z2}

  # need transformation to matrix (k1 x 1 and k2 x 1 matrices) to avoid
  # errors later with crossprod and tcrossprod in computeGammaUn
  Z1tq <- as.matrix(colSums(Z1 * qStart))
  Z2tq <- as.matrix(colSums(Z2 * qStart))

  tStart <- family$t(etaStart, y)

  if (method == "fullSVD") {
    # use SVD to solve eq. system to get gammaUn, applies SVD on complete Z
    # again.
    svdZstart <- svd(cbind(Z1start, Z2start))
    ellStartZ <- crossprod(svdZstart$u, qStart/sqrt(psiStart))

    checkSingularitySVD(svdZstart$d, tol = svdTol, rtol = svdRtol)

    gammaUn <- computeGammaUnSVD(svdZstart$u, svdZstart$v, svdZstart$d,
                                 ellStartZ, gStart, epsilonStart, qStart,
                                 y0Start, tStart, psiStart)
    gammaUn1 <- gammaUn[1:k1]
    gammaUn2 <- gammaUn[(k1 + 1):k]

  } else if (method == "original") {
    Z1y0s <- crossprod(Z1start, y0Start) - (tStart*epsilonStart)*Z1tq
    Z2y0s <- crossprod(Z2start, y0Start) - (tStart*epsilonStart)*Z2tq

    Z1Z1plusinv <- solve(crossprod(Z1start, Z1start) + gStart*epsilonStart*tcrossprod(Z1tq, Z1tq))

    ## Sherman-Morrison-Woodbury for (Z1bar'Z1bar + g*epsilon  * Z1'q*q'Z1)^-1
    # Z1Z1plusinv <- (   Z1startinv -
    #                      ( (Z1startinv %*% tcrossprod(Z1tq, Z1tq) %*% Z1startinv) /
    #                          as.numeric((1 + gStart*epsilonStart*crossprod(Z1tq, Z1startinv %*% Z1tq)))
    #                      )
    # )

    A <- crossprod(Z2start, Z1start) + gStart*epsilonStart*tcrossprod(Z2tq, Z1tq)
    gammaUn2 <- -(A %*% Z1Z1plusinv %*% Z1y0s) + Z2y0s

    Z2Z1s <- crossprod(Z2start, Z1start) + gStart*epsilonStart*tcrossprod(Z2tq, Z1tq)
    gammaUn1 <- (Z1Z1plusinv %*% (  Z1y0s +
                                      crossprod(Z2Z1s, Z2Z1s) %*% Z1Z1plusinv %*% Z1y0s -
                                      crossprod(Z2Z1s, Z2y0s)  ))

    # Preparation for Step 6 & 7
    D <- tcrossprod(Z1Z1plusinv, Z2Z1s)


  }


  ## Step 5: Compute posterior mean of x  ~ N(\gamma, 1)
  gamma2 <- computePosterior(prior, gammaUn2, ...)$postMean


  ## Step 6 & 7: WALS estimates ------------------------------------------------
  if (method == "fullSVD") {
    gamma1 <- computeGamma1(gamma2 = gamma2, Z2start = Z2start,
                            Z2 = Z2, U = svdZ1start$u, V = svdZ1start$v,
                            singularVals = singularVals, ellStart = ellStart,
                            gStart = gStart, epsilonStart = epsilonStart,
                            qStart = qStart, y0Start = y0Start,
                            tStart = tStart, psiStart = psiStart)

    if (keepR) {
      ## TODO: Replace computeGamma1r with computeGammaUnSVD as the computation
      ## of gammaUn and gamma1r is exactly the same except that we need
      ## the SVD of Z1 as input instead of Z.
      gamma1r <- computeGamma1r(svdZ1start$u, svdZ1start$v, singularVals,
                                ellStart, gStart, epsilonStart, qStart, y0Start,
                                tStart, psiStart)
    }

  } else if (method == "original") {
    gamma1r <- gammaUn1 + D %*% gammaUn2
    gamma1 <- gamma1r - D %*% gamma2
  }

  alpha <- family$computeAlpha(gamma1, gamma2, Z1, Z2, y, etaStart, qStart)
  rho <- exp(alpha)

  # convert to betas
  beta1 <- as.vector(Delta1 * gamma1)
  beta2 <- as.vector(outSemiOrt$D2 %*% gamma2)

  if (keepUn) {
    betaUn1 <- as.vector(Delta1 * gammaUn1)
    betaUn2 <- as.vector(outSemiOrt$D2 %*% gammaUn2)
    names(betaUn1) <- X1names
    names(betaUn2) <- X2names
  }


  fit <- list(coef = c(beta1, beta2), beta1 = beta1, beta2 = beta2,
              rho = as.numeric(rho), # otherwise it's 1x1 matrix
              gamma1 = gamma1, gamma2 = gamma2,
              vcovBeta = NULL, vcovGamma = NULL,
              betaStart = c(betaStart1, betaStart2), prior = prior,
              betaUn1 = if (keepUn) betaUn1 else NULL,
              betaUn2 = if (keepUn) betaUn2 else NULL,
              gammaUn1 = if (keepUn) gammaUn1 else NULL,
              gammaUn2 = if (keepUn) gammaUn2 else NULL,
              gamma1r = if (keepR) gamma1r else NULL)

  # contains previous rho, used for obtaining betaStart etc.
  fit$familyStart <- family
  # new family with model averaged rho
  fit$family <- negbinWALS(rho, link = family$link)
  fit$X1names <- X1names
  fit$X2names <- X2names
  fit$betaStart <- c(betaStart1, betaStart2)
  fit$rhoStart <- rhoStart
  fit$fitted.link <- drop(X1 %*% fit$beta1 + X2 %*% fit$beta2)
  fit$fitted.values <- family$linkinv(fit$fitted.link)
  fit$k1 <- k1
  fit$k2 <- k2
  fit$n <- n
  fit$condition <- outSemiOrt$condition

  # assign names to variables
  names(fit$coef) <- Xnames
  names(fit$betaStart) <- Xnames
  names(fit$beta1) <- names(fit$gamma1) <- X1names
  names(fit$beta2) <- names(fit$gamma2) <- X2names
  return(fit)
}

#' Iteratively fitting walsNB, internal function for walsNB.formula and
#' walsNB.vector.
#'
#' See description of \link{walsNB}.
#'
#' @param y Count response as vector
#' @param X1 Design matrix for focus regressors. Usually includes a constant
#' (column full of 1's) and can be generated using model.matrix().
#' @param X2 Design matrix for auxiliary regressors. Usually does not include
#' a constant column and can also be generated using model.matrix().
#' @param link specifies the link function, currently only "log" is supported.
#' @param na.action Not implemented yet.
#' @param weights Not implemented yet.
#' @param offset Not implemented yet.
#' @param prior Object of class \code{familyPrior}, e.g. \link[WALS]{weibull}.
#' @param controlInitNB Controls estimation of starting values for one-step ML,
#' see \link[WALS]{controlNB}.
#' @param keepY If \code{TRUE}, then output keeps response.
#' @param keepX If \code{TRUE}, then output keeps design matrix.
#' @param iterate if TRUE then the WALS algorithm is iterated using the previous
#' estimates as starting values
#' @param tol Only used if iterate = TRUE and nIt = NULL. If the Euclidean distance
#' between the previous beta and current beta falls below tol and the absolute difference between
#' the previous and current rho falls below tol, then the algorithm stops.
#' @param maxIt Only used if iterate = TRUE and nIt = NULL. Aborts iterative fitting
#' when number of iterations exceed maxIt.
#' @param nIt Only used if iterate = TRUE. If this is specified, then tol is ignored
#' and the algorithm iterates nIt times.
#' @param verbose If verbose = TRUE, then it prints the iteration process of
#' internal function walsNB.fitIterate (only relevant if we iterate = TRUE).
#' @param ... Arguments to be passed to the workhorse function walsNB.fit which
#' actually fits the model.
#'
#' @export walsNB.fitIterate
walsNB.fitIterate <- function(y, X1, X2, link = "log", na.action = NULL,
                              weights = NULL, offset = NULL,
                              prior = weibull(), controlInitNB = controlNB(),
                              keepY = TRUE, keepX = FALSE,
                              iterate = FALSE, tol = 1e-6, maxIt = 10000, nIt = NULL,
                              verbose = FALSE, ...) {

  ## Useful quantities
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  ## generate family
  if (is.null(controlInitNB$start$logTheta)) {
    startLogTheta <- 0 ## THIS THETA WILL NOT BE USED FOR INITIAL NB2 FIT!
  } else {
    startLogTheta <- controlInitNB$start$logTheta
  }
  familyStart <- negativeBinomial(theta = exp(startLogTheta),
                                  link = link)

  ## generate starting values via Maximum likelihood

  Xinit <- if (controlInitNB$restricted) X1 else cbind(X1, X2)
  nb2 <- fitNB2(Xinit, y, family = familyStart,
                control = controlInitNB)

  if (nb2$convergence == 1) {
    warning("Iteration limit reached in initial NB2 fit of full model.")
  } else if (nb2$convergence == 10) {
    warning("Degeneracy in Nelder-Mead simplex in initial NB2 fit of full model.")
  } else if (nb2$convergence == 51) {
    warning("Warning in L-BFGS-B in initial NB2 fit of full model.
            See initial fit for more details.")
  } else if (nb2$convergence == 52) {
    message("Error in L-BFGS-B in initial NB2 fit of full model.
        See initial fit for more details. \n")
    out <- list(initialFit = nb2)
    return(out)
  } else if (nb2$convergence == 99) {
    message("Convergence issue in IWLS algo in glm.nb for initial NB2 fit of full model.
            See initial fit for more details.")
  } else if (nb2$convergence != 0) {
    warning("Unknown non-zero exit for optim in initial NB2 fit of full model.
            See initial fit for more details.")
  }

  betaStart <- if (controlInitNB$restricted) {
    c(nb2$coefficients, rep(0, k2))
  } else nb2$coefficients

  betaStart1 <- betaStart[1L:k1]
  betaStart2 <- betaStart[(k1 + 1L):(k1 + k2)]
  rhoStart <- nb2$theta

  # simply reuse iterative code by setting nIt = 1
  if (!iterate) nIt <- 1


  ## generate new family with estimate starting value for rho and beta
  familyStart <- negbinWALS(scale = rhoStart, link = link)
  betaOld <- rep(NA, length(betaStart))
  betaCurrent <- betaStart
  rhoOld <- NA
  rhoCurrent <- rhoStart

  it <- 0
  converged <- FALSE
  if (!is.null(nIt)) {
    maxIt <- nIt
  }

  for (i in 1:maxIt) {
    ## update values
    betaOld <- betaCurrent
    rhoOld <- rhoCurrent


    ## generate new family with estimate starting value for rho and beta
    family <- negbinWALS(scale = rhoCurrent, link = link)

    ## call workhorse
    out <- walsNB.fit(X1 = X1, X2 = X2, y = y, betaStart1 = betaCurrent[1L:k1],
                      betaStart2 = betaCurrent[(k1 + 1L):(k1 + k2)],
                      rhoStart = rhoCurrent,
                      family = family, prior = prior,
                      ...)

    betaCurrent <- out$coef
    rhoCurrent <- out$rho
    it <- it + 1

    if (verbose) cat(paste("\r finished iteration", i))

    if (is.null(nIt) && (norm(betaOld - betaCurrent, type = "2") < tol && abs(rhoOld - rhoCurrent) < tol)) {
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
  out$rhoStart <- rhoStart
  out$familyStart <- familyStart

  # add more elements
  if (keepY) out$y <- y
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$initialFit <- nb2
  out$weights <- weights
  out$converged <- converged
  out$it <- if (iterate) it else NULL

  # deviance & residuals
  wt <- if (is.null(weights)) rep(1, nrow(X1)) else weights
  mu <- out$fitted.values
  out$deviance <- sum(family$dev.resids(out$y, mu, wt))
  out$residuals <- out$y - mu

  return(out)
}



## Class methods ---------------------------------------------------------------

#' @export
vcov.walsNB <- function(object, ...) {
  stop("No method for objects of class ", sQuote(class(object)[1]), " implemented.")
}

#' @export
summary.walsNB <- function(object, ...) {
  object <- summary.wals(object, ...)

  # remove SE estimation from walsGLM
  object$focusCoefs <- object$focusCoefs[,"Estimate", drop = FALSE]
  object$auxCoefs <- object$auxCoefs[,"Estimate", drop = FALSE]

  class(object) <- c("summary.walsNB")
  return(object)
}

#' @export
print.summary.walsNB <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")

  printCallCoefs(x, digits, ...)

  cat(paste0("\nStarting values for estimation: \n"))
  print.default(format(c(x$betaStart, "rho" = x$rhoStart), digits = digits),
                print.gap = 2, quote = FALSE)

  # inspired by print.summary.glm() from stats
  cat(paste0("\n(Dispersion parameter rho for Negative Binomial family estimated as ",
             signif(x$rho, digits),")\n"))

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


## Helper functions ------------------------------------------------------------

#' Define controllable parameters of initial NB fit
#'
#' @param start Optional starting values for \link[WALS]{fitNB2}. Only used if
#' \code{initMASS = FALSE}.
#' @param method Optimization method used in \link[stats]{optim}. Only used if
#' \code{initMASS = FALSE}.
#' @param controlOptim List with parameters controlling optimization process of
#' \link[stats]{optim}. Only used if \code{initMASS = FALSE}.
#' @param initThetaMASS If TRUE, then initial \eqn{\log{\theta}} of
#' \link[WALS]{fitNB2} is estimated using \link[MASS]{theta.ml}
#' (ML-estimation over 1 variable) based on regression coefficients from
#' Poisson regression. If \code{FALSE}, then initial \eqn{\log{theta}} = 0 is used.
#' @param initMASS If TRUE (default), then initial fit in \link[WALS]{fitNB2} is estimated via
#' \link[MASS]{glm.nb} and \code{initThetaMASS} is ignored.
#' @param restricted If TRUE, then initial fit in \link[WALS]{fitNB2} only considers the
#' focus regressors. By default \code{FALSE}, then the unrestricted model is
#' estimated in \link[WALS]{fitNB2} (i.e. all regressors).
#' @param eps Controls argument \code{eps} in \link[WALS]{fitNB2} for generating
#' starting value for \code{logTheta} (\eqn{\log{\theta}}) via \link[MASS]{theta.ml}.
#' @param epsilonMASS Sets epsilon in control argument of \link[MASS]{glm.nb}.
#'
#' @export
controlNB <- function(start = list(mu = NULL, logTheta = NULL), method = "BFGS",
                      controlOptim = list(maxit = 100), initThetaMASS = FALSE,
                      initMASS = TRUE, restricted = FALSE,
                      eps = .Machine$double.eps^0.25,
                      epsilonMASS = 1e-8) {
  return(list(start = start,  method = method, controlOptim = controlOptim,
              restricted = restricted, initThetaMASS = initThetaMASS,
              initMASS = initMASS, eps = eps, epsilonMASS = epsilonMASS))
}

