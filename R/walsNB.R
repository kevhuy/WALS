#' Weighted-Average Least Squares for Negative Binomial Regression
#'
#' Performs model averaging for NB2 regression models using the Weighted-Average
#' Least Squares method of \insertCite{huynhwalsnb;textual}{WALS}.
#'
#' @details
#' Computes WALS estimates when focus regressors (X1) are present in all
#' submodels and model averaging takes place over the auxiliary regressors (X2).
#'
#' @references
#' \insertAllCited{}
#'
#' @export
walsNB <- function(x, ...) UseMethod("walsNB", x)

#' \code{walsNB.formula()} uses formulas to specify the design matrix.
#' @rdname walsNB
#'
#' @inheritParams walsGLM.formula
#' @param link specifies the link function, currently only \code{"log"} is supported.
#' @param controlInitNB Controls estimation of starting values for one-step ML,
#' see \code{\link[WALS]{controlNB}}.
#' @param tol Only used if \code{iterate = TRUE} and \code{nIt = NULL}. If the
#' Euclidean distance between the previous and current coefficient vector divided
#' by the square root of the length of the vector falls below \code{tol} and the
#' absolute difference between the previous and current dispersion parameter
#' falls below \code{tol}, then the algorithm stops.
#' See \code{\link[WALS]{walsNBfitIterate}} for more details.
#' @param verbose If \code{verbose = TRUE}, then it prints the iteration process
#' of internal function \code{\link[WALS]{walsNBfitIterate}}
#' (only relevant if \code{iterate = TRUE}).
#' @param ... Arguments for workhorse \code{\link[WALS]{walsNBfit}}.
#'
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
#' See \code{\link[WALS]{predict.walsGLM}} and \code{\link[WALS]{predict.wals}}
#' for some class methods that the fitted objects inherit from
#' \code{"\link[WALS]{walsGLM}"} and \code{"\link[WALS]{wals}"}, respectively.
#'
#' @returns \code{walsNB.formula()} returns an object of class \code{"walsNB"}
#' which inherits from \code{"walsGLM"} and \code{"wals"}. This is a list that
#' contains all elements returned from \code{\link[WALS]{walsNBfitIterate}} and
#' additionally
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
#' See returns of \code{\link[WALS]{walsNBfit}} and \code{\link[WALS]{walsNBfitIterate}}
#' for more details.
#'
#' @examples
#' ## Example for walsNB.formula()
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
  out <- walsNBfitIterate(Y, X1, X2, link, na.action, weights, offset,
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



#' \code{walsNB.matrix()} uses prespecified design matrices x (focus) and
#' x2 (auxiliary) and response vector y.
#' @rdname walsNB
#' @aliases walsNBmatrix
#'
#' @inheritParams wals.matrix
#' @param y Count response as vector.
#'
#' @returns \code{walsNB.matrix()} returns an object of class \code{"walsNBmatrix"},
#' which inherits from \code{"walsNB"}, \code{"walsGLMmatrix"}, \code{"walsGLM"}
#' and \code{"wals"}. This is a list that contains all elements returned from
#' \code{\link[WALS]{walsNBfitIterate}} and additionally the call in \code{cl}.
#'
#' @examples
#' ## Example for walsNB.matrix()
#' data("NMES1988", package = "AER")
#' X <- model.matrix(visits ~ health + chronic + age + gender + married + region,
#'                   data = NMES1988)
#' X1 <- X[, c("(Intercept)", "healthpoor", "healthexcellent", "chronic",
#'         "age", "gendermale")]
#' X2 <- X[, c("marriedyes", "regionnortheast", "regionmidwest", "regionwest")]
#' y <- NMES1988$visits
#' fit <- walsNB(X1, X2, y, prior = weibull())
#' summary(fit)
#'
#' @export
walsNB.matrix <- function(x, x2, y, link = "log", subset = NULL,
                          na.action = NULL, weights = NULL, offset = NULL,
                          prior = weibull(), controlInitNB = controlNB(),
                          model = TRUE, keepY = TRUE, keepX = FALSE,
                          iterate = FALSE, tol = 1e-6, maxIt = 50, nIt = NULL,
                          verbose = FALSE, ...) {
  cl <- match.call()
  X1 <- x
  X2 <- x2
  if (!is.null(subset)) {
    X1 <- X1[subset,]; X2 <- X2[subset,]; y <- y[subset]
  }

  out <- walsNBfitIterate(y, X1, X2, link, na.action, weights, offset,
                           prior, controlInitNB, keepY, keepX,
                           iterate, tol, maxIt, nIt,
                           verbose, ...)

  out$call <- cl
  class(out) <- c("walsNBmatrix", "walsNB", "walsGLMmatrix", "walsGLM",
                  "walsMatrix", "wals")
  return(out)
}

#' @rdname walsNB
#'
#' @details
#' \code{walsNB.default()} raises an error if \code{x} is not an object of class
#' \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise
#' it calls \code{walsNB.matrix()}. It is a modified version of \code{glmboost.default}
#' from the \code{mboost} package version 2.9-8 (2023-09-06) \insertCite{mboost}{WALS}.
#'
#' @returns \code{walsNB.default()} raises an error if \code{x} is not an object
#' of class \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise
#' returns an object of class \code{"walsNBmatrix"}. See above for more details.
#'
#' @export
walsNB.default <- function(x, ...) {
  # inspired by glmboost.default in mboost.
  if (extends(class(x), "matrix")) {
    return(walsNB.matrix(x, ...))
  }
  stop("No method for objects of class ", sQuote(class(x)), " implemented.")
}



#' Fitter function for Weighted Average Least Squares estimation of NB2 regression model
#'
#' Workhorse function behind \code{\link{walsNB}} and used internally in
#' \code{\link{walsNBfitIterate}}.
#'
#' @inheritParams walsGLMfit
#' @param y Count response as vector.
#' @param rhoStart Starting value for log-dispersion parameter of NB2
#' @param family Object of class \code{"\link[WALS]{familyNBWALS}"}. Currently only supports
#' \code{\link[WALS]{negbinWALS}}.
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
#' to compute the eigendecomposition of \eqn{\bar{\Xi}} instead of \code{eigen()}.
#' In this case, the tolerances of \code{svdTol} and \code{svdRtol} are used to
#' determine whether \eqn{\bar{\Xi}} is of full rank (need it for \eqn{\bar{\Xi}^{-1/2}}).
#' @param ... Arguments for internal function \code{\link[WALS]{computePosterior}}.
#'
#' @details The method to be specified in \code{method} mainly differ in the way
#' they compute the fully restricted and unrestricted estimators for the
#' transformed regressors \eqn{Z}, i.e. \eqn{\tilde{\gamma}_{1r}},
#' and \eqn{\tilde{\gamma}_{u}}.
#'
#' \describe{
#' \item{"fullSVD"}{Recommended approach. First applies an SVD to \eqn{\bar{Z}_{1}}
#' to compute \eqn{\bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2}}:
#' It is used for computing the inverse of
#'
#' \deqn{\bar{X}_{1}^{\top}\bar{X}_{1}
#'  + \bar{g} \bar{\epsilon} X_{1}^{\top}\bar{q} \bar{q}^{\top} X_{1},}
#'
#' when using the Sherman-Morrison-Woodbury formula. We further leverage the
#' SVD of \eqn{\bar{Z}_1} and additionally \eqn{\bar{Z}} to compute the
#' unrestricted estimator \eqn{\tilde{\gamma}_{u}} and the fully restricted
#' estimator \eqn{\tilde{\gamma}_{r}}. For \eqn{\tilde{\gamma}_{u}}, we simply
#' use the SVD of \eqn{\bar{Z}} to solve the full equation system derived from
#' the one-step ML problem for more details. The SVD of \eqn{\bar{Z}_1} is further
#' used in computing the model averaged estimator for the focus regressors
#' \eqn{\hat{\gamma}_1}.
#'
#' Described in more detail in the appendix of \insertCite{huynhwals;textual}{WALS}.}
#'
#'
#' \item{"original"}{Computes all inverses directly using \code{\link[base]{solve}}
#' and does not make use of the Sherman-Morrison-Woodbury formula for certain
#' inverses. Specifically, it directly inverts the matrix
#' \eqn{\bar{Z}_{1}^{\top} \bar{Z}_{1}} using \code{\link[base]{solve}}
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
#' All variables in the code that contain "start" in their name are computed
#' using the starting values of the one-step ML estimators. See section
#' "One-step ML estimator" of \insertCite{huynhwalsnb}{WALS} for details.
#'
#' @returns A list containing
#' \item{coef}{Model averaged estimates of all coefficients.}
#' \item{beta1}{Model averaged estimates of the coefficients of the focus regressors.}
#' \item{beta2}{Model averaged estimates of the coefficients of the auxiliary regressors.}
#' \item{rho}{Model averaged estimate of the log-dispersion parameter of the
#' NB2 distribution.}
#' \item{gamma1}{Model averaged estimates of the coefficients of the transformed
#' focus regressors.}
#' \item{gamma2}{Model averaged estimates of the coefficients of the transformed
#' auxiliary regressors.}
#' \item{condition}{Condition number of the matrix
#' \eqn{\bar{\Xi} = \bar{\Delta}_{2} \bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2} \bar{\Delta}_{2}}.}
#' \item{vcovBeta}{\code{NULL}, not implemented yet, placeholder for estimated
#' covariance matrix of the regression coefficients.}
#' \item{vcovGamma}{\code{NULL}, not implemented yet, placeholder for estimated
#' covariance matrix of the coefficients of the transformed regressors.}
#' \item{betaStart}{Starting values of the regression coefficients for the
#' one-step ML estimators.}
#' \item{rhoStart}{Starting values of the dispersion parameter for the
#' one-step ML estimators.}
#' \item{method}{Stores \code{method} used from the arguments.}
#' \item{prior}{\code{familyPrior}. The \code{prior} specified in the arguments.}
#' \item{betaUn1}{If \code{keepUn = TRUE}, contains the unrestricted one-step ML
#' estimators of the coefficients of the focus regressors. Else \code{NULL}.}
#' \item{betaUn2}{If \code{keepUn = TRUE}, contains the unrestricted one-step ML
#' estimators of the coefficients of the auxiliary regressors. Else \code{NULL}.}
#' \item{gammaUn1}{If \code{keepUn = TRUE}, contains the unrestricted one-step ML
#' estimators of the coefficients of the transformed focus regressors. Else \code{NULL}.}
#' \item{gammaUn2}{If \code{keepUn = TRUE}, contains the unrestricted one-step ML
#' estimators of the coefficients of the transformed auxiliary regressors. Else \code{NULL}.}
#' \item{gamma1r}{If \code{keepR = TRUE}, contains the fully restricted one-step
#' ML estimator for the transformed regressors (only focus regressors).
#' Else \code{NULL}.}
#' \item{k1}{Number of focus regressors.}
#' \item{k2}{Number of auxiliary regressors.}
#' \item{n}{Number of observations.}
#' \item{X1names}{Names of the focus regressors.}
#' \item{X2names}{Names of the auxiliary regressors.}
#' \item{familyStart}{The family object of class \code{"\link{familyNBWALS}"} used for the
#' estimation of the starting values.}
#' \item{family}{The family object of class \code{"\link{familyNBWALS}"} used later for predictions.}
#' \item{fitted.link}{Linear link fitted to the data.}
#' \item{fitted.values}{Estimated conditional mean for the data. Lives on the
#' scale of the response.}
#'
#' @seealso [walsNB], [walsNBfitIterate].
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data("NMES1988", package = "AER")
#' NMES1988 <- na.omit(NMES1988)
#' form <- (visits ~ health + chronic + age + insurance + adl + region + gender
#'          + married + income + school + employed)
#' X <- model.matrix(form, data = NMES1988)
#' focus <- c("(Intercept)", "healthpoor", "healthexcellent", "chronic", "age",
#'         "insuranceyes")
#' aux <- c("adllimited", "regionnortheast", "regionmidwest", "regionwest",
#'          "gendermale", "marriedyes", "income", "school", "employedyes")
#' X1 <- X[, focus]
#' X2 <- X[, aux]
#' y <- NMES1988$visits
#'
#' # starting values from glm.nb() from MASS
#' startFit <- MASS::glm.nb(y ~ X[,-1])
#' betaStart <- coef(startFit)
#' rhoStart <- startFit$theta
#' k1 <- ncol(X1)
#' k2 <- ncol(X2)
#'
#' str(walsNBfit(X1, X2, y, rhoStart, family = negbinWALS(scale = rhoStart, link = "log"),
#'               betaStart1 = betaStart[1:k1],
#'               betaStart2 = betaStart[(k1 + 1):(k1 + k2)],
#'               prior = weibull(), method = "fullSVD"))
#'
#' @export
walsNBfit <- function(X1, X2, y, betaStart1, betaStart2, rhoStart, family,
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
  etaStart <- drop(X1 %*% betaStart1 + X2 %*% betaStart2)
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
  qStart <- family$q(etaStart, y)
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
              gamma1 = gamma1, gamma2 = gamma2, condition = outSemiOrt$condition,
              vcovBeta = NULL, vcovGamma = NULL,
              betaStart = c(betaStart1, betaStart2), rhoStart = rhoStart,
              method = method, prior = prior,
              betaUn1 = if (keepUn) betaUn1 else NULL,
              betaUn2 = if (keepUn) betaUn2 else NULL,
              gammaUn1 = if (keepUn) gammaUn1 else NULL,
              gammaUn2 = if (keepUn) gammaUn2 else NULL,
              gamma1r = if (keepR) gamma1r else NULL,
              k1 = k1, k2 = k2, n = n, X1names = X1names, X2names = X2names)

  # contains previous rho, used for obtaining betaStart etc.
  fit$familyStart <- family
  # new family with model averaged rho
  fit$family <- negbinWALS(rho, link = family$link)
  fit$fitted.link <- drop(X1 %*% fit$beta1 + X2 %*% fit$beta2)
  fit$fitted.values <- family$linkinv(fit$fitted.link)

  # assign names to variables
  names(fit$coef) <- Xnames
  names(fit$betaStart) <- Xnames
  names(fit$beta1) <- names(fit$gamma1) <- X1names
  names(fit$beta2) <- names(fit$gamma2) <- X2names
  return(fit)
}

#' Iteratively fitting walsNB, internal function for walsNB.formula and
#' walsNB.matrix.
#'
#' Wrapper around \code{\link[WALS]{walsNBfit}} that allows iteratively
#' (re-)fitting \code{\link[WALS]{walsNB}} models.
#'
#' @param y Count response as vector.
#' @inheritParams walsGLMfitIterate
#' @param link specifies the link function, currently only "log" is supported.
#' @param controlInitNB Controls estimation of starting values for one-step ML,
#' see \code{\link[WALS]{controlNB}}.
#' @param tol Only used if \code{iterate = TRUE} and \code{nIt = NULL}. If the
#' Euclidean distance between the previous and current coefficient vector divided
#' by the square root of the length of the vector falls below \code{tol} and the
#' absolute difference between the previous and current dispersion parameter
#' falls below \code{tol}, then the algorithm stops. See below for more details.
#' @param ... Arguments to be passed to the workhorse function \code{\link[WALS]{walsNBfit}}.
#'
#' @returns A list containing all elements returned from \code{\link[WALS]{walsNBfit}}
#' and additionally the following elements:
#' \item{y}{If \code{keepY = TRUE}, contains the response vector.}
#' \item{x}{list. If \code{keepX = TRUE}, then it is a list with elements
#' \code{x1} and \code{x2} containing the design matrices of the focus and
#' auxiliary regressors, respectively.}
#' \item{initialFit}{List containing information (e.g. convergence) on the
#' estimation of the starting values for \code{\link[WALS]{walsNBfit}}.
#' See return of \code{\link[WALS]{fitNB2}} for more information.}
#' \item{weights}{returns the argument \code{weights}.}
#' \item{offset}{returns the argument \code{offset}.}
#' \item{converged}{Logical. Only relevant if \code{iterate = TRUE}. Equals
#' \code{TRUE} if iterative fitting converged, else \code{FALSE}. Is \code{NULL}
#' if \code{iterate = FALSE}.}
#' \item{it}{Number of iterations run in the iterative fitting algorithm.
#' \code{NULL} if \code{iterate = FALSE}.}
#' \item{deviance}{Deviance of the fitted (conditional) NB2 regression model.}
#' \item{residuals}{Raw residuals, i.e. response - fitted mean.}
#'
#' @seealso [walsNB], [walsNBfit].
#'
#' @details
#' The parameter \code{tol} is used to control the convergence of the iterative
#' fitting algorithm. Let \eqn{i} be the current iteration step for the
#' coefficient vector \eqn{\beta_{i} = (\beta_{i,1}, \ldots, \beta_{i,k})'},
#' \eqn{k > 0}, and dispersion parameter \eqn{\rho_{i}}. If
#' \deqn{\frac{||\beta_{i} - \beta_{i-1}||_{2}}{\sqrt{k}}
#' = \sqrt{\frac{\sum_{j = 1}^{k} (\beta_{i,j} - \beta_{i-1,j})^{2}}{k}} < \texttt{tol},}
#' and
#' \deqn{|\rho_{i} - \rho_{i-1}| < \texttt{tol},}
#' then the fitting process is assumed to have converged and stops.
#'
#'
#' @examples
#' data("NMES1988", package = "AER")
#' NMES1988 <- na.omit(NMES1988)
#' form <- (visits ~ health + chronic + age + insurance + adl + region + gender
#'          + married + income + school + employed)
#' X <- model.matrix(form, data = NMES1988)
#' focus <- c("(Intercept)", "healthpoor", "healthexcellent", "chronic", "age",
#'         "insuranceyes")
#' aux <- c("adllimited", "regionnortheast", "regionmidwest", "regionwest",
#'          "gendermale", "marriedyes", "income", "school", "employedyes")
#' X1 <- X[, focus]
#' X2 <- X[, aux]
#' y <- NMES1988$visits
#'
#' str(walsNBfitIterate(y, X1, X2, prior = weibull(), link = "log",
#'                      method = "fullSVD", iterate = TRUE))
#'
#' @export
walsNBfitIterate <- function(y, X1, X2, link = "log", na.action = NULL,
                              weights = NULL, offset = NULL,
                              prior = weibull(), controlInitNB = controlNB(),
                              keepY = TRUE, keepX = FALSE,
                              iterate = FALSE, tol = 1e-6, maxIt = 50, nIt = NULL,
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
    warning("Iteration limit reached in initial NB2 fit of full model. ",
            "Try raising maxIt in controlOptim or try using ",
            "controlInitNB = controlNB(initMASS = TRUE). ",
            "See component $initialFit for more details.\n")
  } else if (nb2$convergence == 10) {
    warning("Degeneracy in Nelder-Mead simplex in initial NB2 fit of full model. ",
            "Try using controlInitNB = controlNB(initMASS = TRUE) ",
            "See component $initialFit for more details.\n")
  } else if (nb2$convergence == 51) {
    warning("Warning in L-BFGS-B in initial NB2 fit of full model. ",
            "Try using controlNB(initMASS = TRUE) ",
            "See component $initialFit$message for more details.\n")
  } else if (nb2$convergence == 52) {
    stop("Error in L-BFGS-B in initial NB2 fit of full model. ",
         "No model fitted.", paste0("Error message: ", nb2$message),
         " Try using controlInitNB = controlNB(initMASS = TRUE). \n")
  } else if (nb2$convergence == 99) {
    stop("Convergence issue in IWLS algo in glm.nb() for initial NB2 fit of full model. ",
         "No model fitted. Try using controlInitNB = controlNB(initMASS = TRUE). ",
         "See ?controlNB() and ?glm.nb() for more info.\n")
  } else if (nb2$convergence != 0) {
    warning("Unknown non-zero exit for optim in initial NB2 fit of full model. ",
            "See initial fit for more details.\n")
  }

  betaStart <- if (controlInitNB$restricted) {
    c(nb2$coefficients, rep(0, k2))
  } else nb2$coefficients

  betaStart1 <- betaStart[1L:k1]
  betaStart2 <- betaStart[(k1 + 1L):(k1 + k2)]
  rhoStart <- nb2$theta

  # reuse iterative code by setting nIt = 1
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
    out <- walsNBfit(X1 = X1, X2 = X2, y = y, betaStart1 = betaCurrent[1L:k1],
                      betaStart2 = betaCurrent[(k1 + 1L):(k1 + k2)],
                      rhoStart = rhoCurrent,
                      family = family, prior = prior,
                      ...)

    betaCurrent <- out$coef
    rhoCurrent <- out$rho
    it <- it + 1

    if (verbose) cat(paste("\rfinished iteration", i))

    if (is.null(nIt) &&
        ((norm(betaOld - betaCurrent, type = "2") / sqrt(length(betaCurrent))) < tol
         && abs(rhoOld - rhoCurrent) < tol)
    ) {
      converged <- TRUE
      if (verbose) cat("\nalgorithm converged\n")
      break
    }

  }

  if (!is.null(nIt)) {
    converged <- NULL
  } else if (!converged) warning("algorithm failed to converge\n")


  # replace starting values with original starting values
  out$betaStart <- betaStart
  out$rhoStart <- rhoStart
  out$familyStart <- familyStart

  # add more elements
  if (keepY) out$y <- y
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$initialFit <- nb2
  out$weights <- weights
  out$offset <- offset
  out$converged <- if (iterate) converged else NULL
  out$it <- if (iterate) it else NULL

  # deviance & residuals
  wt <- if (is.null(weights)) rep(1, nrow(X1)) else weights
  mu <- out$fitted.values
  out$deviance <- sum(family$dev.resids(out$y, mu, wt))
  out$residuals <- out$y - mu

  return(out)
}



## Class methods ---------------------------------------------------------------

#' Calculate Variance-Covariance Matrix for a \code{"walsNB"} object
#'
#' This method always raises an error because the covariance matrix of the
#' walsNB estimator has not been derived yet.
#'
#' @param object An object of class \code{"walsNB"}.
#' @param ... For expansion in the future.
#'
#' @returns No return value, only raises error because no covariance matrix
#' estimator exists yet.
#'
#' @export
vcov.walsNB <- function(object, ...) {
  stop("No method for objects of class ", sQuote(class(object)[1]), " implemented.")
}

#' @rdname predict.walsGLM
#'
#' @returns \code{summary.walsNB()} returns an object of class
#' \code{"summary.walsNB"} which contains the necessary fields for printing the
#' summary in \code{print.summary.walsNB()}.
#'
#' @export
summary.walsNB <- function(object, ...) {
  object <- summary.wals(object, ...)

  # remove SE estimation from walsGLM
  object$focusCoefs <- object$focusCoefs[,"Estimate", drop = FALSE]
  object$auxCoefs <- object$auxCoefs[,"Estimate", drop = FALSE]

  class(object) <- c("summary.walsNB")
  return(object)
}

#' @rdname predict.walsGLM
#'
#' @returns \code{print.summary.walsNB()} invisibly returns its input argument
#' \code{x}, i.e. an object of object of class \code{"summary.walsNB"}.
#'
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

#' Control function for initial NB fit
#'
#' Defines controllable parameters of initial NB fit in \code{\link[WALS]{walsNB}}.
#'
#' @param start Optional starting values for \code{\link[WALS]{fitNB2}}. Only used if
#' \code{initMASS = FALSE}.
#' @param method Optimization method used in \code{\link[stats]{optim}}. Only used if
#' \code{initMASS = FALSE}.
#' @param controlOptim List with parameters controlling optimization process of
#' \code{\link[stats]{optim}}. Only used if \code{initMASS = FALSE}.
#' @param initThetaMASS If TRUE, then initial \eqn{\log{\theta}} of
#' \code{\link[WALS]{fitNB2}} is estimated using \code{\link[MASS]{theta.ml}}
#' (ML-estimation over 1 variable) based on regression coefficients from
#' Poisson regression. If \code{FALSE}, then initial \eqn{\log{\theta}} = 0 is used.
#' @param initMASS If \code{TRUE} (default), then initial fit in \code{\link[WALS]{fitNB2}}
#' is estimated via \code{\link[MASS]{glm.nb}} and \code{initThetaMASS} is ignored.
#' If \code{FALSE}, then the initial fit is estimated by minimizing the
#' log-likelihood using \code{\link[stats]{optim}}.
#' @param restricted If \code{TRUE}, then initial fit in \code{\link[WALS]{fitNB2}}
#' only considers the focus regressors. By default \code{FALSE}, then the unrestricted
#' model is estimated in \code{\link[WALS]{fitNB2}} (i.e. all regressors).
#' @param eps Controls argument \code{eps} in \code{\link[WALS]{fitNB2}} for generating
#' starting value for \code{logTheta} (\eqn{\log{\theta}}) via \code{\link[MASS]{theta.ml}}.
#' @param epsilonMASS Sets epsilon in control argument of \code{\link[MASS]{glm.nb}}.
#'
#' @returns Returns a list containing the parameters specified in the arguments
#' to be used in \code{\link[WALS]{walsNB}} (and \code{\link[WALS]{walsNBfitIterate}}).
#'
#' @examples
#' data("NMES1988", package = "AER")
#' walsNB(visits ~ health + chronic + age + gender | I((age^2)/10) +
#'        married + region, data = NMES1988, prior = weibull(),
#'        controlInitNB = controlNB(initMASS = FALSE, restricted = TRUE))
#'
#' @seealso [walsNB], [walsNBfitIterate].
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

