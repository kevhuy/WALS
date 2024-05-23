#' Weighted Average Least Squares for Generalized Linear Models
#'
#' Performs model averaging of generalized linear models (GLMs) using the
#' Weighted-Average Least Squares method described in \insertCite{deluca2018glm;textual}{WALS}.
#'
#' @details
#' Computes WALS estimates when focus regressors (X1) are present in all
#' submodels and model averaging takes place over the auxiliary regressors (X2).
#'
#' @references
#' \insertAllCited{}
#'
#' @export
walsGLM <- function(x, ...) UseMethod("walsGLM", x)

#' \code{walsGLM.formula()} uses formulas to specify the design matrix.
#' @rdname walsGLM
#'
#' @inheritParams wals.formula
#' @param family Object of class \code{"\link[WALS]{familyWALS}"}.
#' @inheritParams walsGLMfitIterate
#' @param tol Only used if \code{iterate = TRUE} and \code{nIt = NULL}.
#' If the Euclidean distance between the previous and current coefficient vector
#' divided by the square root of the length of the vector falls below \code{tol},
#' then the algorithm stops. See \code{\link[WALS]{walsGLMfitIterate}} for more details.
#' @param nIt Only used if \code{iterate = TRUE}. If this is specified, then
#' \code{tol} is ignored and the algorithm iterates \code{nIt} times. This option
#' should not be used unless the user has a specific reason to run the algorithm
#' \code{nIt} times, e.g. for replication purposes.
#' @param ... Arguments for workhorse \code{\link[WALS]{walsGLMfit}}.
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
#' **WARNING:** Interactions in formula do work work properly yet.
#' It is recommended to manually create the interactions beforehand and then
#' to insert them as 'linear terms' in the formula.
#'
#' @returns \code{walsGLM.formula()} returns an object of class \code{"walsGLM"}
#' which inherits from \code{"\link[WALS]{wals}"}. This is a list that contains
#' all elements returned from \code{\link[WALS]{walsGLMfitIterate}} and additionally
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
#' See returns of \code{\link[WALS]{walsGLMfit}} and \code{\link[WALS]{walsGLMfitIterate}}
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
                            prior = weibull(), controlInitGLM = controlGLM(),
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
                           controlInitGLM, keepY, keepX, iterate, tol, maxIt,
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
#' @aliases walsGLMmatrix
#'
#' @inheritParams wals.matrix
#'
#' @returns \code{walsGLM.matrix()} returns an object of class
#' \code{"walsGLMmatrix"}, which inherits from \code{"walsGLM"}, \code{"walsMatrix"}
#' and \code{"wals"}. This is a list that contains all elements returned from
#' \code{\link[WALS]{walsGLMfitIterate}} and additionally the call in \code{cl}.
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
                           prior = weibull(), controlInitGLM = controlGLM(),
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
                           controlInitGLM, keepY, keepX, iterate, tol, maxIt,
                           nIt, verbose, ...)

  out$call <- cl

  class(out) <- c("walsGLMmatrix", "walsGLM", "walsMatrix", "wals")
  return(out)
}

#' @rdname walsGLM
#'
#' @details
#' \code{walsGLM.default()} raises an error if \code{x} is not an object of class
#' \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise it calls
#' \code{walsGLM.matrix()}. It is a modified version of \code{glmboost.default}
#' from the \code{mboost} package version 2.9-8 (2023-09-06) \insertCite{mboost}{WALS}.
#'
#' @returns \code{walsGLM.default()} raises an error if \code{x} is not an object
#' of class \code{"matrix"} or a class that extends \code{"matrix"}. Otherwise
#' returns an object of class \code{"walsGLMmatrix"}. See above for more details.
#'
#'
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
#' Workhorse function behind \code{\link[WALS]{walsGLM}} and used internally in
#' \code{\link[WALS]{walsGLMfitIterate}}.
#'
#' @inheritParams walsFit
#' @param betaStart1 Starting values for coefficients of focus regressors X1.
#' @param betaStart2 Starting values for coefficients of auxiliary regressors X2.
#' @param family Object of class \code{"\link[WALS]{familyWALS}"}.
#' @param postmult If \code{TRUE} (default), then it computes
#' \deqn{\bar{Z}_{2} = \bar{X}_{2} \bar{\Delta}_{2} \bar{T} \bar{\Lambda}^{-1/2} \bar{T}^{\top},}
#' where \eqn{\bar{T}} contains the eigenvectors and \eqn{\bar{\Lambda}} the
#' eigenvalues from the eigenvalue decomposition
#' \deqn{\bar{\Xi} = \bar{T} \bar{\Lambda} \bar{T}^{\top},}
#' instead of
#' \deqn{\bar{Z}_{2} = \bar{X}_{2} \bar{\Delta}_{2} \bar{T} \bar{\Lambda}^{-1/2}.}
#' See \insertCite{huynhwals;textual}{WALS} for more details. The latter is used
#' in the original MATLAB code for WALS in the linear regression model,
#' see eq. (12) of \insertCite{magnus2016wals;textual}{WALS}.
#' The first form is required in eq. (9) of \insertCite{deluca2018glm;textual}{WALS}.
#' **Thus, it is not recommended to set \code{postmult = FALSE}.**
#' @param ... Further arguments passed to \code{\link[WALS]{walsFit}}.
#'
#' @details
#' Uses \code{\link[WALS]{walsFit}} under the hood after transforming the regressors
#' \code{X1} and \code{X2} and the response \code{y}. For more details, see
#' \insertCite{huynhwals}{WALS} and \insertCite{deluca2018glm;textual}{WALS}.
#'
#' @returns A list containing all elements returned by \code{\link[WALS]{walsFit}},
#' except for \code{residuals}, and additionally (some fields are replaced)
#' \item{condition}{Condition number of the matrix
#' \eqn{\bar{\Xi} = \bar{\Delta}_{2} \bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2} \bar{\Delta}_{2}}.}
#' \item{family}{Object of class \code{"\link[WALS]{familyWALS}"}. The family used.}
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
#' @examples
#' data("HMDA", package = "AER")
#' X <- model.matrix(deny ~ pirat + hirat + lvrat + chist + mhist + phist + selfemp + afam,
#'                   data = HMDA)
#' X1 <- X[,c("(Intercept)", "pirat", "hirat", "lvrat", "chist2", "chist3",
#'         "chist4", "chist5", "chist6", "mhist2", "mhist3", "mhist4", "phistyes")]
#' X2 <- X[,c("selfempyes", "afamyes")]
#' y <- HMDA$deny
#'
#' # starting values from glm.fit()
#' betaStart <- glm.fit(X, y, family = binomialWALS())$coefficients
#' k1 <- ncol(X1)
#' k2 <- ncol(X2)
#'
#' str(walsGLMfit(X1, X2, y,
#'                betaStart1 = betaStart[1:k1],
#'                betaStart2 = betaStart[(k1 + 1):(k1 + k2)],
#'                family = binomialWALS(), prior = weibull()))
#'
#'
#' @export
walsGLMfit <- function(X1, X2, y, betaStart1, betaStart2,
                       family, prior = weibull(), postmult = TRUE, ...) {
  X1names <- colnames(X1)
  X2names <- colnames(X2)
  Xnames <- c(X1names, X2names)

  etaStart <- drop(X1 %*% betaStart1 + X2 %*% betaStart2)
  X1start <- family$transformX(X1, etaStart, y)
  X2start <- family$transformX(X2, etaStart, y)
  yStart <- family$transformY(y, X1start, X2start, betaStart1, betaStart2, etaStart)


  # use generic WALS algo for linear models
  fit <- walsFit(X1 = X1start, X2 = X2start, y = yStart, sigma = 1,
                  prior = prior, prescale = TRUE, postmult = postmult, ...)

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
#' Wrapper around \code{\link[WALS]{walsGLMfit}} that allows iteratively
#' (re-)fitting \code{\link[WALS]{walsGLM}} models.
#'
#' @inheritParams walsGLMfit
#' @param na.action Not implemented yet.
#' @param weights Not implemented yet.
#' @param offset Not implemented yet.
#' @param controlInitGLM Controls estimation of starting values for one-step ML,
#' see \code{\link[WALS]{controlGLM}}.
#' @param keepY If \code{TRUE}, then output keeps response.
#' @param keepX If \code{TRUE}, then output keeps the design matrices.
#' @param iterate if \code{TRUE} then the WALS algorithm is iterated using the previous
#' estimates as starting values.
#' @param tol Only used if \code{iterate = TRUE} and \code{nIt = NULL}.
#' If the Euclidean distance between the previous and current coefficient vector
#' divided by the square root of the length of the vector falls below \code{tol},
#' then the algorithm stops. See below for more details.
#' @param maxIt Only used if \code{iterate = TRUE} and \code{nIt = NULL}. Aborts
#' iterative fitting when number of iterations exceed \code{maxIt}.
#' @param nIt Only used if \code{iterate = TRUE}. If this is specified, then
#' \code{tol} is ignored and the algorithm iterates \code{nIt} times.
#' @param verbose If \code{verbose = TRUE}, then it prints the iteration process
#' (only relevant if \code{iterate = TRUE}).
#' @param ... Arguments to be passed to the workhorse function \code{\link[WALS]{walsGLMfit}}.
#'
#' @returns A list containing all elements returned from \code{\link[WALS]{walsGLMfit}}
#' and additionally the following elements:
#' \item{y}{If \code{keepY = TRUE}, contains the response vector.}
#' \item{x}{list. If \code{keepX = TRUE}, then it is a list with elements
#' \code{x1} and \code{x2} containing the design matrices of the focus and
#' auxiliary regressors, respectively.}
#' \item{initialFit}{List containing information (e.g. convergence) on the
#' estimation of the starting values for \code{\link[WALS]{walsGLMfit}}.
#' See \code{\link[stats]{glm.fit}} for more information.}
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
#' @details
#' The parameter \code{tol} is used to control the convergence of the iterative
#' fitting algorithm. Let \eqn{i} be the current iteration step for the
#' coefficient vector \eqn{\beta_{i} = (\beta_{i,1}, \ldots, \beta_{i,k})', k > 0}.
#' If
#' \deqn{\frac{||\beta_{i} - \beta_{i-1}||_{2}}{\sqrt{k}}
#' = \sqrt{\frac{\sum_{j = 1}^{k} (\beta_{i,j} - \beta_{i-1,j})^{2}}{k}} < \texttt{tol},}
#' then the fitting process is assumed to have converged and stops.
#'
#' @seealso [walsGLM], [walsGLMfit].
#'
#' @examples
#' data("HMDA", package = "AER")
#' X <- model.matrix(deny ~ pirat + hirat + lvrat + chist + mhist + phist + selfemp + afam,
#'                   data = HMDA)
#' X1 <- X[,c("(Intercept)", "pirat", "hirat", "lvrat", "chist2", "chist3",
#'         "chist4", "chist5", "chist6", "mhist2", "mhist3", "mhist4", "phistyes")]
#' X2 <- X[,c("selfempyes", "afamyes")]
#' y <- HMDA$deny
#'
#' str(walsGLMfitIterate(y, X1, X2, family = binomialWALS(), prior = weibull(),
#'                       iterate = TRUE))
#'
#' @export
walsGLMfitIterate <- function(y, X1, X2, family, na.action = NULL,
                              weights = NULL, offset = NULL,
                              prior = weibull(), controlInitGLM = controlGLM(),
                              keepY = TRUE, keepX = FALSE, iterate = FALSE,
                              tol = 1e-6, maxIt = 50, nIt = NULL,
                              verbose = FALSE, ...) {
  # Useful quantities
  k1 <- ncol(X1)
  k2 <- ncol(X2)

  # check if X1 and X2 contain the same variables
  if (any(colnames(X1) %in% colnames(X2))) stop("X1 and X2 contain the same variables")

  # generate starting values
  Xinit <- if (controlInitGLM$restricted) X1 else cbind(X1, X2)
  initialFit <- glm.fit(Xinit, y, family = family,
                        control = controlInitGLM$controlGLMfit)

  if (!initialFit$converged) {
    warning("Convergence issue in IWLS algo in glm.fit for initial fit of full model. ",
            "See initial fit for more details.")

  }

  if (controlInitGLM$restricted) {
    betaStart <- c(coef(initialFit), rep(0, k2))
    names(betaStart) <- c(colnames(X1), colnames(X2))
  } else {
    betaStart <- coef(initialFit)
  }

  # reuse iterative code by setting nIt = 1
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

    if (verbose) cat(paste("\rfinished iteration", it))

    if (is.null(nIt)
        && ((norm(betaOld - betaCurrent, type = "2") / sqrt(length(betaCurrent))) < tol)
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

  # add more elements
  if (keepY) out$y <- family$initializeY(y) # e.g. convert logical to 0s and 1s.
  if (keepX) out$x <- list(focus = X1, aux = X2)
  out$initialFit <- initialFit
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
#' Methods for walsGLM, walsGLMmatrix, walsNB and walsNBmatrix Objects
#'
#' Methods for extracting information from fitted model-averaging objects of
#' classes \code{"walsGLM"}, \code{"walsGLMmatrix"}, \code{"walsNB"} and
#' \code{"walsNBmatrix"}.
#'
#' @param object,x An object of class \code{"walsGLM"}, \code{"walsGLMmatrix"},
#' \code{"walsNB"}, \code{"walsNBmatrix"}, \code{"summary.walsGLM"} or
#' \code{"summary.walsNB"}.
#' @inheritParams predict.wals
#' @param type Character specifying the type of prediction, residual or model
#' part to be returned. For details see below.
#' @param at  Optional. Only available if a family of class \code{"\link[WALS]{familyWALScount}"}
#' was used for fitting. If \code{type = "prob"}, a numeric vector at which
#' the probabilities are evaluated. By default \code{0:max(y)} is used
#' where \code{y} is the original observed response.
#' @param log Logical. If \code{TRUE}, then returns the log-density. If
#' \code{FALSE} (default), then returns density. Only relevant if
#' \code{type = "density"}.
#'
#' @details
#' As the \code{"-matrix"} classes \code{"walsGLMmatrix"} and \code{"walsNBmatrix"}
#' inherit from the "non-matrix" classes, i.e. \code{"walsGLM"} and \code{"walsNB"},
#' respectively, the following text will treat them as equivalent because
#' they inherit all methods but \code{predict} from their "non-matrix" versions.
#' Thus, when \code{"walsGLM"} or \code{"walsNB"} are mentioned, we also refer to
#' their \code{"-matrix"} versions, except when explicitly stated. Moreover,
#' note that \code{"walsNB"} and \code{"walsNBmatrix"} inherit most methods from
#' \code{"walsGLM"} and \code{"walsGLMmatrix"}.
#'
#' A set of standard extractor functions for fitted model objects is available
#' for objects of class \code{"walsGLM"} and \code{"walsNB"}, including methods to
#' the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}
#' which print the model-averaged estimation of the coefficients along with some
#' further information.
#'
#' The \code{\link[base]{summary}} methods returns an object of
#' class \code{"summary.walsGLM"} for objects of class \code{"walsGLM"} and an
#' object of class \code{"summary.walsNB"} for objects of class \code{"walsNB"}.
#' They contain the relevant summary statistics which can then be printed using
#' the associated \code{print()} methods.
#' Inspired by \insertCite{deluca2011stata;textual}{WALS}, the summary statistics
#' also show \code{Kappa} which is an indicator for the numerical stability of
#' the method, i.e. it shows the square root of the condition number of the
#' matrix \eqn{\bar{\Xi} = \bar{\Delta}_{2} \bar{X}_{2}^{\top} \bar{M}_{1}
#' \bar{X}_{2} \bar{\Delta}_{2}}. The summary further shows the deviance and
#' provides information on the prior and family used.
#'
#' A \code{\link[stats]{logLik}} method is provided that returns the log-likelihood
#' given the family used and the model-averaged estimates of the coefficients.
#'
#' \code{"walsGLM"} inherits from \code{"wals"}, while \code{"walsNB"} inherits from
#' both, \code{"walsGLM"} and \code{"wals"}. Thus, see \code{\link[WALS]{predict.wals}}
#' for more methods.
#'
#' The \code{\link[stats]{predict}} and \code{\link[stats]{residuals}} methods,
#' especially the different types of predictions/residuals controlled by
#' \code{type}, are inspired by the corresponding methods in \code{countreg}
#' version 0.2-1 (2023-06-13) \insertCite{countreg,countreghurdle}{WALS}, see
#' e.g. \code{predict.hurdle()} from \code{countreg}, and \code{\link[stats]{stats}}
#' version 4.3.1 (2023-06-16) \insertCite{R2023}{WALS}, see e.g.
#' \code{\link[stats]{residuals.glm}}. The \code{summary()}, \code{print.summary()},
#' \code{print()} and \code{logLik()} methods are also inspired by the corresponding
#' methods for objects of class \code{"glm"} in \code{\link[stats]{stats}}, see
#' e.g. \code{\link[stats]{print.summary.glm}}.
#'
#' \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} are inherited from
#' \code{"wals"} (see \code{\link[WALS]{predict.wals}} for more), except for
#' objects of class \code{"walsNB"} (see \code{\link[WALS]{vcov.walsNB}}).
#' The \code{type} argument specifies which part of the coefficient
#' vector/covariance matrix of the estimates should be returned.
#' For \code{type = "all"}, they return the complete vector/matrix.
#' For \code{type = "focus"} and \code{type = "aux"} they return only the part
#' corresponding to the focus and auxiliary regressors, respectively.
#' Additionally, the user can choose whether to return the
#' estimated coefficients/covariance matrix for the original regressors \eqn{X}
#' (\eqn{\beta} coefficients) or of the transformed regressors \eqn{Z}
#' (\eqn{\gamma} coefficients).
#'
#' The extractors \code{\link[stats]{terms}} and \code{\link[stats]{model.matrix}}
#' are also inherited from \code{"wals"}. They only allow \code{type = "focus"}
#' and \code{type = "aux"} and extract the corresponding component of the model.
#'
#' @returns \code{predict.walsGLM()} and \code{predict.walsGLMmatrix()} return
#' different types of predictions depending on the argument \code{type}:
#' * \code{type = "response"}: vector. Predicted mean
#' * \code{type = "link"}: vector. Predicted linear link
#' * \code{type = "variance"}: vector. Predicted variance
#' * \code{type = "prob"}: matrix. Only available if a family of class
#' \code{"\link[WALS]{familyWALScount}"} was used for fitting or for objects of
#' class \code{"walsNB"} or \code{"walsNBmatrix"}. Returns the probability at
#' counts specified by \code{at}.
#' * \code{type = "density"}: vector. Predicted density
#' * \code{type = "logDens"}: vector. For convenience, returns predicted log-density.
#' Equivalent to setting \code{type = "density"} and \code{log = TRUE}.
#'
#' If \code{type = "prob"}, \code{type = "density"} or \code{type = "logDens"},
#' then \code{newdata} must contain the response or \code{newY} must be
#' specified depending on the class of the object.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [walsGLM], [walsNB], [predict.wals].
#'
#' @examples
#' ## Example for walsGLM objects
#' data("HMDA", package = "AER")
#' fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
#'                         selfemp + afam, family = binomialWALS(), data = HMDA,
#'                        prior = weibull())
#' summary(fitBinomial)
#' vcov(fitBinomial, type = "focus")
#' logLik(fitBinomial)
#' predict(fitBinomial, newdata = HMDA[1:10,], type = "response")
#' familyWALS(fitBinomial)
#'
#' ## Example for walsNB objects
#' data("NMES1988", package = "AER")
#'
#' fWals <- (visits ~ chronic + age + I((age^2)/10) + insurance + medicaid |
#'            adl + gender + married + income + school + afam + employed)
#' fitNB <- walsNB(fWals, data = NMES1988, link = "log", prior = weibull(),
#'                 method = "fullSVD")
#' summary(fitNB)
#' coef(fitNB, type = "aux")
#' residuals(fitNB, type = "pearson")
#' predict(fitNB, newdata = NMES1988[1:10,], type = "prob")
#' terms(fitNB, type = "aux")
#'
#' @export
predict.walsGLM <- function(object, newdata,
                            type = c("response", "link", "variance", "prob",
                                     "density", "logDens"),
                            at = NULL,
                            na.action = na.pass, log = FALSE, ...) {
  # TODO: include offsets
  type <- match.arg(type)

  # convenience type
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

#' @rdname predict.walsGLM
#' @param newY Response vector to be used in predictions. Only relevant when
#' \code{type = "prob"}, \code{type = "density"} or \code{type = "logDens"}.
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

#' @rdname predict.walsGLM
#'
#' @returns \code{residuals.walsGLM()} returns different types of residuals
#' depending on the argument \code{type}:
#' * \code{type = "deviance"}: deviance residuals
#' * \code{type = "pearson"}: Pearson residuals (raw residuals scaled by
#' square root of variance function)
#' * \code{type = "response"}: raw residuals (observed - fitted)
#'
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

#' @rdname predict.walsGLM
#'
#' @returns \code{print.walsGLM()} invisibly returns its input argument \code{x},
#' i.e. an object of object of class \code{"walsGLM"}.
#'
#' @export
print.walsGLM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print.wals(x, digits, ...)
  cat(paste0("\nResidual Deviance: ", signif(x$deviance, digits), "\n"))
  invisible(x)
}

#' @rdname predict.walsGLM
#'
#' @returns \code{summary.walsGLM()} returns an object of class
#' \code{"summary.walsGLM"} which contains the necessary fields for printing the
#' summary in \code{print.summary.walsGLM()}.
#'
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

#' @rdname predict.walsGLM
#'
#' @returns \code{print.summary.walsGLM()} invisibly returns its input argument
#' \code{x}, i.e. an object of object of class \code{"summary.walsGLM"}.
#'
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

#' @rdname predict.walsGLM
#'
#' @returns \code{logLik.walsGLM()} returns the log-likelihood of the fitted
#' model.
#'
#' @export
logLik.walsGLM <- function(object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  y <- residuals(object, type = "response") + fitted(object)
  return(sum(object$family$density(y, object$fitted.link, log = TRUE)))
}

#' @rdname familyWALS
#' @export
familyWALS.walsGLM <- function(object, ...) return(object$family)


## Helper functions ------------------------------------------------------------

#' Control function for initial GLM fit
#'
#' Defines controllable parameters of initial GLM fit in \code{\link[WALS]{walsGLM}}.
#'
#' @param restricted If \code{TRUE}, then initial fit in \code{\link[stats]{glm.fit}}
#' only considers the focus regressors. By default \code{FALSE}, then the unrestricted
#' model is estimated in \code{\link[stats]{glm.fit}} (i.e. all regressors).
#' @param controlGLMfit List. Arguments to be passed to \code{control} argument
#' of \code{\link[stats]{glm.fit}}. See also \code{\link[stats]{glm.control}}.
#'
#' @returns Returns a list containing the parameters specified in the arguments
#' to be used in \code{\link[WALS]{walsGLM}} (and \code{\link[WALS]{walsGLMfitIterate}}).
#'
#' @examples
#' data("HMDA", package = "AER")
#' fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
#'                        selfemp + afam, data = HMDA,
#'                        family = binomialWALS(),
#'                        prior = weibull(),
#'                        controlInitGLM = controlGLM(restricted = TRUE,
#'                                                    controlGLMfit = list(trace = TRUE)))
#'
#' @seealso [walsGLM], [walsGLMfitIterate], [glm.fit], [glm.control].
#'
#' @export
controlGLM <- function(restricted = FALSE, controlGLMfit = list()) {
  return(list(restricted = restricted, controlGLMfit = controlGLMfit))
}
