## This script contains helper functions

# Transformations for one-step estimator of GLMs (DeLuca et al. 2018, Eq. 5)----

# Xbar
# Does nothing with etaStart, just passes it. The argument etaStart is included
# so that it has the same arguments as negbinFixedWALS()$transformX
# (see in familyWALS.R).
transformX <- function(X, etaStart, psiBar) sqrt(as.vector(psiBar)) * X

# muBar
transformY <- function(y, X1bar, X2bar, beta1, beta2, etaStart, vBar, muBar,
                       psiBar) {
  # potentially y - muFun can explode, e.g. if muFun = exp()
  uBar <- (psiBar^(-0.5)) * vBar * (y - muBar)

  as.vector(X1bar %*% beta1  + X2bar %*% beta2 + uBar)
}


#' Multiplies all rows of matrix X with vector y
#' @noRd
multAllRows <- function(X, y) {
  t(t(X) * y)
}


#' Internal function: Semiorthogonal-type transformation of X2 to Z2
#'
#' Uses the matrix Z2s (called \eqn{\bar{\Xi}} in eq. (9) of
#' \insertCite{deluca2018glm;textual}{WALS}) to transform \eqn{\bar{X}_2} to
#' \eqn{\bar{Z}_2}, i.e. to perform \eqn{\bar{Z}_2 = \bar{X}_2 \bar{\Delta}_2 \bar{\Xi}^{-1/2}}.
#' For WALS in the linear regression model, the variables do not have a "bar".
#'
#' @param Z2s Matrix for which we take negative square root in
#' \eqn{X2 * Delta2 * Z2s^{1/2}}.
#' @param X2 Design matrix of auxiliary regressors to be transformed to Z2
#' @param Delta2 Scaling matrix such that diagonal of
#' \eqn{\bar{\Delta}_2 \bar{X}_2^{\top} \bar{M}_1 \bar{X}_2 \Delta_{2}} is one
#' (ignored scaling by \eqn{n} because not needed in code).
#' See \insertCite{deluca2018glm;textual}{WALS}
#' @param SVD If \code{TRUE}, uses \code{\link[base]{svd}} to compute eigendecomposition
#' of \code{Z2s}, otherwise uses \code{\link[base]{eigen}}.
#' @param postmult If \code{TRUE}, then it uses
#' \eqn{Z2s^{-1/2} = T \Lambda^{-1/2} T^{\top}}, where \eqn{T} contains
#' the eigenvectors of \eqn{Z2s} in its columns and \eqn{\Lambda} the corresponding
#' eigenvalues. If \code{FALSE} it uses \eqn{Z2s^{-1/2} = T \Lambda^{-1/2}}.
#'
#' @section On the "semiorthogonal-type" transformation:
#' For WALS GLM (and WALS in the linear regression model),
#' the transformation is semiorthogonal (ignored scaling by \eqn{n} for clarity
#' and because it is not needed in the code)
#' in the sense that \eqn{\bar{M}_{1} \bar{Z}_{2}} is semiorthogonal since
#' \deqn{\bar{Z}_{2}^{\top} \bar{M}_1 \bar{Z}_{2} =
#' (\bar{Z}_{2}^{\top} \bar{M}_1) (\bar{M}_{1} \bar{Z}_{2}) = I_{k_2},}
#' where \eqn{\bar{M}_1} is an idempotent matrix.
#'
#' For WALS in the NB2 regression model, \eqn{\bar{M}_{1} \bar{Z}_{2}} is not
#' semiorthogonal anymore due to the rank-1 perturbation in \eqn{\bar{M}_1} which
#' causes \eqn{\bar{M}_1} to not be idempotent anymore, see
#' the section "Transformed model" in \insertCite{huynhwalsnb;textual}{WALS}.
#'
#'
#' @section On the use of \code{postmult = TRUE}:
#' The transformation of the auxiliary regressors \eqn{Z_2} for linear WALS in
#' eq. (12) of \insertCite{magnus2016wals;textual}{WALS} differs from the
#' transformation for WALS GLM (and WALS NB) in eq. (9) of
#' \insertCite{deluca2018glm;textual}{WALS}:
#'
#' In \insertCite{magnus2016wals;textual}{WALS} the transformed auxiliary
#' regressors are
#'
#' \deqn{Z_{2} = X_2 \Delta_2 T \Lambda^{-1/2},}
#'
#' where \eqn{T} contains the eigenvectors of
#' \eqn{\Xi = \Delta_2 X_{2}^{\top} M_{1} X_{2} \Delta_2} in the columns and
#' \eqn{\Lambda} the respective eigenvalues. This definition is used when
#' \code{postmult = FALSE}.
#'
#' In contrast, \insertCite{deluca2018glm;textual}{WALS} defines
#'
#' \deqn{Z_2 = X_2 \Delta_2 T \Lambda^{-1/2} T^{\top},}
#'
#' where we ignored scaling by \eqn{n} and the notation with "bar" for easier
#' comparison. This definition is used when \code{postmult = TRUE} and is
#' strongly preferred for \code{\link[WALS]{walsGLM}} and \code{\link[WALS]{walsNB}}.
#'
#' See \insertCite{huynhwals;textual}{WALS} for more details.
#'
#' @references
#' \insertAllCited{}
#'
semiorthogonalize <- function(Z2s, X2, Delta2, SVD = TRUE, postmult = FALSE) {
  if (SVD) { # use SVD to get eigenvectors and -values
    svdZ2s <- svd(Z2s)
    eigenVecs <- svdZ2s$u
    eigenVals <- svdZ2s$d
  } else {
    eigenZ2s <- eigen(Z2s)
    eigenVecs <- eigenZ2s$vectors
    eigenVals <- eigenZ2s$values
  }

  # check "numerical" rank of Z2s
  order <- max(dim(eigenVecs))
  rank <- sum(eigenVals > .Machine$double.eps)

  if (rank < order) stop("Z2s matrix is close to multicollinearity")

  # Computes "sort of" Delta2 * Xi^{-1/2} (Z2s in code) without postmultiplication
  # by eigenVecs. Xi^{-1/2} is not complete, instead computes
  # Delta2 * T * Lambda^{-1/2}, where Xi = T * Lambda * T' (EVD).
  # This is used in  Magnus & DeLuca, 2016, WALS Survey, p.126, eq. (12).
  D2 <- multAllRows(Delta2*eigenVecs, 1.0 / sqrt(eigenVals))

  # postmult uses proper Delta2 * Xi^{-1/2} computed by EVD, i.e.
  # Xi^{-1/2} = T * Lambda^{-1/2} T'.
  if (postmult) D2 <- D2 %*% t(eigenVecs)

  out <- list(Z2 = X2 %*% D2, D2 = D2, condition = max(eigenVals) / min(eigenVals))
  return(out)
}


#' Internal function: Transform gammas back to betas
#'
#' Transforms posterior means \eqn{\hat{\gamma}_2} and variances corresponding
#' to transformed auxiliary regressors \eqn{Z_2} back to regression coefficients
#' \eqn{\hat{\beta}} of original regressors \eqn{X_1} and \eqn{X_2}.
#'
#' @param posterior Object returned from \code{\link[WALS]{computePosterior}}.
#' @param y Response \eqn{y}.
#' @param Z1 Transformed focus regressors \eqn{Z_1}.
#' @param Z2 Transformed auxiliary regressors \eqn{Z_1}.
#' @param Delta1 \eqn{\Delta_1} or \eqn{\bar{\Delta}_1}.
#' @param D2 From \code{\link[WALS]{semiorthogonalize}}, if \code{postmult = FALSE}
#' was used, then D2 = \eqn{\Delta_2 T \Lambda^{-1/2}}, where \eqn{T} are the
#' eigenvectors of \eqn{\Xi} and \eqn{\Lambda} the diagonal matrix containing
#' the corresponding eigenvalues. If \code{postmult = TRUE} was used, then
#' D2 = \eqn{\Delta_2 T \Lambda^{-1/2} T^{\top} = \Delta_2 \Xi^{-1/2}}.
#' @param sigma Prespecified or estimated standard deviation of the error term.
#' @param Z1inv \eqn{(Z_{1}^{\top} Z_{1})^{-1}}.
#' @param method Character. \eqn{\hat{\gamma}_1} is obtained from a linear
#' regression of \eqn{Z_1} on pseudo-responses \eqn{y - Z_2 \hat{\gamma}_2}.
#' If \code{method = original}, then we use \code{\link[stats]{lm.fit}} to perform
#' the linear regression, if \code{method = "svd"}, then reuse the SVD of
#' \eqn{Z_1} in \code{svdZ1} to perform the regression.
#' @param svdZ1 Optional, only needed if \code{method = "svd"}. SVD of \eqn{Z_1}
#' computed using \code{\link[base]{svd}}.
#'
#' @details
#' The same transformations also work for GLMs, where we replace \eqn{X_1},
#' \eqn{X_2}, \eqn{Z_1} and \eqn{Z_2} with \eqn{\bar{X}_1}, \eqn{\bar{X}_2},
#' \eqn{\bar{Z}_1} and \eqn{\bar{Z}_2}, respectively. Generally, we need to
#' replace all variables with their corresponding "bar" version. Further,
#' for GLMs \code{sigma} is always 1.
#'
#' See \insertCite{magnus2016wals;textual}{WALS}, \insertCite{deluca2018glm;textual}{WALS}
#' and \insertCite{huynhwals;textual}{WALS} for the definitions of the variables.
#'
#' @references
#' \insertAllCited{}
#'
gammaToBeta <- function(posterior, y, Z1, Z2, Delta1, D2, sigma, Z1inv,
                        method = "original", svdZ1) {

  # Step 6: WALS estimates
  gamma2 <- sigma * posterior$postMean

  # c1 = (Z_1' Z_1)^(-1) Z_1' (y - Z_2 c_2) --> linear regression of
  # Z_1 on response (y - Z_2 c_2). P.136 of Magnus and De Luca WALS survey.
  pseudoY <- y - Z2 %*% gamma2

  if (method == "original") {
    # Can reuse QR-decomp of lm.fit from earlier in wals.default...
    # can recover (Z'Z)^-1 from the decomp instead of running
    # lm.fit again for another QR decomp and do not need to pass Z1inv
    outlm <- lm.fit(Z1, pseudoY, singular.ok = FALSE)
    gamma1 <- outlm$coefficients
  } else if (method == "svd") {
    # as.vector to stay consistent with output from lm.fit
    gamma1 <- as.vector(svdZ1$v %*% ((1.0 / svdZ1$d) * t(svdZ1$u)) %*% pseudoY)
  }

  beta1 <- Delta1*gamma1
  beta2 <- as.vector(D2 %*% gamma2)

  # Step 7: WALS precision
  Q <- Z1inv %*% crossprod(Z1, Z2)
  varGamma2 <- diag((sigma^2.0) * posterior$postVariance, nrow = ncol(D2),
                    ncol = ncol(D2))
  varBeta2 <- D2 %*% varGamma2 %*% t(D2)
  varGamma1 <- (sigma^2.0) * Z1inv + Q %*% varGamma2 %*% t(Q)
  varBeta1 <- multAllRows(Delta1*varGamma1, Delta1)

  covGamma1Gamma2 <- -Q %*% varBeta2
  covBeta1Beta2 <- multAllRows(covGamma1Gamma2, Delta1) %*% t(D2)

  vc <- rbind(cbind(varBeta1, covBeta1Beta2), cbind(t(covBeta1Beta2), varBeta2))
  vcGamma <- rbind(cbind(varGamma1, covGamma1Gamma2),
                   cbind(t(covGamma1Gamma2), varGamma2))

  return(list(coef = c(beta1, beta2), beta1 = beta1, beta2 = beta2,
              gamma1 = gamma1, gamma2 = gamma2,
              vcovBeta = vc, vcovGamma = vcGamma))
}

#' Internal function: Check singularity of SVDed matrix
#'
#' Checks whether matrix is singular based on singular values of SVD.
#'
#' @param singularValues Vector of singular values.
#' @param tol Absolute tolerance, singular if \code{min(singularValues) < tol}
#' @param rtol Relative tolerance, singular if
#'        \code{min(singularValues) / max(singularValues) < rtol}
#' @param digits The number significant digits to show in case a
#' warning is triggered by singularity.
#'
checkSingularitySVD <- function(singularValues, tol, rtol, digits = 5) {
  if (min(singularValues) < tol) {
    warning(paste("minimum singular value of",
                  signif(min(singularValues), digits),
                  ", is smaller than tol, possibly singular Z1"))
  }

  ratio <- min(singularValues) / max(singularValues)
  if (ratio < rtol) {
    warning(paste("Largest singular value is", signif(1/ratio, digits),
                  "times larger than smallest singular value, possibly singular Z1"))
  }
}


## offsets
expand_offset <- function(offset, n) {
  if (is.null(offset)) offset <- 0
  if (length(offset) == 1) offset <- rep.int(offset, n)
  as.vector(offset)
}

getOffset <- function(formula, mf, modelCall, n) {
  ## in mean part of formula
  offsetX <- expand_offset(model.offset(
    Formula::model.part(formula, data = mf, rhs = 1L, terms = TRUE)
    ), n)
  ## in precision part of formula
  offsetZ <- expand_offset(model.offset(
    Formula::model.part(formula, data = mf, rhs = 2L, terms = TRUE)
    ), n)
  ## in offset argument (used for mean)
  if (!is.null(modelCall$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"], n)
  ## collect
  return(list(mean = offsetX, precision = offsetZ))
}

# process weights
processWeights <- function(weights, mf, n) {
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- 1
  if (length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  return(weights)
}

# creates matrices for prediction from terms and contrasts and newdata
genNewdata <- function(terms, contrasts, newdata, na.action, xlev) {
  mf <- model.frame(delete.response(terms$full), newdata, na.action = na.action,
                    xlev = xlev$full)
  newdata <- newdata[rownames(mf), , drop = FALSE] # why do we need this?
  X1 <- model.matrix(delete.response(terms$focus), mf,
                     contrasts = contrasts$focus)
  X2 <- model.matrix(delete.response(terms$aux), mf,
                     contrasts = contrasts$aux)
  # y <- model.response(mf)
  ## HACK ##
  # drop constant from X2...
  # TODO: Can we solve this more elegantly?
  X2 <- X2[, -1L, drop = FALSE]

  return(list(X1 = X1, X2 = X2))
}

# get response from terms and newdata
getY <- function(terms, newdata, na.action = NULL) {
  mf <- model.frame(terms, newdata, na.action = na.action)
  return(model.response(mf))
}


extractModel <- function(formula, mf, data) {
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX1 <- terms(formula, data = data, rhs = 1L)
  mtX2 <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "any") # also allow factors for family = binomialWALS()
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
  X2 <- X2[,-1L, drop = FALSE] # intercept is always first column

  return(list(Y = Y, X1 = X1, X2 = X2, mt = mt, mtX1 = mtX1, mtX2 = mtX2,
              cont = cont))
}

#' Internal function: Computes X2M1X2 for walsNB when SVD is applied to Z1
#'
#' Exploits the SVD of \eqn{\bar{Z}_1} to compute
#' \eqn{\bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2}} to avoid directly inverting
#' \eqn{\bar{Z}_{1}^{\top} \bar{Z}_{1}}.
#'
#'
#' @param X2 Design matrix for auxiliary regressors
#' @param X2start Transformed design matrix for auxiliary regressors. Refers to
#' \eqn{\bar{X}_{2} = \bar{\Psi}^{1/2} X_{2}}.
#' @param qStart Vector \eqn{\bar{q}}, see section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param U \eqn{U} of SVD of \eqn{Z_1}. See details.
#' @param UellStart Vector \eqn{U \bar{\ell}}, see details.
#' @param ellStart Vector \eqn{\bar{\ell}} see details.
#' @param psiStart Diagonal matrix \eqn{\bar{\Psi}}, see section
#' "One-step ML estimator" of \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param gStart Derivative of dispersion parameter \eqn{\rho} of NB2 with
#' respect to \eqn{\alpha = \log(\rho)} evaluated at starting values of
#' one-step ML. \code{gStart} is a scalar.
#' See section "ML estimation" of  \insertCite{huynhwalsnb;textual}{WALS}.
#' @param epsilonStart Scalar \eqn{\bar{\epsilon}}, see section
#' "One-step ML estimator" of \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param geB \eqn{\bar{g} \bar{\epsilon} / (1 + B)}. In code
#' \code{gStart*epsilonStart / (1+B)}. See details for definition of \eqn{B}.
#' \code{gStart} is \eqn{\bar{g}} and \code{epsilonStart} is \eqn{\bar{\epsilon}}.
#'
#'
#' @details See section
#' "Simplification for computing \eqn{\bar{X}_{2}^{\top} \bar{M}_{1} \bar{X}_{2}}"
#' in the appendix of \insertCite{huynhwals;textual}{WALS} for details of the
#' implementation and for the definitions of arguments \code{Uellstart},
#' \code{ellStart}, and \code{geB}.
#'
#' All parameters that contain "start" feature the starting values for the
#' one-step ML estimation of submodels. See section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for details.
#'
#'
#' @references
#' \insertAllCited{}
#'
computeX2M1X2 <- function(X2, X2start, qStart, U, UellStart, ellStart, psiStart,
                          gStart, epsilonStart, geB) {
  X2tq <- colSums(X2 * qStart)
  X2plus <- (X2start +
    (gStart*epsilonStart) * tcrossprod((qStart/sqrt(psiStart)), X2tq))

  temp <- crossprod(U, X2plus)

  # see p.(22) of WALS-NegBin log rho derivation
  Z1invZ1X2plus <- (U %*% temp -
    geB * (UellStart %*% crossprod(ellStart, temp)))

  return(crossprod(X2start, X2start) +
    gStart*epsilonStart * tcrossprod(X2tq, X2tq) - crossprod(X2plus, Z1invZ1X2plus))
}


## computes gamma (coefs of transformed regressors Z) for walsNB
computeGammaUn <- function(U, V, singularVals,Z2, Z2start, Z2tq, Z1y0s, Z2y0s,
                           y0Start, UellStart, gStart, ellStart, geB, B,
                           qStart, psiStart, tStart, epsilonStart) {
  Z1startE <- (# see eq. 19.1.0 in WALS -NegBin variable log rho
    U %*% crossprod(U, y0Start) -
      geB * (UellStart %*% crossprod(UellStart, y0Start)) -
      (tStart * epsilonStart) * (U %*% crossprod(U, qStart / sqrt(psiStart)))
  )

  Z2plus <- Z2start + (gStart*epsilonStart) * tcrossprod(qStart/sqrt(psiStart), Z2tq)
  gammaUn2 <- Z2y0s - crossprod(Z2plus, Z1startE)


  temp <- sum(ellStart * crossprod(U, y0Start)) # ell'U'y_0
  Vpsiq <- V %*% (crossprod(U, qStart/sqrt(psiStart)) / singularVals)
  E <- (V %*% (crossprod(U, y0Start) / singularVals) -
          geB * (V %*%  (ellStart * temp / singularVals)) -
          (tStart*epsilonStart) * Vpsiq)

  # Compute D
  UtZ2start <- crossprod(U, Z2start)
  term1 <- (V %*% (UtZ2start/singularVals) -
              geB * (Vpsiq %*% crossprod(ellStart, UtZ2start)))
  term2 <- ((gStart * epsilonStart) * tcrossprod(Vpsiq, Z2tq) -
              (((gStart * epsilonStart)^2)/(1 + B)) *
              (V %*% tcrossprod(ellStart,sum(ellStart^2)*Z2tq))/singularVals)

  D <- term1 + term2

  # see eq. (20.1.1) in WALS-NegBin variable log rho on p.(20.1)
  gammaUn1 <- E + D %*% (crossprod(Z2plus, Z1startE) - Z2y0s)

  return(list(gammaUn1 = gammaUn1, gammaUn2 = gammaUn2, D = D))
}

#' Internal function: Computes unrestricted one-step ML estimator for transformed
#' regressors in walsNB
#'
#' Computes one-step ML estimator for the unrestricted model in walsNB
#' (coefs of transformed regressors \eqn{\bar{Z}})
#' by using SVD on entire transformed design matrix \eqn{\bar{Z}}.
#' The matrix \eqn{\bar{Z}} should have full column rank.
#'
#' @param U Left singular vectors of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param V Right singular vectors of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param singularVals Singular values of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param ellStart Vector \eqn{\bar{\ell}} see details.
#' @param gStart Derivative of dispersion parameter \eqn{\rho} of NB2 with
#' respect to \eqn{\alpha = \log(\rho)} evaluated at starting values of
#' one-step ML. \code{gStart} is a scalar.
#' See section "ML estimation" of  \insertCite{huynhwalsnb;textual}{WALS}.
#' @param epsilonStart Scalar \eqn{\bar{\epsilon}}, see section
#' "One-step ML estimator" of \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param qStart Vector \eqn{\bar{q}}, see section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param y0Start Vector \eqn{\bar{y}_0}, see section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param tStart Scalar \eqn{\bar{t}}, see section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for definition.
#' @param psiStart Diagonal matrix \eqn{\bar{\Psi}}, see section
#' "One-step ML estimator" of \insertCite{huynhwalsnb;textual}{WALS} for definition.
#'
#' @details
#' See section "Simplification for computing \eqn{\tilde{\gamma}_{u}}"
#' in the appendix of \insertCite{huynhwals;textual}{WALS} for details of the
#' implementation and for the definitions of argument \code{ellStart}.
#'
#' All parameters that contain "start" feature the starting values for the
#' one-step ML estimation of submodels. See section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for details.
#'
#' Uses \code{\link[WALS]{svdLSplus}} under-the-hood.
#'
#' @references
#' \insertAllCited{}
#'
#'
computeGammaUnSVD <- function(U, V, singularVals, ellStart, gStart, epsilonStart, qStart,
                              y0Start, tStart, psiStart) {

  B <- gStart * epsilonStart * sum(ellStart^2.0)
  geB <- (gStart * epsilonStart / (1 + B))
  yDiff <- y0Start - tStart * epsilonStart * (qStart / sqrt(psiStart))

  svdLSplus(U, V, singularVals, yDiff, ellStart, geB)
}


#' Internal function: Compute model-averaged estimator of focus regressors in walsNB
#'
#' Exploits the SVD of the design matrix of the focus regressors \eqn{\bar{Z}_1},
#' the model-averaged estimator for the auxiliary regressors
#' \eqn{\hat{\gamma}_{2}} and the Sherman-Morrison-Woodbury
#' formula for computing the model-averaged estimator of the focus regressors
#' in walsNB.
#'
#' @param gamma2 Model-averaged estimate for auxiliary regressors
#' from \code{\link[WALS]{computePosterior}}.
#' @param Z2start Transformed design matrix of auxiliary regressors \eqn{\bar{Z}_2}.
#' See details.
#' @param Z2 Another transformed design matrix of auxiliary regressors \eqn{Z_2}.
#' See details.
#' @param U Left singular vectors of \eqn{\bar{Z}_1} from \code{\link[base]{svd}}.
#' @param V Right singular vectors of \eqn{\bar{Z}_1} from \code{\link[base]{svd}}.
#' @param singularVals Singular values of \eqn{\bar{Z}_1} from \code{\link[base]{svd}}.
#' @inheritParams computeGammaUnSVD
#'
#' @details
#' See section "Simplification for computing \eqn{\hat{\gamma}_{1}}"
#' in the appendix of \insertCite{huynhwals;textual}{WALS} for details of the
#' implementation and for the definitions of argument \code{ellStart}.
#'
#' All parameters that contain "start" feature the starting values for the
#' one-step ML estimation of submodels. See section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for details.
#'
#' The argument \code{Z2start} is defined as \insertCite{huynhwalsnb}{WALS}
#'
#' \deqn{
#' \bar{Z}_{2} := \bar{X}_{2} \bar{\Delta}_{2} \bar{\Xi}^{-1/2},
#' }
#'
#' and \code{Z2} is defined as
#'
#' \deqn{
#' Z_{2} := X_{2} \bar{\Delta}_{2} \bar{\Xi}^{-1/2}.
#' }
#'
#' Uses \code{\link[WALS]{svdLSplus}} under-the-hood.
#'
#' @references
#' \insertAllCited{}
#'
computeGamma1 <- function(gamma2, Z2start, Z2, U, V, singularVals, ellStart, gStart,
                          epsilonStart, qStart, y0Start, tStart, psiStart) {

  B <- gStart * epsilonStart * sum(ellStart^2.0)
  geB <- (gStart * epsilonStart / (1 + B))

  qScaled <- qStart/sqrt(psiStart)
  yDiff <- y0Start - tStart * epsilonStart * qScaled
  Z2toY <- (Z2start %*% gamma2 + (gStart * epsilonStart *
               (qScaled %*% (crossprod(qStart, Z2) %*% gamma2)) ))
  pseudoY <- yDiff - Z2toY

  svdLSplus(U, V, singularVals, pseudoY, ellStart, geB)

}


#' Internal function: Computes fully restricted one-step ML estimator for
#' transformed regressors in walsNB
#'
#' Computes one-step ML estimator of fully restricted model
#' (coefs of transformed regressors of \eqn{\bar{Z}_1})
#' in walsNB by using SVD on transformed design matrix of the focus regressors
#' \eqn{\bar{Z}_1}. The matrix \eqn{\bar{Z_1}} should have full column rank.
#'
#' @inheritParams computeGamma1
#'
#' @details
#' See section "Simplification for computing \eqn{\tilde{\gamma}_{1r}}"
#' in the appendix of \insertCite{huynhwals;textual}{WALS} for details of the
#' implementation and for the definitions of argument \code{ellStart}.
#'
#' All parameters that contain "start" feature the starting values for the
#' one-step ML estimation of submodels. See section "One-step ML estimator" of
#' \insertCite{huynhwalsnb;textual}{WALS} for details.
#'
#' Uses \code{\link[WALS]{svdLSplus}} under-the-hood.
#'
#' @references
#' \insertAllCited{}
#'
computeGamma1r <- function(U, V, singularVals, ellStart, gStart,
                           epsilonStart, qStart, y0Start, tStart, psiStart) {
  ## TODO: computeGamma1R does exactly the same thing as computeGammaUnSVD
  ## --> do not need this function!
  ## Compute only fully restricted gamma using SVD of Z1start
  B <- gStart * epsilonStart * sum(ellStart^2.0)
  geB <- (gStart * epsilonStart / (1 + B))

  qScaled <- qStart/sqrt(psiStart)
  yDiff <- y0Start - tStart * epsilonStart * qScaled
  svdLSplus(U, V, singularVals, yDiff, ellStart, geB)

}

#' Internal function: Uses SVD components to compute final estimate via
#' Sherman-Morrison-Woodbury formula.
#'
#' Solves the equation system in walsNB via Sherman-Morrison-Woodbury formula
#' for the unrestricted estimator \eqn{\hat{\gamma}_{u}}.
#'
#' @param U Left singular vectors of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param V Right singular vectors of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param singularVals Singular values of \eqn{\bar{Z}} or \eqn{\bar{Z}_{1}}
#' from \code{\link[base]{svd}}.
#' @param y "Pseudo"-response, see details.
#' @param ell Vector \eqn{\bar{\ell}} from section
#' "Simplification for computing \eqn{\tilde{\gamma}_{u}}"
#' \insertCite{huynhwals;textual}{WALS}
#' @param geB Scalar \eqn{\bar{g} \bar{\epsilon} / (1 + B)}. See section
#' "Simplification for computing \eqn{\tilde{\gamma}_{u}}"
#' \insertCite{huynhwals;textual}{WALS} for definition of
#' \eqn{\bar{g}}, \eqn{\bar{\epsilon}} and \eqn{B}.
#'
#'
#' @details
#' The function can be reused for the computation of the fully restricted
#' estimator \eqn{\tilde{\gamma}_{1r}} and the model averaged estimator
#' \eqn{\hat{\gamma}_{1}}.
#'
#' For \eqn{\tilde{\gamma}_{1r}} and \eqn{\hat{\gamma}_{1}} use
#' \code{U}, \code{V} and \code{singularVals} from SVD of \eqn{\bar{Z}_{1}}.
#'
#' For \eqn{\hat{\gamma}_{u}} and \eqn{\tilde{\gamma}_{1r}} use same
#' pseudo-response \eqn{\bar{y_{0}} - \bar{t} \bar{\epsilon} \bar{\Psi}^{-1/2} \bar{q}}
#' in argument \code{y}.
#'
#' For \eqn{\hat{\gamma}_{1}} use pseudo-response
#' \eqn{\bar{y_{0}} - \bar{t} \bar{\epsilon} \bar{\Psi}^{-1/2} \bar{q} -
#' (\bar{Z}_{2} + \bar{g} \bar{\epsilon} \bar{\Psi}^{-1/2} \bar{q} \bar{q}^{\top}
#' Z_{2}) \hat{\gamma}_{2}}.
#'
#' See section "Note on function svdLSplus from WALS"
#' in \insertCite{huynhwals;textual}{WALS}.
#'
#' @references
#' \insertAllCited{}
#'
svdLSplus <- function(U, V, singularVals, y, ell, geB) {
  A <- V %*% (crossprod(U, y)/singularVals)
  ellUy <- sum(ell * crossprod(U, y))
  E <- geB * ellUy * (V %*% (ell/singularVals))

  return(A - E)
}

## Helpers for summary methods
printCallCoefs <- function(x, digits, ...) {
  cat(paste0("\nCoefficients of k1 = ", x$k1, " focus regressors: \n"))
  printCoefmat(x$focusCoefs, digits = digits,...)

  cat(paste0("\nCoefficients of k2 = ", x$k2, " auxiliary regressors: \n"))
  printCoefmat(x$auxCoefs, digits = digits,...)
}

printPriorNKappa <- function(x, digits) {
  priorPars <- paste(names(x$prior$printPars), signif(x$prior$printPars, digits),
                     sep = " = ", collapse = ", ")
  cat(paste0("\nPrior: ", x$prior$prior, "(", priorPars, ")"))

  cat(paste0("\nNumber of observations: ", x$n))

  cat(paste0("\nKappa: ", signif(sqrt(x$condition), 3), "\n"))
}
