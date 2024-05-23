test_that("walsNB estimation converges", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)

  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6,
                   verbose = TRUE)

  expect_true(nbWals$converged)
})

test_that("walsNB limiting nIt works and returns finite coefficients", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  # maximum value of L2 norm of coefficients
  # For NMES1988 dataset and given regressors, it should not exceed this value!
  limit <- 10

  # Check if limiting number of iteration works
  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6, nIt = 50,
                   verbose = TRUE)
  coefs <- coef(nbWals)
  coefsNorm <- norm(coefs, type = "2")

  # expect that L2 norm of coefficients is smaller than limit and finite
  expect_true(is.finite(coefsNorm))
  expect_true(coefsNorm < limit)

})

test_that("Different walsNB class methods yield same result", {
  tol <- 1e-8
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6,
                   verbose = TRUE, keepY = TRUE, keepX = TRUE)

  nbWalsMatrix <- walsNB(nbWals$x$focus, x2 = nbWals$x$aux, y = nbWals$y,
                         link = "log", prior = weibull(), method = "fullSVD",
                         iterate = TRUE, tol = 1e-6,
                         verbose = TRUE)

  nbWalsDefault <- walsNB.default(nbWals$x$focus, x2 = nbWals$x$aux,
                                  y = nbWals$y, link = "log", prior = weibull(),
                                  method = "fullSVD", iterate = TRUE, tol = 1e-6,
                                  verbose = TRUE)

  # check coefs
  expect_equal(coef(nbWals), coef(nbWalsMatrix), tolerance = tol)
  expect_equal(coef(nbWals), coef(nbWalsDefault), tolerance = tol)
})

test_that("walsNBmatrix predictions equal walsNB predictions", {
  tol <- 1e-8
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6,
                   verbose = TRUE, keepY = TRUE, keepX = TRUE)

  nbWalsMatrix <- walsNB(nbWals$x$focus, x2 = nbWals$x$aux, y = nbWals$y,
                         link = "log", prior = weibull(), method = "fullSVD",
                         iterate = TRUE, tol = 1e-6,
                         verbose = TRUE)

  # check predictions
  pred1 <- predict(nbWals, newdata = dd, type = "density")
  pred2 <- predict(nbWalsMatrix, newX1 = nbWals$x$focus, newX2 = nbWals$x$aux,
                   newY = nbWals$y, type = "density")

  expect_equal(pred1, pred2, tolerance = tol)
})

test_that("residuals of walsNB return correct values", {
  tol <- 1e-6
  data("NMES1988", package = "AER")
  NMES1988 <- na.omit(NMES1988)

  fitNB <- walsNB(visits ~ health + chronic + age + gender |
                    I((age^2)/10) + married + region,
                  data = NMES1988, prior = laplace(), verbose = TRUE)

  y <- NMES1988$visits
  mu <- predict(fitNB, newdata = NMES1988)
  rho <- fitNB$rho

  pearson <- (y - mu) / sqrt(mu + (mu^2)/rho)
  devres <- pmax(2 * (y*log(pmax(1, y)/mu) - (y + rho)*log((y + rho) / (mu + rho))), 0)
  devres <- sqrt(devres)
  devres <- ifelse(y - mu > 0, devres, -devres)
  resids <- y - mu

  expect_equal(pearson, residuals(fitNB, type = "pearson"), tolerance = tol)
  expect_equal(devres, residuals(fitNB, type = "deviance"), tolerance = tol)
  expect_equal(resids, residuals(fitNB, type = "response"), tolerance = tol)
})

test_that("Predictions use correct link values", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  fitFormula <- walsNB(fWals, data = dd, link = "log", prior = weibull(),
                       method = "fullSVD", keepY = TRUE, keepX = TRUE)
  linkFormula <- as.vector(model.matrix(fitFormula, "focus") %*% coef(fitFormula, "focus")
                           + model.matrix(fitFormula, "aux") %*% coef(fitFormula, "aux"))

  fitMatrix <- walsNB(fitFormula$x$focus, fitFormula$x$aux, fitFormula$y,
                      prior = weibull())
  linkMatrix <- as.vector(fitFormula$x$focus %*% coef(fitMatrix, "focus")
                          + fitFormula$x$aux %*% coef(fitMatrix, "aux"))

  expect_identical(linkFormula, as.vector(predict(fitFormula, newdata = dd,
                                                  type = "link")))
  expect_identical(linkMatrix, as.vector(predict(fitMatrix, newX1 = fitFormula$x$focus,
                                                 newX2 = fitFormula$x$aux, type = "link")))
})

test_that("Probability prediction works", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd, link = "log",prior = weibull(),
                   method = "fullSVD")

  # expect no error in probability prediction.
  expect_error(predict(nbWals, newdata = dd, type = "prob"), regexp = NA)
})

test_that("Initialization with MASS::glm.nb runs", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  # expect no error in estimation.
  expect_error(walsNB(fWals, data = dd, link = "log",prior = weibull(),
                      method = "fullSVD", iterate = TRUE, tol = 1e-6,
                      verbose = TRUE, keepY = TRUE, keepX = TRUE,
                      controlInitNB = controlNB(initMASS = TRUE)),
               regexp = NA)

})

test_that("Initialization with restricted estimator", {
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  # expect no error in estimation.
  expect_error(walsNB(fWals, data = dd, link = "log",prior = weibull(),
                      method = "fullSVD", iterate = TRUE, tol = 1e-6,
                      verbose = TRUE, keepY = TRUE, keepX = TRUE,
                      controlInitNB = controlNB(initMASS = TRUE,
                                                restricted = TRUE)
                      ),
               regexp = NA)

  expect_error(walsNB(fWals, data = dd, link = "log",prior = weibull(),
                      method = "fullSVD", iterate = TRUE, tol = 1e-6,
                      verbose = TRUE, keepY = TRUE, keepX = TRUE,
                      controlInitNB = controlNB(initMASS = FALSE,
                                                restricted = TRUE)
                      ),
               regexp = NA)

  ## only constant as focus
  fWals <- (visits ~ 1 | health + chronic + age + I((age^2)/10) + insurance + medicaid +
              adl + region + gender + married + income + school + afam + employed)

  # expect no error in estimation.
  expect_error(walsNB(fWals, data = dd, link = "log",prior = weibull(),
                      method = "fullSVD", iterate = TRUE, tol = 1e-6,
                      verbose = TRUE, keepY = TRUE, keepX = TRUE,
                      controlInitNB = controlNB(initMASS = TRUE,
                                                restricted = TRUE)
  ),
  regexp = NA)

  expect_error(walsNB(fWals, data = dd, link = "log",prior = weibull(),
                      method = "fullSVD", iterate = TRUE, tol = 1e-6,
                      verbose = TRUE, keepY = TRUE, keepX = TRUE,
                      controlInitNB = controlNB(initMASS = FALSE,
                                                restricted = TRUE)
  ),
  regexp = NA)

})

test_that("Different methods for walsNB yield same results", {
  ## Check if estimated regression coefficients from different methods
  ## yield same results.
  tol <- 1e-06 # relative tolerance for deviations
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  # Turn off iteration in all models as iterating WALS algo will magnify
  # small numerical differences between the methods.

  nbWalsfSVD <- walsNB(fWals, data = dd,
                       link = "log",
                       prior = weibull(), method = "fullSVD",
                       iterate = FALSE, tol = 1e-6,
                       verbose = TRUE, keepY = TRUE, keepX = TRUE)

  nbWalsOriginal <- walsNB(fWals, data = dd,
                          link = "log",
                          prior = weibull(), method = "original",
                          iterate = FALSE, tol = 1e-6,
                          verbose = TRUE, keepY = TRUE, keepX = TRUE)

  expect_equal(coef(nbWalsfSVD), coef(nbWalsOriginal), tolerance = tol)
})

test_that("One-part formula uses all variables as auxiliary regressors", {
  ## Check if y ~ x1 + x2 is the same as y ~ 1 | x1 + x2
  tol <- 1e-08 # relative tolerance for deviations
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)
  fOne <- visits ~ health + chronic + age + insurance + medicaid
  fTwo <- visits ~ 1 | health + chronic + age + insurance + medicaid

  nbOne <- walsNB(fOne, data = dd, link = "log", prior = laplace(),
                  method = "fullSVD")
  nbTwo <- walsNB(fTwo, data = dd, link = "log", prior = laplace(),
                  method = "fullSVD")

  expect_equal(coef(nbOne), coef(nbTwo), tolerance = tol)
})
