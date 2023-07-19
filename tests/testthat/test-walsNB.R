test_that("walsNB estimation converges", {
  ## Test on NMES1988 dataset
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
  ## Test on NMES1988 dataset
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

test_that("walsNBmatrix estimation returns same coefs as with formula", {
  ## Test on NMES1988 dataset
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)

  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6,
                   verbose = TRUE, keepY = TRUE, keepX = TRUE)

  nbWalsMatrix <- walsNB(nbWals$x$focus, X2 = nbWals$x$aux, y = nbWals$y,
                         link = "log", prior = weibull(), method = "fullSVD",
                         iterate = TRUE, tol = 1e-6,
                         verbose = TRUE)

  # check coefs
  expect_equal(coef(nbWals), coef(nbWalsMatrix))
})

test_that("walsNBmatrix predictions equal walsNB predictions", {
  ## Test on NMES1988 dataset
  data("NMES1988", package = "AER")
  dd <- na.omit(NMES1988)

  fWals <- (visits ~ health + chronic + age + I((age^2)/10) + insurance + medicaid |
              adl + region + gender + married + income + school + afam + employed)

  nbWals <- walsNB(fWals, data = dd,
                   link = "log",
                   prior = weibull(), method = "fullSVD",
                   iterate = TRUE, tol = 1e-6,
                   verbose = TRUE, keepY = TRUE, keepX = TRUE)

  nbWalsMatrix <- walsNB(nbWals$x$focus, X2 = nbWals$x$aux, y = nbWals$y,
                         link = "log", prior = weibull(), method = "fullSVD",
                         iterate = TRUE, tol = 1e-6,
                         verbose = TRUE)

  # check predictions
  pred1 <- predict(nbWals, newdata = dd, type = "density")
  pred2 <- predict(nbWalsMatrix, newX1 = nbWals$x$focus, newX2 = nbWals$x$aux,
                   newY = nbWals$y, type = "density")

  expect_equal(pred1, pred2)
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

  ## Test on NMES1988 dataset
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
