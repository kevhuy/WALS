test_that("walsNB.formula with only constant focus works", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ 1 | read + stratio + english + lunch + expenditure

  # expect no error
  expect_error(wals(fWals, data = dd, prior = weibull(), method = "svd",
                    eigenSVD = TRUE),
               regexp = NA)
})

test_that("some class methods of wals", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  # check coefs
  expect_length(coef(walsEst), 6L)

  # check vcov
  expect_true(all(dim(vcov(walsEst)) == c(6L,6L)))

  # nobs
  expect_true(nobs(walsEst) == nrow(dd))
})

test_that("walsMatrix coefs and covmat equal to wals", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  walsEstMatrix <- wals(walsEst$x$focus, X2 = walsEst$x$aux, y = walsEst$y,
                        method = "original", eigenSVD = TRUE, prior = weibull())

  # check coefs & covariance matrix
  expect_equal(coef(walsEst), coef(walsEstMatrix))
  expect_equal(vcov(walsEst), vcov(walsEstMatrix))
})

test_that("walsMatrix predictions equal to wals", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  # add artificial factors for testing
  dd$englishFactor <- as.factor(dd$english > 20)
  dd$incomeFactor <- as.factor(dd$income > 17)

  fWals <- math ~ read + stratio | englishFactor + lunch + expenditure + incomeFactor

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  walsEstMatrix <- wals(walsEst$x$focus, X2 = walsEst$x$aux, y = walsEst$y,
                        prior = weibull(), method = "original", keepX = TRUE)

  # check predictions
  pred1 <- predict(walsEst, newdata = dd)
  pred2 <- predict(walsEstMatrix, newX1 = walsEstMatrix$x$focus,
                   newX2 = walsEstMatrix$x$aux)
  expect_equal(pred1, pred2)
})

test_that("Predictions are correct", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  # add artificial factors for testing
  dd$englishFactor <- as.factor(dd$english > 20)
  dd$incomeFactor <- as.factor(dd$income > 17)

  fWals <- math ~ read + stratio | englishFactor + lunch + expenditure + incomeFactor

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)
  meanFormula <- as.vector(model.matrix(walsEst, "focus") %*% coef(walsEst, "focus")
                           + model.matrix(walsEst, "aux") %*% coef(walsEst, "aux"))

  walsEstMatrix <- wals(walsEst$x$focus, X2 = walsEst$x$aux, y = walsEst$y,
                        prior = weibull(), method = "original", keepX = TRUE)
  meanMatrix <- as.vector(walsEst$x$focus %*% coef(walsEstMatrix, "focus")
                          + walsEst$x$aux %*% coef(walsEstMatrix, "aux"))

  # check predictions
  predFormula <- as.vector(predict(walsEst, newdata = dd))
  predMatrix <- as.vector(predict(walsEstMatrix, newX1 = walsEstMatrix$x$focus,
                                  newX2 = walsEstMatrix$x$aux))

  expect_identical(predFormula, meanFormula)
  expect_identical(predMatrix, meanMatrix)
})


test_that("wals.matrix and wals.default identical", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  y <- dd$math
  X1 <- as.matrix(cbind(1, CASchools[, c("read", "stratio")]))
  X2 <- as.matrix(CASchools[, c("english", "lunch", "expenditure")])

  walsMatrix <- wals(X1, X2 = X2, y = y, method = "original",
                     prior = weibull())
  walsDefault <- wals.default(X1, X2 = X2, y = y, method = "original",
                              prior = weibull())

  # check coefs & covariance matrix
  expect_identical(coef(walsMatrix), coef(walsDefault))
  expect_identical(vcov(walsMatrix), vcov(walsDefault))
})


test_that("Different methods for wals yield same results", {
  ## Check if estimated regression coefficients from different methods
  ## yield same results.

  tol <- 1e-06 # relative tolerance for deviations

  ## Test on CASchools dataset
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsfSVD <- wals(fWals, data = dd, prior = weibull(), method = "svd",
                   eigenSVD = TRUE)

  walsOriginal <- wals(fWals, data = dd, prior = weibull(), method = "original",
                       eigenSVD = FALSE)

  expect_equal(coef(walsfSVD), coef(walsOriginal), tolerance = tol)
})

test_that("Fitted values and prediction on same dataset are identical", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, prior = laplace(), method = "svd",
                  eigenSVD = TRUE)

  expect_identical(predict(walsEst, newdata = dd), fitted(walsEst))
})
