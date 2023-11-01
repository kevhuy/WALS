test_that("wals.formula with only constant focus works", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)
  fWals <- math ~ 1 | read + stratio + english + lunch + expenditure

  # expect no error
  expect_error(wals(fWals, data = dd, prior = weibull(), method = "svd",
                    eigenSVD = TRUE, keepUn = TRUE),
               regexp = NA)
})

test_that("some class methods of wals", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)
  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  # expect no error in familyPrior extraction.
  expect_error(familyPrior(walsEst), regexp = NA)

  expect_length(coef(walsEst), 6L) # check coefs
  expect_true(all(dim(vcov(walsEst)) == c(6L,6L))) # check vcov
  expect_true(nobs(walsEst) == nrow(dd)) # check nobs
})

test_that("Estimates match Magnus et al. (2010), Journal of Econometrics", {
  tol <- 1e-10
  # original values from Table 2 in Magnus et al. (2010)
  coefVals <- c("(Intercept)" = .0594,
                "lgdp60" = -.0156,
                "equipinv" = .1555,
                "school60" = .0175,
                "life60" = .0009,
                "popgrowth" = .2651,
                "law" = .0147,
                "tropics" = -.0055,
                "avelf" = -.0053,
                "confucian" = .0443)
  seVals <- c("(Intercept)" = 0.0221,
              "lgdp60" = .0033,
              "equipinv" = .0551,
              "school60" = .0097,
              "life60" = .0004,
              "popgrowth" = .2487,
              "law" = .0065,
              "tropics" = .0037,
              "avelf" = .0048,
              "confucian" = .0163)

  # Important: prescale = FALSE, still used old version of WALS in Magnus et al. (2010)
  fitWals <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
                    law + tropics + avelf + confucian, data = GrowthMPP,
                  prior = laplace(), prescale = FALSE)

  expect_equal(round(coef(fitWals), 4), coefVals, tolerance = tol)
  expect_equal(round(sqrt(diag(vcov(fitWals))), 4), seVals, tolerance = tol)
})

test_that("Estimates match De Luca & Magnus (2011), The Stata Journal", {
  tol <- 1e-10
  # original values from table on p. 534 of De Luca & Magnus (2011)
  coefVals <- c("(Intercept)" = .0617514,
                "lgdp60" = -.0156501,
                "equipinv" = .1582128,
                "school60" = .0166758,
                "life60" = .0008515,
                "popgrowth" = .2713869,
                "law" = .0134105,
                "tropics" = -.0059973,
                "avelf" = -.0076757,
                "confucian" = .046455)
  seVals <- c("(Intercept)" = .0217909,
              "lgdp60" = .0031439,
              "equipinv" = .054421,
              "school60" = .009667,
              "life60" = .0003505,
              "popgrowth" = .2425285,
              "law" = .0058037,
              "tropics" = .0034556,
              "avelf" = .0050657,
              "confucian" = .0142765)

  # Difference to Magnus et al. (2010): prescale = TRUE, use recent version of
  # WALS which is scale independent thanks to prescaling of regressors.
  fitWals <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
                    law + tropics + avelf + confucian, data = GrowthMPP,
                  prior = laplace(), prescale = TRUE)

  expect_equal(round(coef(fitWals), 7), coefVals, tolerance = tol)
  expect_equal(round(sqrt(diag(vcov(fitWals))), 7), seVals, tolerance = tol)
})

test_that("walsMatrix coefs and covmat equal to wals", {
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)
  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  walsEstMatrix <- wals(walsEst$x$focus, x2 = walsEst$x$aux, y = walsEst$y,
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

  walsEstMatrix <- wals(walsEst$x$focus, x2 = walsEst$x$aux, y = walsEst$y,
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

  walsEstMatrix <- wals(walsEst$x$focus, x2 = walsEst$x$aux, y = walsEst$y,
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

  walsMatrix <- wals(X1, x2 = X2, y = y, method = "original",
                     prior = weibull())
  walsDefault <- wals.default(X1, x2 = X2, y = y, method = "original",
                              prior = weibull())

  # check coefs & covariance matrix
  expect_identical(coef(walsMatrix), coef(walsDefault))
  expect_identical(vcov(walsMatrix), vcov(walsDefault))
})

test_that("Different methods for wals yield same results", {
  ## Check if estimated regression coefficients from different methods
  ## yield same results.
  tol <- 1e-06 # relative tolerance for deviations
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
