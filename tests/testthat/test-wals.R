test_that("walsNB.formula with only constant focus works", {
  ## Test on CASchools dataset
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ 1 | read + stratio + english + lunch + expenditure

  # expect no error
  expect_error(wals(fWals, data = dd, prior = weibull(), method = "svd",
                    eigenSVD = TRUE),
               regexp = NA)
})


test_that("walsNB.matrix coefs and covmat equal to walsNB.formula", {
  ## Test on CASchools dataset
  data("CASchools", package = "AER")
  CASchools$stratio <- CASchools$students / CASchools$teachers
  dd <- na.omit(CASchools)

  fWals <- math ~ read + stratio | english + lunch + expenditure

  walsEst <- wals(fWals, data = dd, method = "original", eigenSVD = TRUE,
                  prior = weibull(), keepY = TRUE, keepX = TRUE)

  walsEstMatrix <- wals(walsEst$x$focus, X2 = walsEst$x$aux, y = walsEst$y,
                         prior = weibull(), method = "original")

  # check coefs & covariance matrix
  expect_equal(coef(walsEst), coef(walsEstMatrix))
  expect_equal(walsEst$vcovBeta, walsEstMatrix$vcovBeta)
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

