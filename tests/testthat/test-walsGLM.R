test_that("walsGLM estimation converges", {
  ## Test Binomial on HMDA dataset
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, data = HMDA, family = binomialWALS(),
                         prior = weibull(), iterate = TRUE)
  expect_true(fitBinomial$converged)

  ## Test Poisson on NMES1988 dataset
  data("NMES1988", package = "AER")
  NMES1988 <- na.omit(NMES1988)

  fitPoisson <- walsGLM(emergency ~ health + chronic + age + gender |
                          I((age^2)/10) + married + region, data = NMES1988,
                        family = poissonWALS(), prior = laplace(), iterate = TRUE)
  expect_true(fitPoisson$converged)
})


test_that("walsGLM limiting nIt works and returns finite coefficients", {
  ## Test on HMDA dataset
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  # maximum value of L2 norm of coefficients
  # For HMDA dataset and given regressors, it should not exceed this value!
  limit <- 10

  fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, data = HMDA, family = binomialWALS(),
                         prior = weibull(), iterate = TRUE)

  expect_true(fitBinomial$converged)
  coefs <- coef(fitBinomial)
  coefsNorm <- norm(coefs, type = "2")

  # expect that L2 norm of coefficients is smaller than limit and finite
  expect_true(is.finite(coefsNorm))
  expect_true(coefsNorm < limit)
})


test_that("Different methods for walsNB yield same results", {
  ## Check if estimated regression coefficients from different methods
  ## yield same results.

  tol <- 1e-06 # relative tolerance for deviations

  ## Test on HMDA dataset
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  fWals <- deny ~ pirat + hirat + lvrat + chist + mhist + phist | selfemp + afam

  # Turn off iteration in all models as iterating WALS algo will magnify
  # small numerical differences between the methods.

  glmWALSsvd <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, data = HMDA, family = binomialWALS(),
                         prior = weibull(), iterate = TRUE, method = "svd")

  glmWALSoriginal <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                              selfemp + afam, data = HMDA, family = binomialWALS(),
                            prior = weibull(), iterate = TRUE, method = "original")

  expect_equal(coef(glmWALSsvd), coef(glmWALSoriginal), tolerance = tol)
})

