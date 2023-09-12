test_that("walsGLM estimation converges", {
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)
  fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, family = binomialWALS(), data = HMDA,
                         prior = weibull(), iterate = TRUE)
  expect_true(fitBinomial$converged)

  data("NMES1988", package = "AER")
  NMES1988 <- na.omit(NMES1988)
  fitPoisson <- walsGLM(emergency ~ health + chronic + age + gender |
                          I((age^2)/10) + married + region, family = poissonWALS(),
                        data = NMES1988, prior = laplace(), iterate = TRUE)
  expect_true(fitPoisson$converged)
})

test_that("walsGLM limiting nIt works and returns finite coefficients", {
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  # maximum value of L2 norm of coefficients
  # For HMDA dataset and given regressors, it should not exceed this value!
  limit <- 10

  fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, family = binomialWALS(), data = HMDA,
                         prior = weibull(), iterate = TRUE)

  expect_true(fitBinomial$converged)
  coefs <- coef(fitBinomial)
  coefsNorm <- norm(coefs, type = "2")

  # expect that L2 norm of coefficients is smaller than limit and finite
  expect_true(is.finite(coefsNorm))
  expect_true(coefsNorm < limit)
})

test_that("Different class methods of walsGLM yield same results", {
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  fitFormula <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, data = HMDA, family = binomialWALS(),
                         prior = weibull(), keepX = TRUE, keepY = TRUE)
  fitMatrix <- walsGLM(fitFormula$x$focus, fitFormula$x$aux, fitFormula$y,
                       family = binomialWALS(), prior = weibull())
  fitDefault <- walsGLM.default(fitFormula$x$focus, fitFormula$x$aux, fitFormula$y,
                                family = binomialWALS(), prior = weibull())

  # expect no error in familyWALS extraction.
  expect_error(familyWALS(fitFormula), regexp = NA)

  expect_identical(coef(fitFormula), coef(fitMatrix), coef(fitDefault))
  expect_identical(vcov(fitFormula), vcov(fitMatrix), vcov(fitDefault))
})

test_that("Different methods for walsGLM yield same results", {
  ## Check if estimated regression coefficients from different methods
  ## yield same results.
  tol <- 1e-06 # relative tolerance for deviations
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)
  fWals <- deny ~ pirat + hirat + lvrat + chist + mhist + phist | selfemp + afam

  # Turn off iteration in all models as iterating WALS algo will magnify
  # small numerical differences between the methods.
  glmWALSsvd <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, family = binomialWALS(), data = HMDA,
                         prior = weibull(), iterate = TRUE, method = "svd")

  glmWALSoriginal <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                              selfemp + afam, family = binomialWALS(), data = HMDA,
                            prior = weibull(), iterate = TRUE, method = "original")

  expect_equal(coef(glmWALSsvd), coef(glmWALSoriginal), tolerance = tol)
})

test_that("residuals.walsGLM return correct values", {
  tol <- 1e-6 # relative tolerance for deviations
  data("NMES1988", package = "AER")
  NMES1988 <- na.omit(NMES1988)

  fitPoisson <- walsGLM(emergency ~ health + chronic + age + gender |
                          I((age^2)/10) + married + region, family = poissonWALS(),
                        data = NMES1988, prior = laplace(), iterate = TRUE,
                        verbose = TRUE)
  y <- NMES1988$emergency
  mu <- predict(fitPoisson, newdata = NMES1988)
  pearson <- (y - mu) / sqrt(mu)
  devres <- mu # for y = 0
  devres[y > 0] <- (y * (log(y) - log(mu)) - (y - mu))[y > 0]
  devres <- 2 * devres
  devres <- sqrt(devres)
  devres <- ifelse(y - mu > 0, devres, -devres)
  resids <- y - mu

  expect_equal(pearson, residuals(fitPoisson, type = "pearson"), tolerance = tol)
  expect_equal(devres, residuals(fitPoisson, type = "deviance"), tolerance = tol)
  expect_equal(resids, residuals(fitPoisson, type = "response"), tolerance = tol)
})

test_that("Predictions use correct link values", {
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  fitFormula <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                          selfemp + afam, data = HMDA, family = binomialWALS(),
                        prior = weibull(), keepX = TRUE, keepY = TRUE)
  linkFormula <- as.vector(model.matrix(fitFormula, "focus") %*% coef(fitFormula, "focus")
                           + model.matrix(fitFormula, "aux") %*% coef(fitFormula, "aux"))

  fitMatrix <- walsGLM(fitFormula$x$focus, fitFormula$x$aux, fitFormula$y,
                       family = binomialWALS(), prior = weibull())
  linkMatrix <- as.vector(fitFormula$x$focus %*% coef(fitMatrix, "focus")
                          + fitFormula$x$aux %*% coef(fitMatrix, "aux"))

  expect_identical(linkFormula, as.vector(predict(fitFormula, newdata = HMDA,
                                                  type = "link")))
  expect_identical(linkMatrix, as.vector(predict(fitMatrix, newX1 = fitFormula$x$focus,
                                                 newX2 = fitFormula$x$aux, type = "link")))
})

test_that("logLik.walsGLM returns correct value", {
  tol <- 1e-6
  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)

  fitFormula <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                          selfemp + afam, data = HMDA, family = binomialWALS(),
                        prior = weibull(), keepX = TRUE, keepY = TRUE)

  X <- cbind(model.matrix(fitFormula, "focus"), model.matrix(fitFormula, "aux"))
  mu <- plogis(X %*% coef(fitFormula))
  ll <- sum(dbinom(fitFormula$y, size = 1, prob = mu, log = TRUE))

  expect_equal(logLik(fitFormula), ll, tolerance = tol)
})

test_that("Probability predictions work as intended", {
  data("NMES1988", package = "AER")
  NMES1988 <- na.omit(NMES1988)
  fitPoisson <- walsGLM(emergency ~ health + chronic + age + gender |
                          I((age^2)/10) + married + region, family = poissonWALS(),
                        data = NMES1988, prior = laplace())
  # expect no error in probability prediction.
  expect_error(predict(fitPoisson, newdata = NMES1988, type = "prob"),
               regexp = NA)

  data("HMDA", package = "AER")
  HMDA <- na.omit(HMDA)
  fitBinomial <- walsGLM(deny ~ pirat + hirat + lvrat + chist + mhist + phist |
                           selfemp + afam, data = HMDA, family = binomialWALS(),
                         prior = subbotin())
  # expect error in probability prediction
  expect_error(predict(fitBinomial, newdata = HMDA, type = "prob"),
               regexp = "Probability predictions of counts not supported for")
})
