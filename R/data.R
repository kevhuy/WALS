#' Determinants of Economic Growth
#'
#' Growth regression data used in \insertCite{magnus2010growth;textual}{WALS}.
#'
#' @format A data frame with 72 observations on 11 variables:
#' \describe{
#' \item{country}{factor. Name of the country.}
#' \item{gdpgrowth}{Average growth rate of GDP per capita from 1960 - 1996 at
#' purchasing power parity.}
#' \item{lgdp60}{Logarithm of GDP per capita in 1960.}
#' \item{equipinv}{Average real equipment investment share of GDP from
#' 1960 - 1985 comprising investments in electrical and nonelectrical machinery
#' (in relative prices constant across countries).}
#' \item{school60}{Enrollment rate for primary education in 1960.}
#' \item{life60}{Life expectancy at age 0 in 1960.}
#' \item{popgrowth}{Average growth rate of population from 1960 - 1996.}
#' \item{law}{Index for the overall maintenance of the rule of law ('law and
#' order tradition').}
#' \item{tropics}{Proportion of country's land area within geographical tropics.}
#' \item{avelf}{Average of five different indices of ethnolinguistic
#' fragmentation which is measured as the probability of two random people
#' in a country not sharing the same language.}
#' \item{confucian}{Fraction of Confucian population in 1970 and 1980.}
#' }
#'
#' @details The dataset is used in \insertCite{magnus2010growth;textual}{WALS}
#' to illustrate the WALS model averaging approach and combines the data used in
#' \insertCite{sala2004bace;textual}{WALS} and \insertCite{sala1997reg;textual}{WALS}.
#' See the references for more detailed descriptions and original sources of the
#' variables.
#'
#' @source WALS package for MATLAB (and Stata) provided on Jan Magnus' personal
#' website.
#' \url{https://www.janmagnus.nl/items/WALS.pdf}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ## Replicate Table 2 in Magnus et al. (2010)
#' # NOTE: prescale = FALSE, still used old version of WALS in Magnus et al. (2010).
#' # Not recommended anymore!
#' fitMPP <- wals(gdpgrowth ~ lgdp60 + equipinv + school60 + life60 + popgrowth |
#'                   law + tropics + avelf + confucian, data = GrowthMPP,
#'                 prior = laplace(), prescale = FALSE)
#' tableMPP <- cbind("coef" = coef(fitMPP), "se" = sqrt(diag(vcov(fitMPP))))
#' print(round(tableMPP, 4))
#'
"GrowthMPP"

#' Determinants of Economic Growth
#'
#' Growth regression data used in \insertCite{masanjala2008growth;textual}{WALS}.
#'
#' @format A data frame with 37 observations on 25 variables:
#' \describe{
#' \item{gdpgrowth}{Average growth rate of GDP per capita from 1960 - 1992 at
#' purchasing power parity.}
#' \item{lgdp60}{Logarithm of GDP per capita in 1960.}
#' \item{yrsopen}{Fraction of years economy open from 1960 - 1990.}
#' \item{mining}{Fraction of GDP in mining.}
#' \item{primexp70}{Share of exports of primary products in GDP in 1970.}
#' \item{invest}{Ratio of real domestic investment (public and private) to real GDP.}
#' \item{rerd}{Real exchange rate distortion.}
#' \item{school60}{Average years of primary schooling for population over 25
#' years of age in 1960.}
#' \item{life60}{Life expectancy at age 0 in 1960.}
#' \item{popgrowth}{Average growth rate of population from 1960 - 1990.}
#' \item{war}{factor. \code{"yes"} if country participates in at least one external
#' war from 1960 to 1985. \code{"no"} else.}
#' \item{revcoup}{Average number of revolutions and coups per year from 1960 - 1990.}
#' \item{rights}{Index of political rights ranging from 1 (most restrictive)
#' to 7 (most freedom)}
#' \item{civil}{Index of civil liberties ranging from 1 (most restrictive)
#' to 7 (most freedom)}
#' \item{out}{Index of outward orientation.}
#' \item{capitalism}{Degree of capitalism.}
#' \item{colony}{factor. Shows if the country used to be \code{"british"} or
#' \code{"french"} colony. If neither of them applies, then \code{"none"}.}
#' \item{english}{Fraction of English speakers.}
#' \item{foreign}{Fraction speaking foreign language.}
#' \item{frac}{Probability that two random people are from different
#' ethnolinguistic groups.}
#' \item{protestant}{Fraction of population Protestant.}
#' \item{catholic}{Fraction of population Catholic.}
#' \item{muslim}{Fraction of population Muslim.}
#' \item{area}{Size of country in millions of square kilometers.}
#' \item{abslat}{Distance from the equator.}
#' }
#'
#' @details
#' The dataset of \insertCite{masanjala2008growth;textual}{WALS} is a subset of
#' sub-Sahara African countries from the data used in
#' \insertCite{sala1997reg;textual}{WALS}. See Table A2. in
#' \insertCite{masanjala2008growth;textual}{WALS} for the original sources of the
#' variables. This dataset is also used for replication purposes in
#' \insertCite{amini2012replication;textual}{WALS}.
#'
#' To replicate the WALS estimates in \insertCite{amini2012replication;textual}{WALS},
#' use all variables except for a constant as auxiliary regressors and divide all
#' regressors by their in-sample maximum before running
#' \code{wals(..., prescale = FALSE)} (\bold{NOTE: It is not recommended to use
#' \code{prescale = FALSE} as this runs an old version of the WALS estimator,
#' \code{prescale = FALSE} should only be used for replication purposes}).
#' The resulting coefficients and standard errors have to be divided by the maximum
#'of the regressors again to get the values presented in Table I of the paper.
#'
#' @source Journal of Applied Econometrics Data Archive.
#' The data was taken from the archive entry of
#' \insertCite{amini2012replication;textual}{WALS} for replication purposes but
#' they can also be found in the archive entry of
#' \insertCite{masanjala2008growth;textual}{WALS}.
#'
#' \url{https://journaldata.zbw.eu/dataset/comparison-of-model-averaging-techniques-assessing-growth-determinants}
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ## Replicate second panel of Table I in Amini & Parmeter (2012)
#' ## NOTE: Authors manually scale data, then rescale the resulting coefs and se.
#' X <- model.matrix(gdpgrowth ~ ., data = GrowthMP)
#' scaleVector <- apply(X, MARGIN = 2, max)
#' Xscaled <- apply(X, MARGIN = 2, function(x) x/max(x))
#' Xscaled <- Xscaled[,-1]
#' datscaled <- as.data.frame(cbind(gdpgrowth = GrowthMP$gdpgrowth, Xscaled))
#'
#' fitMP <- wals(gdpgrowth ~ 1 | ., data = datscaled, prescale = FALSE,
#'               prior = laplace(), eigenSVD = FALSE)
#' tableMP <- cbind("coef" = coef(fitMP)/scaleVector,
#'                  "se" = sqrt(diag(vcov(fitMP)))/scaleVector)
#' printVars <- c("(Intercept)", "lgdp60", "yrsopen", "mining", "primexp70",
#'                "invest")
#' print(round(tableMP[printVars,], 4))
#'
"GrowthMP"
