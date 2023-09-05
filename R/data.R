#' Determinants of Economic Growth
#'
#' Growth regression data used in \insertCite{magnus2010growth;textual}{WALS}.
#'
#' @usage data("GrowthMPP")
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
#' \item{confuc}{Fraction of Confucian population in 1970 and 1980.}
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
"GrowthMPP"