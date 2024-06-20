# WALS 0.2.5

* Allow one-part formulas in all fitting functions: Considers all specified
regressors as auxiliary regressors.
* Exporting all S3 methods due to update to roxygen2 7.3.1:
Internal function `predictCounts()` is now also exported.
* Updated references, CITATION, DESCRIPTION and README.md.


# WALS 0.2.4

* Fixed bug where prediction with only a single auxiliary regressor was not possible: 
Do not drop columns of `X2` in `genNewdata()`.
* Renamed argument `controlGLMfit` in `walsGLM()` to `controlInitGLM` and allow
initialization with restricted model. `controlInitGLM` should be specified using the new 
function `controlGLM()`.
* Updated docs.
* Updated references, CITATION, DESCRIPTION and README.md.
* Added option `postmult` for `walsGLMfit()`. Corresponds to eq. (9) of De Luca et al. (2018).
Had to add `postmult` option for `walsFit()` as well, but set to `FALSE` by default.
* Fixed GitHub issue #1: Fixed tests comparing results to literature for PowerPC.
Replaced `expect_identical()` with `expect_equal()` using small `tolerance`.


# WALS 0.2.3

* First CRAN release
