# WALS 0.2.4

* Added option `postmult` for `walsGLMfit()`. Corresponds to eq. (9) of De Luca et al. (2018).
Had to add `postmult` option for `walsFit()` as well, but set to `FALSE` by default.

* Fixed GitHub issue #1: Fixed tests comparing results to literature for PowerPC.
Replaced `expect_identical()` with `expect_equal()` using small `tolerance`.


# WALS 0.2.3

* First CRAN release
