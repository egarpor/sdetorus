# sdetorus 0.2.0

* Bug fixes:
  * `mleOptimWrapper()` now computes the upper-bound penalty against `upper`
    instead of `lower`.
  * `approxMleWn2D()` now actually passes its positive-definiteness `region()`
    to `mleOptimWrapper()` (it was previously ignored), and the feasibility
    check no longer contains a term that evaluated to `NA`. The check is now
    robust to fixing a subset of the parameters.
  * `driftMixVm()` now honors its `expTrc` argument (it was hard-coded to 30).
  * `scoreMatchWnBvm()` now returns `c(0, 0, 0)` (as documented) instead of
    `NA` when the score-matching matrix is singular.
  * `mlePde2D()` now wraps the second coordinate with `My` instead of `Mx` when
    computing the closest grid bins, fixing incorrect binning when `Mx != My`.
  * `approxMleWnPairs()` default `start`, `lower` and `upper` are now consistent
    with the `(alpha, mu, sigma, rho)` parameter order (they were written for a
    different order, and `start` was missing `rho`), so the defaults yield a
    valid, correctly parametrized optimization.
  * `dPsTpd()` now normalizes the Shoji--Ozaki density in the multivariate
    (`p >= 2`) branch used for per-observation start points (a
    `p * log(2 * pi)` factor was missing), so it integrates to one as the
    single-start branch already did. This does not affect `psMle()` estimates
    (the missing factor was constant in the parameters).
  * `dTpdWou()` now accepts a single evaluation point supplied as a length-`p`
    vector instead of erroring.
* Efficiency:
  * `rTrajMou()` precomputes the square root of the (constant) conditional
    covariance once instead of recomputing an eigendecomposition at every step
    (about an 8x speed-up for long trajectories). The sampled trajectory is
    identical to the previous implementation for a given RNG state.
* Documentation:
  * Fixed several inconsistencies (e.g. the `weightsLinearInterp1D()` return
    description, an incomplete sentence in the quadrature rules, the
    `rTrajLangevin()` return type, the arXiv id in `scoreMatchWnBvm()`, and the
    "Fokker-Planck" spelling) and added `\seealso` cross-references throughout.
  * Modernized the package-level documentation to the `"_PACKAGE"` sentinel.
* Added a `testthat` unit-test suite.

# sdetorus 0.1.6

* Initial version.

# sdetorus 0.1.7

* Fix the redefinition of arma::mat forwardSweepPeriodicTridiag() as arma::vec raising an "additional issue" on CRAN.

# sdetorus 0.1.8

* Support forthcoming `Rcpp`'s STRICT_R_HEADERS.

# sdetorus 0.1.9

* Drop C++11 requirement to adhere to new CRAN policies.
* Drop `personList()` and `citEntry()`.

# sdetorus 0.1.10

* Drop dependency on `colorRamps` given its possible archival (BR email 2024-03-01).
