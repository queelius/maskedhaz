# cran-comments

## Submission notes

This is an initial CRAN submission of `maskedhaz` 0.1.0.

`maskedhaz` provides likelihood-based inference for series systems with masked
component cause of failure, using arbitrary dynamic failure rate component
distributions from the `flexhaz` and `serieshaz` packages. It implements the
`series_md` protocol defined in the `maskedcauses` package (already on CRAN)
for cases that have no closed-form likelihood, via numerical integration and
numerical differentiation.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3. 0 errors, 0 warnings, 0 notes
  (other than the sandbox "unable to verify current time" note that does not
  appear on CRAN servers).
* win-builder (R 4.5.3 release): 1 NOTE, the expected "New submission" item.
* rhub v2 (linux and windows, both R-devel): 0 errors, 0 warnings, 0 notes.

## About the NOTE

The R-release win-builder check returns 1 NOTE:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Alexander Towell <lex@metafunctor.com>'
    New submission

This is expected for a first-time submission. The "DFR" term that earlier
win-builder runs flagged as possibly misspelled has been rephrased as
"dynamic failure rate" in the Description field. The `Remotes` field that
earlier win-builder runs noted has been removed.

## Note on R-devel transient binary availability

The win-builder R-devel run (R 4.6.0 beta, r89874) currently errors with
"Package required but not available: 'serieshaz'" even though `serieshaz` is
on CRAN. This is a transient upstream issue: R 4.6 beta had an ABI change
(affecting the `SETLENGTH` symbol, and hence the `rlang` package) that
temporarily broke several Windows binaries on the R-devel mirror, including
`serieshaz`'s. `maskedhaz` itself builds cleanly wherever its Imports are
available (confirmed on linux R-devel and windows R-release). This resolves
as the CRAN R-devel mirror rebuilds.

## Reverse dependencies

None; this is a new package.

## Notes on dependencies

All Imports (`serieshaz`, `flexhaz`, `algebraic.dist`, `likelihood.model`,
`maskedcauses`, `generics`, `numDeriv`, `stats`) are on CRAN.

## Vignettes

Four vignettes build cleanly in about 10 seconds total on a standard dev
machine. They demonstrate: the package overview and ecosystem layering,
mixed-distribution series systems, the four censoring types with
cross-validation against the closed-form reference implementation in
`maskedcauses`, and hypothesis testing on fitted models using the
`hypothesize` package.
