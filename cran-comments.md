# cran-comments

## Submission notes

This is an initial CRAN submission of `maskedhaz` 0.1.0.

`maskedhaz` provides likelihood-based inference for series systems with masked
component cause of failure, using arbitrary dynamic failure rate component
distributions from the `flexhaz` / `serieshaz` packages. It implements the
`series_md` protocol defined in the `maskedcauses` package (already on CRAN)
for cases that have no closed-form likelihood, via numerical integration and
numerical differentiation.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3 (x86_64-pc-linux-gnu)
* GitHub Actions CI (tested via r-universe build chain)

## R CMD check results

0 ERRORs, 0 WARNINGs, 0 NOTEs (excluding the sandbox-specific
"unable to verify current time" note that disappears on CRAN servers).

## Reverse dependencies

None; this is a new package.

## Notes on dependencies

All Imports (serieshaz, flexhaz, algebraic.dist, likelihood.model,
maskedcauses, generics, numDeriv, stats) are on CRAN. The `Remotes:` field
in DESCRIPTION is retained to aid non-CRAN dev installs from r-universe
and is stripped by CRAN during build.

## Vignettes

Three vignettes build cleanly in under 10 seconds total on a standard dev
machine. They demonstrate: the package overview and ecosystem layering,
mixed-distribution series systems, and the four censoring types with
cross-validation against the closed-form reference implementation in
`maskedcauses`.
