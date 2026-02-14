# Decode candidate set matrix from Boolean columns

Extracts Boolean columns matching the pattern `prefix + digits` from a
data frame and assembles them into a logical matrix. This replaces
[`md.tools::md_decode_matrix()`](https://queelius.github.io/md.tools/reference/md_decode_matrix.html)
to avoid the dependency.

## Usage

``` r
decode_candidate_matrix(df, prefix = "x")
```

## Arguments

- df:

  data frame containing candidate set columns

- prefix:

  column name prefix (default `"x"`)

## Value

logical matrix with one column per component, or NULL if no matching
columns found
