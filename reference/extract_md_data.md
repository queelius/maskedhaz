# Extract and validate masked data from a data frame

Shared validation logic for all likelihood model methods. Checks that
the data frame is non-empty, required columns exist, decodes the
candidate set matrix, and validates observation types.

## Usage

``` r
extract_md_data(df, lifetime, omega, candset, lifetime_upper = NULL)
```

## Arguments

- df:

  masked data frame

- lifetime:

  column name for system lifetime

- omega:

  column name for observation type

- candset:

  column prefix for candidate set indicators

- lifetime_upper:

  column name for interval upper bound

## Value

list with components: t, t_upper, omega, C, m, n
