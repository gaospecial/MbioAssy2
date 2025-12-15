# Run Normalized Stochasticity Ratio (NST) Analysis

Run Normalized Stochasticity Ratio (NST) Analysis

## Usage

``` r
run_nst(
  comm,
  group,
  dist_method = "jaccard",
  abundance_weighted = TRUE,
  rand = 20,
  output_file = "NST.output.txt"
)
```

## Arguments

- comm:

  Community abundance matrix (rows are samples, columns are OTUs).

- group:

  Group information (data frame with row names matching samples).

- dist_method:

  Distance method, default "jaccard".

- abundance_weighted:

  Logical, default TRUE.

- rand:

  Number of randomizations, default 20.

- output_file:

  Output file path for the summary.
