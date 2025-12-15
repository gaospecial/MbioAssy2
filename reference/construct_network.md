# Construct Co-occurrence Network

Construct Co-occurrence Network

## Usage

``` r
construct_network(
  comm,
  cor_cutoff = 0.8,
  p_cutoff = 0.05,
  use_rmt = FALSE,
  method = "spearman",
  rmt_interval = c(0.6, 0.99)
)
```

## Arguments

- comm:

  Abundance matrix (rows are samples, columns are OTUs).

- cor_cutoff:

  Correlation coefficient cutoff (optional if use_rmt is TRUE).

- p_cutoff:

  P-value cutoff.

- use_rmt:

  Logical, whether to use RMT to find cutoff.

- method:

  Correlation method, default "spearman".

- rmt_interval:

  Interval for RMT threshold search.
