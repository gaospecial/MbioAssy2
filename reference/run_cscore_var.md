# Run C-score-var Analysis

Run C-score-var Analysis

## Usage

``` r
run_cscore_var(
  comm,
  nReps = 500,
  save_plot = FALSE,
  output_file = "c-score-var.summary.txt",
  plot_file = "c-score-var.hist.emf"
)
```

## Arguments

- comm:

  Abundance matrix (rows are samples, columns are OTUs).

- nReps:

  Number of replications, default 500.

- save_plot:

  Logical, whether to save the histogram plot.

- output_file:

  Output file path for summary.

- plot_file:

  Output file path for plot (EMF format).
