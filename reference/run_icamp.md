# Run iCAMP Analysis

Run iCAMP Analysis

## Usage

``` r
run_icamp(
  comm,
  tree,
  treat,
  env = NULL,
  clas,
  prefix = "iCAMP_Test",
  save_dir = tempdir(),
  rand = 1000,
  nworker = 4,
  memory_G = 50,
  bin_size_limit = 24
)
```

## Arguments

- comm:

  Community abundance matrix (rows are samples, columns are OTUs).

- tree:

  Phylogenetic tree object (phylo).

- treat:

  Group information (data frame).

- env:

  Environmental variables (data frame, optional).

- clas:

  Taxonomy information (data frame).

- prefix:

  Prefix for output files.

- save_dir:

  Directory to save outputs and intermediate files.

- rand:

  Randomization times, default 1000.

- nworker:

  Number of threads, default 4.

- memory_G:

  Memory limit in GB, default 50.

- bin_size_limit:

  Minimum bin size, default 24.

## Value

List containing iCAMP results (icres and icbin).
