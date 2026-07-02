# ppindices

`ppindices` is an R package for quantifying phenotypic plasticity — the extent to
which a genotype's trait values change across environments — from reaction-norm
data. It collects a broad range of plasticity indices from the literature under
a single, consistently documented interface, covering variance-based, slope-based,
distance-based, and index-based measures, as well as convenience wrappers for
computing indices across multiple genotypes stored in data frames.

The package currently provides 43 exported plasticity indices.

## Installation

```r
# install.packages("remotes")
remotes::install_local("~/ppindices")
```

## Usage

```r
library(ppindices)

calculate_CVt(c(2, 4, 6, 8))
```

See `vignette("computing-plasticity-indices", package = "ppindices")` for a
worked example computing several indices on a toy reaction norm, and
`help(package = "ppindices")` for the full list of available indices.

## License

MIT
