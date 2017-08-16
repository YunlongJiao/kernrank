# kernrank

[![Travis-CI Build Status](https://travis-ci.org/YunlongJiao/kernrank.svg?branch=master)](https://travis-ci.org/YunlongJiao/kernrank)

`kernrank` is an `R` package implementing kernel functions and kernel methods for analyzing rank data, typically total rankings (or permutations), interleaving and top-k partial rankings, multivariate rankings.

## Installation

```r
# install.packages("devtools")
devtools::install_github("YunlongJiao/kernrank")
```

## Demo

Refer to github repo [kendallkernel_demo](https://github.com/YunlongJiao/kendallkernel_demo) for plenty of examples of using the package.

## References

> Yunlong Jiao, Jean-Philippe Vert. The Kendall and Mallows Kernels for Permutations. 2016. [hal-01279273](https://hal.archives-ouvertes.fr/hal-01279273) 

## Note

This package is built upon the CRAN package [`RMallow`](https://cran.r-project.org/web/packages/RMallow/index.html) and serves as a stable update and mostly a significant extension of that package.

## Misc

- Author: [Yunlong Jiao](http://cbio.ensmp.fr/~yjiao/)
- License: GPL-3
