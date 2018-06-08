# kernrank

[![Travis-CI Build Status](https://travis-ci.org/YunlongJiao/kernrank.svg?branch=master)](https://travis-ci.org/YunlongJiao/kernrank)

`kernrank` is an `R` package implementing kernel functions and kernel methods for analyzing rank data, typically total rankings (or permutations), interleaving and top-k partial rankings, multivariate rankings.

## Installation

```r
# install.packages("devtools")
devtools::install_github("YunlongJiao/kernrank")
```

## Demo

Please refer to GitHub repos [kendallkernel](https://github.com/YunlongJiao/kendallkernel) and [weightedkendall](https://github.com/YunlongJiao/weightedkendall) for various examples of using the package.

## References

> Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. [DOI:10.1109/TPAMI.2017.2719680](https://doi.org/10.1109/TPAMI.2017.2719680)
> 
> Yunlong Jiao, Jean-Philippe Vert. "The Weighted Kendall and High-order Kernels for Permutations." arXiv preprint arXiv:1802.08526, 2018. [arXiv:1802.08526](https://arxiv.org/abs/1802.08526)

## Note

This package is built upon and is a stable and significant extension of the R CRAN package [`RMallow`](https://cran.r-project.org/web/packages/RMallow/index.html).

## Misc

- Author: [Yunlong Jiao](https://yunlongjiao.github.io/)
- License: GPL-3
