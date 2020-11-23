# Multi-snpnet-Cox (mrcox) Efficient Group-Sparse Lasso solver for multi-response Cox model
This package provide functions that solve the following multi-response, group-sparse regularized Cox regression:
<img src="https://latex.codecogs.com/gif.latex?\min_{\beta_1,\cdots,&space;\beta_K}&space;\sum_{k=1}^K&space;\frac{1}{n_k}&space;\left[\sum_{i:O_i^k&space;=&space;1}&space;-\beta_k^T&space;X_i&space;&plus;&space;\log&space;\left(\sum_{j:T^k_j&space;\ge&space;T^k_i}&space;\exp(\beta_k^T&space;X_j)\right)\right]&space;&plus;&space;\lambda&space;\left(&space;\sum_{j=1}^d&space;\|\beta^j\|_1&space;&plus;&space;\alpha&space;\|\beta^j\|_2&space;\right)." title="\min_{\beta_1,\cdots, \beta_K} \sum_{k=1}^K \frac{1}{n_k} \left[\sum_{i:O_i^k = 1} -\beta_k^T X_i + \log \left(\sum_{j:T^k_j \ge T^k_i} \exp(\beta_k^T X_j)\right)\right] + \lambda \left( \sum_{j=1}^d \|\beta^j\|_1 + \alpha \|\beta^j\|_2 \right)." />,

where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta_k" title="\beta_k" /> is the coefficient vector of the kth (out of K) responses, and  <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta^j" title="\beta^j" /> are the coefficients of the jth variable on the K responses.

For genetics data in [Plink2](https://www.cog-genomics.org/plink/2.0/) format, we provide a screening procedure similar to the one in [this paper](https://journals.plos.org/plosgenetics/article?rev=2&id=10.1371/journal.pgen.1009141).

## Installation
Currently mrcox only supports linux and intel processors. It requires Intel's Math Kernel Library (MKL). I suspect you can still get it to run on AMD processors but performance could be significantly worse. To install MKL:
1. Register and download MKL from https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html
2. Under Choose Product to Download, select Intel Math Kernel Library for linux
3. Untar the downloaded file, run `./install.sh`. After installation is done, `source intel/mkl/bin/mklvars.sh intel64`. 
Now install the R dependencies of mrcox (Rcpp, RcppEigen). 
### Additional dependencies for genetics data
1.[zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd)
2. [Plink2](https://www.cog-genomics.org/plink/2.0/)
3. 
```r
library(devtools)
install_github("chrchang/plink-ng", subdir="/2.0/cindex")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```

