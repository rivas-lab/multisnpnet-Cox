# Multi-snpnet-Cox (mrcox) Efficient Group-Sparse Lasso solver for multi-response Cox model
This package provide functions that solve the following multi-response, group-sparse regularized Cox regression:

![equation](https://latex.codecogs.com/svg.latex?%5Cmin_%7B%5Cbeta_1%2C%5Ccdots%2C%26space%3B%5Cbeta_K%7D%26space%3B%5Csum_%7Bk%3D1%7D%5EK%26space%3B%5Cfrac%7B1%7D%7Bn_k%7D%26space%3B%5Cleft%5B%5Csum_%7Bi%3AO_i%5Ek%26space%3B%3D%26space%3B1%7D%26space%3B-%5Cbeta_k%5ET%26space%3BX_i%26space%3B%26plus%3B%26space%3B%5Clog%26space%3B%5Cleft%28%5Csum_%7Bj%3AT%5Ek_j%26space%3B%5Cge%26space%3BT%5Ek_i%7D%26space%3B%5Cexp%28%5Cbeta_k%5ET%26space%3BX_j%29%5Cright%29%5Cright%5D%26space%3B%26plus%3B%26space%3B%5Clambda%26space%3B%5Cleft%28%26space%3B%5Csum_%7Bj%3D1%7D%5Ed%26space%3B%5C%7C%5Cbeta%5Ej%5C%7C_1%26space%3B%26plus%3B%26space%3B%5Calpha%26space%3B%5C%7C%5Cbeta%5Ej%5C%7C_2%26space%3B%5Cright%29.)

Where <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta_k" title="\beta_k" /> is the coefficient vector of the kth (out of K) responses, and  <img src="https://latex.codecogs.com/gif.latex?\inline&space;\beta^j" title="\beta^j" /> are the coefficients of the jth variable on the K responses.

For genetics data in [PLINK2](https://www.cog-genomics.org/plink/2.0/) format, we provide a screening procedure similar to the one in [this paper](https://journals.plos.org/plosgenetics/article?rev=2&id=10.1371/journal.pgen.1009141).

## Installation
Currently mrcox only supports linux and 64-bit intel processors. It requires Intel's Math Kernel Library (MKL). I suspect you can still get it to run on AMD processors but performance could be significantly worse. To install MKL:
1. Register and download MKL from https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html
2. Under Choose Product to Download, select Intel Math Kernel Library for linux
3. Untar the downloaded file, run `./install.sh`. After installation is done, `source intel/mkl/bin/mklvars.sh intel64`. 

You will also need to install the R dependencies of mrcox (Rcpp, RcppEigen). 
### Additional dependencies for genetics data
1. [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd)
2. [PLINK2](https://www.cog-genomics.org/plink/2.0/)
3. 
```r
library(devtools)
install_github("chrchang/plink-ng", subdir="/2.0/cindex")
install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
```

Once dependencies are installed, run the following in R
```r
devtools::install_github("rivas-lab/multisnpnet-Cox")
```
### Example of Usage
```r
library(mrcox)
# Simulate some data
n = 1000
p = 5000
X = matrix(rnorm(n*p), n, p)
y1 = rexp(n) * exp(X %*% rbinom(p, 1, 0.1))
y2 = rexp(n) * exp(X %*% rbinom(p, 1, 0.1))
s1 = rbinom(n, 1, 0.3)
s2 = rbinom(n, 1, 0.3)
y_list = list(y1, y2)
s_list = list(s1, s2)

# Initialize coefficient matrix at 0
B = matrix(0.0, p, 2)

# Compute residual at B
res = get_residual(X, y_list, s_list, B)

# Compute the gradient of B
g = t(X) %*% res

# Get the dual norm and lambda sequence
alpha = sqrt(2)
dnorm = get_dual_norm(g, alpha)
lambda_max = max(dnorm)
lambda_min = 0.001 * lambda_max
lambda_seq = exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = 100))

# Fit a model
fit = solve_aligned(X, y_list, s_list, lambda_seq, lambda_seq * alpha)
fit = fit[['result']] # fit also contains the residuals at each solution


```
