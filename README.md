# Multi-snpnet-Cox (mrcox) Efficient Group-Sparse Lasso solver for multiresponse Cox model

### Installation guide for SCG and Sherlock Cluster
To install this R package on SCG, first load the dependencies through `ml msc/0.1.0`. The install the package in `R`:
```r
install.packages("<where the package is downloaded>", repo=NULL, type="source")
```
After installation, everytime you run `ml purge` or log out from the node you will need to run `ml msc/0.1.0` to get the right version of plink2, zstd.

On Sherlock, one extra step is needed before installing the package. After running `ml msc/0.1.0`, you will also need to load one of the Sherlock provided gcc compilers. I've tested `ml gcc/7.1.0` works fine. This only needs to be done everytime you install the package.

For better performance, use multiple Intel CPUs (all Sherlock and SCG compute modes use Intel CPUs). To ask for multiple CPUs, set the flag --cpus-per-task when submitting the job.

### General Installation guide
Coming soon!