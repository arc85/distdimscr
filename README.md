
<!-- README.md is generated from README.Rmd. Please edit that file -->

# distdimscr

<!-- badges: start -->
<!-- badges: end -->

The general goal of distdimscr is to quantify the differences in cell
populations between conditions in single-cell RNAseq analysis. For
example, one may want to quantify the differences between B cells and
CD4+ T cells in peripheral blood versus B cells and CD4+ T cells in a
tonsil. distdimscr quantifies these differences by measuring the
distance between cells populations in a high-dimensional space (e.g. in
principal component space).

## Installation

distdimscr is currently only available on [GitHub](https://github.com/)
and can be install with:

``` r
# install.packages("devtools")
devtools::install_github("arc85/distdimscr")
```

## Example use

Here is a super basic example as a placeholder until I can cook up
something better. This creates a few matrices and calculate the
Bhattacharrya distance between subsets of cells in each embedding.

By analogy, let’s pretend that these are principal componets derived
from scRNAseq data. Will include an example later…

``` r
library(distdimscr)

set.seed("0222")
mat1 <- matrix(data=rnorm(100000,mean=1,sd=1),nrow=2000,ncol=50)
mat2 <- matrix(data=rnorm(100000,mean=2,sd=1),nrow=2000,ncol=50)
mat3 <- matrix(data=rnorm(100000,mean=4,sd=1),nrow=2000,ncol=50)

## Compare mat1 to mat2
dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
#>          [,1]
#> [1,] 1.859134

# Setting random_sample=TRUE let's you compute a background distance that is independent of the actual conditions.
dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
#>           [,1]
#> [1,] 0.2293587
# Notice that it's much less than the result with random_sample=TRUE

## Comparing mat1 vs mat3 - distance is larger compared with mat1 vs mat2
dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
#>          [,1]
#> [1,] 12.01026
dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
#>           [,1]
#> [1,] 0.1552183
# Notice with random_sample=TRUE the result is similar to the initial result

## Comparing mat2 vs mat3 - distance will be smaller than mat1 vs mat3
dim_dist(embed_mat_x=mat2,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
#>          [,1]
#> [1,] 4.782157
dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
#>           [,1]
#> [1,] 0.2201744
# Notice with random_sample=TRUE the result is similar once again
```

## Future use

If we were actually using this function, we would want to do many
subsamples of cells across conditions to get a sense for the true
Bhattacharrya distance. This will be updated soon with a small use cause
from scRNAseq data.
