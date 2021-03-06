
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

distdimscr is available on [GitHub](https://github.com/) and can be
install with:

``` r
# install.packages("devtools")
# devtools::install_github("arc85/distdimscr")
```

## Example use

Why are we interested in using the Bhattacharrya distance to measure
differences in high-dimensional space? In high-dimensional space, the
intuitive notion of Euclidean distance between points breaks down,
necessiting a different metric to measure distances. The Bhattacharrya
distance can overcome this problem by measuring the distance between two
non-normal probability distributions in high-dimensional space. See
[Aggarwal et
al](https://link.springer.com/chapter/10.1007/3-540-44503-X_27) for
further reading.

Here, we outline a basic use case for distdimscr. We have two sets of 3
samples, with the first set derived from the peripheral blood of healthy
donors and the second set derived from tonsil tissues from patients
undergoing tonsillectomy. Let’s say we want to quantify the difference
in transcriptional signatures between the same immune cells in
peripheral blood and tonsils (e.g. how different are CD4+ T cells in
peripheral blood versus tonsil). distdimistscr lets us readily quantify
these similiatires and differences between populations in these
different tissues as outline below.

The Bhattacharrya distance approach has been implemented in several
single-cell RNAseq papers, first by [Azizi et al, Cell
2018](https://pubmed.ncbi.nlm.nih.gov/29961579/) and also by [Cillo et
al, Immunity 2020](https://pubmed.ncbi.nlm.nih.gov/31924475/).

``` r
# Load distdimscr
library(distdimscr)
#> Loading required package: Seurat
library(ggplot2)

# Check out UMAP of peripheral blood and tonsil with cell types identified
# Have a look at data-raw for sample acquisition and pre-processing

overall.data <- cbind(overall.umap,overall.metadata)

ggplot(overall.data,aes(x=UMAP_1,y=UMAP_2,colour=cell_types)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~sample_type)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
# Check out cell numbers in each sample
knitr::kable(table(overall.data$sample_type,overall.data$cell_types))
```

|        | B cells | CD14 monocytes | CD1C DCs | CD4 T cells | CD8 T cells | NK cells | pDCs | Plasmablasts | RBCs |
|--------|--------:|---------------:|---------:|------------:|------------:|---------:|-----:|-------------:|-----:|
| PBMC   |     279 |           1173 |       73 |        3049 |        1514 |      493 |   37 |            7 |   23 |
| Tonsil |    4573 |              4 |       41 |        3725 |         267 |       57 |   31 |          176 |    2 |

``` r
# We should only compare cells that are present in both samples
# We will keep B cells, CD4 cells, and CD8 cells
b.cells.tonsil <- rownames(overall.data)[overall.data$cell_types=="B cells" & overall.data$sample_type=="Tonsil"]
b.cells.pbmc <- rownames(overall.data)[overall.data$cell_types=="B cells" & overall.data$sample_type=="PBMC"]

# We have pre-extracted the PCA embeddings from our pre-processed Seurat object
# Let's subset to the cell types identified above
tonsil.b.cells.pca <- overall.pca[b.cells.tonsil,]
pbmc.b.cells.pca <- overall.pca[b.cells.pbmc,]

# Compare tonsil B cells and PBMC B cells - subsample 100 times

bhatt.dist <- bhatt.dist.rand <- vector("logical",length=100)
set.seed("0222")

for (i in 1:100) {

  bhatt.dist[[i]] <- dim_dist(embed_mat_x=tonsil.b.cells.pca,embed_mat_y=pbmc.b.cells.pca,dims_use=1:10,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE)

  bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x=tonsil.b.cells.pca,embed_mat_y=pbmc.b.cells.pca,dims_use=1:10,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=TRUE)

}

# Combine the results and plot
bhatt.dist <- data.frame(B.cells.distance=bhatt.dist,comparison="real")
bhatt.dist.rand <- data.frame(B.cells.distance=bhatt.dist.rand,comparison="random")

bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)

ggplot(bhatt.res,aes(x=comparison,y=B.cells.distance)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=0.5) +
  theme_bw() +
  xlab("Comparison type") +
  ylab("Bhattacharrya distance")
```

<img src="man/figures/README-example-2.png" width="100%" />

## Recommendations for use

When selecting principal components for inclusison, it is best to select
those that explain a signifcant amount of the variance. Here, we
selected 10 as a simple use case. Also the number of cells to subset per
sample can be thought of as a hyperparameter. While not strictly
necessary to subsample, it gives a sense of the underlying distributions
that contribute to the high-dimensional differences in distance between
the samples.

## Future features

In further iterations of this package, we could include additional
functions for direct interaction with Seurat objects and easy ways to
measure the distances between multiple cell types.
