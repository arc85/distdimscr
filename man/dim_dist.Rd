% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dim_dist.R
\name{dim_dist}
\alias{dim_dist}
\title{Calculate the distance between a subset of cells in a high-dimensional embedding}
\usage{
dim_dist(
  embed_mat_x,
  embed_mat_y,
  dims_use = 1:10,
  num_cells_sample = 100,
  distance_metric = "bhatt_dist",
  random_sample = FALSE
)
}
\arguments{
\item{embed_mat_x}{A matrix or data.frame containing high-dimensional embeddings for each cell (e.g. PCs) in condition x. Assumes cells are in rows and embedding dimensions are in columns.}

\item{embed_mat_y}{A matrix or data.frame containing high-dimensional embeddings for each cell (e.g. PCs) in condition y. Assumes cells are in rows and embedding dimensions are in columns.}

\item{dims_use}{Dimensions of the high-dimensional embeddings to use for calculating the distance. Defaults to 10 dimensions.}

\item{num_cells_sample}{Number of cells to subset from the overall embedding matrices. Defaults to 100 cells.}

\item{distance_metric}{Distance metric to calculate distances between cells. Defaults to "bhatt_dist", that is the Bhattacharyya, which is currently the only high dimensional distance metric implemented.}

\item{random_sample}{Whether to sample cells from each condition or to sample cells irrespective of condition to calculate a background distribution. Defaults to random_sample=FALSE to calculate a background distribution.}
}
\value{
A numeric distance value.
}
\description{
Calculate the distance between a subset of cells in a high-dimensional embedding
}
\examples{
# Generate two matrices of 1000 cells each with different distributions
set.seed("0222")

mat1 <- matrix(data=rnorm(100000,mean=1,sd=1),nrow=2000,ncol=50)
mat2 <- matrix(data=rnorm(100000,mean=2,sd=1),nrow=2000,ncol=50)
mat3 <- matrix(data=rnorm(100000,mean=3,sd=1),nrow=2000,ncol=50)
dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)

}
