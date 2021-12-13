#' Calculate the distance between a subset of cells in a high-dimensional embedding
#'
#' @param embed_mat_x A matrix or data.frame containing high-dimensional embeddings for each cell (e.g. PCs) in condition x. Assumes cells are in rows and embedding dimensions are in columns.
#' @param embed_mat_y A matrix or data.frame containing high-dimensional embeddings for each cell (e.g. PCs) in condition y. Assumes cells are in rows and embedding dimensions are in columns.
#' @param dims_use Dimensions of the high-dimensional embeddings to use for calculating the distance. Defaults to 10 dimensions.
#' @param num_cells_sample Number of cells to subset from the overall embedding matrices. Defaults to 100 cells.
#' @param distance_metric Distance metric to calculate distances between cells. Defaults to "bhatt_dist", that is the Bhattacharyya, which is currently the only high dimensional distance metric implemented.
#' @param random_sample Whether to sample cells from each condition or to sample cells irrespective of condition to calculate a background distribution. Defaults to random_sample=FALSE to calculate a background distribution.
#'
#' @return A numeric distance value.
#' @export
#'
#' @examples
#' # Generate two matrices of 1000 cells each with different distributions
#' set.seed("0222")
#'
#' mat1 <- matrix(data=rnorm(100000,mean=1,sd=1),nrow=2000,ncol=50)
#' mat2 <- matrix(data=rnorm(100000,mean=2,sd=1),nrow=2000,ncol=50)
#' mat3 <- matrix(data=rnorm(100000,mean=3,sd=1),nrow=2000,ncol=50)
#' dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
#' dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
#' dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=FALSE)
#' dim_dist(embed_mat_x=mat1,embed_mat_y=mat3,dims_use=1:10,num_cells_sample=100,random_sample=TRUE)
#'

dim_dist <-  function(embed_mat_x,embed_mat_y,dims_use=1:10,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE) {

	if (random_sample==FALSE) {

		samp1 <- sample(1:nrow(embed_mat_x),num_cells_sample)
		samp2 <- sample(1:nrow(embed_mat_y),num_cells_sample)

		mat.x <- embed_mat_x[samp1,dims_use]
		mat.y <- embed_mat_y[samp2,dims_use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- stats::cov(mat.x)
		cov.mat2 <- stats::cov(mat.y)

		comb.cov <- (cov.mat1 + cov.mat2)/2

		if (distance_metric=="bhatt_dist") {

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))

		} else {

			message <- "Only bhatt_dist is currently available."
			fail(message)

		}

	} else {

		comb.mat <- rbind(embed_mat_x,embed_mat_y)

		samp1 <- sample(1:nrow(comb.mat),num_cells_sample)
		samp2 <- sample(1:nrow(comb.mat),num_cells_sample)

		mat.x <- comb.mat[samp1,dims_use]
		mat.y <- comb.mat[samp2,dims_use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- stats::cov(mat.x[,dims_use])
		cov.mat2 <- stats::cov(mat.y[,dims_use])

		comb.cov <- (cov.mat1 + cov.mat2)/2

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))

		if (distance_metric=="bhatt_dist") {

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))

		} else {

			message <- "Only bhatt_dist is currently available."
			fail(message)

		}

	}

}
