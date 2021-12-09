dim_dist <-  function(embed_mat_x,embed_mat_y,dims_use=1:50,num_cells_sample=100,random_sample=FALSE) {

	if (random_sample==FALSE) {

		samp1 <- sample(1:nrow(embed_mat_x),num_cells_sample)
		samp2 <- sample(1:nrow(embed_mat_y),num_cells_sample)

		mat.x <- embed_mat_x[samp1,dims_use]
		mat.y <- embed_mat_y[samp2,dims_use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- cov(mat.x)
		cov.mat2 <- cov(mat.y)

		comb.cov <- (cov.mat1 + cov.mat2)/2

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))

	} else {

		comb.mat <- rbind(embed_mat_x,embed_mat_y)

		samp1 <- sample(1:nrow(comb.mat),num_cells_sample)
		samp2 <- sample(1:nrow(comb.mat),num_cells_sample)

		mat.x <- comb.mat[samp1,dims_use]
		mat.y <- comb.mat[samp2,dims_use]

		mean.mat1 <- apply(mat.x,2,mean)
		mean.mat2 <- apply(mat.y,2,mean)

		cov.mat1 <- cov(mat.x[,dims_use])
		cov.mat2 <- cov(mat.y[,dims_use])

		comb.cov <- (cov.mat1 + cov.mat2)/2

		(1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))


	}

}
