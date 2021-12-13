test_that("dim_dist calculates a high dimensional distance", {

	set.seed("0222")
	mat1 <- matrix(data=rnorm(100000,mean=1,sd=1),nrow=2000,ncol=50)
	mat2 <- matrix(data=rnorm(100000,mean=2,sd=1),nrow=2000,ncol=50)

	res <- dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE)
	res <- as.numeric(round(res,3))

	set.seed("0222")
	res2 <- dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,dims_use=1:10,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=TRUE)
	res2 <- as.numeric(round(res2,3))

  expect_equal(res, 1.585)
	expect_equal(res2, 0.118)
	expect_failure(dim_dist(embed_mat_x=mat1,embed_mat_y=mat2,distance_metric="euclidean"))

})
