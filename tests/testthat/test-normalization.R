

print(ls())

testthat::test_that('norm_tmm()', {
	mat_tmm <- norm_tmm(mat1000)

	testthat::expect_true(is.matrix(mat_tmm))
	testthat::expect_identical(typeof(mat_tmm), 'integer')
	testthat::expect_identical(dim(mat_raw), dim(mat1000))
})