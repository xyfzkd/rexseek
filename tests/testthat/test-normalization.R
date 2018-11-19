



testthat::test_that('norm_tmm()', {
	mat_tmm <- norm_tmm(mat1000)

	testthat::expect_true(is.matrix(mat_tmm))
	testthat::expect_identical(typeof(mat_tmm), 'numer')
	testthat::expect_identical(ncol(mat_raw), 191L)
})