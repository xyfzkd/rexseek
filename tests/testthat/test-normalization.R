testthat::context('Testing normalization')
if (basename(getwd()) == 'testthat') setwd('../..')

# norm_tmm ----------------------
testthat::test_that('norm_tmm()', {
	mat_tmm <- norm_tmm(mat1000)

	testthat::expect_true(is.matrix(mat_tmm))
	testthat::expect_identical(dim(mat_tmm), dim(mat1000))
})

# norm_rle ----------------------
testthat::test_that('norm_rle()', {
	mat_rle <- norm_rle(mat1000)

	testthat::expect_true(is.matrix(mat_rle))
	testthat::expect_identical(dim(mat_rle), dim(mat1000))
})

# norm_cpm ----------------------
testthat::test_that('norm_cpm()', {
	mat_cpm <- norm_cpm(mat1000)

	testthat::expect_true(is.matrix(mat_cpm))
	testthat::expect_identical(dim(mat_cpm), dim(mat1000))
})

# norm_cpm_top ----------------------
testthat::test_that('norm_cpm_top()', {
	mat_cpm_top <- norm_cpm_top(mat1000, 20)

	testthat::expect_true(is.matrix(mat_cpm_top))
	testthat::expect_identical(dim(mat_cpm_top), dim(mat1000))
})

testthat::test_that('norm_cpm_top() error', {
	norm_cpm_top(mat1000, 1000)

	testthat::expect_error(
		norm_cpm_top(mat1000, 1001),
		'two few feature for CPM top k normalization'
	)
})

# norm_cpm_rm -----------------
testthat::test_that('norm_cpm_rm()', {
	mat_cpm_rm <- norm_cpm_rm(mat1000, c('miRNA', 'piRNA'))

	testthat::expect_true(is.matrix(mat_cpm_rm))
	testthat::expect_identical(dim(mat_cpm_rm), dim(mat1000))
})


testthat::test_that('norm_cpm_rm()', {
	norm_cpm_rm(mat1000, 'lncRNA')

	testthat::expect_error(
		norm_cpm_rm(mat1000, 'noRNA'),
		'unknown transcript type to remove for CPM normalization'
	)
})

# norm_cpm_refer -----------------
testthat::test_that('norm_cpm_refer()', {
	mat_cpm_refer <- norm_cpm_refer(mat1000, suggest_refer_tran_id)

	testthat::expect_true(is.matrix(mat_cpm_refer))
	testthat::expect_identical(dim(mat_cpm_refer), dim(mat1000))
})


testthat::test_that('norm_cpm_refer()', {
	norm_cpm_refer(mat1000, 'lncRNA')

	testthat::expect_error(
		norm_cpm_refer(mat1000, 'noRNA'),
		'unknown transcript type to remove for CPM normalization'
	)
})

