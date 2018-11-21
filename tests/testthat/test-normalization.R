testthat::context('Testing normalization')
if (basename(getwd()) == 'testthat') setwd('../..')


# norm_mat ------------------------------------

testthat::test_that('norm_mat() check args', {
	testthat::expect_error(
		norm_mat(tempfile(), 'CPM_top'),
		'top_n must be specified for CPM_top method'
	)
	testthat::expect_error(
		norm_mat(tempfile(), 'CPM_rm'),
		'rm_gene_type must be specified for CPM_rm method'
	)
	testthat::expect_error(
		norm_mat(tempfile(), 'CPM_refer'),
		'refer_gene_id must be specified for CPM_refer method'
	)
})



# SCnorm ----------------------
testthat::test_that('norm_SCnorm()', {
	mat_SCnorm <- suppressMessages(norm_SCnorm(sim_mat*10))

	testthat::expect_true(is.matrix(mat_SCnorm))
	testthat::expect_identical(dim(mat_SCnorm), dim(sim_mat))
})

# norm_tmm ----------------------
testthat::test_that('norm_tmm()', {
	mat_tmm <- norm_tmm(sim_mat)
norm_SCnorm(sim_mat*10)
	testthat::expect_true(is.matrix(mat_tmm))
	testthat::expect_identical(dim(mat_tmm), dim(sim_mat))
})

# norm_rle ----------------------
testthat::test_that('norm_rle()', {
	mat_rle <- norm_rle(sim_mat)

	testthat::expect_true(is.matrix(mat_rle))
	testthat::expect_identical(dim(mat_rle), dim(sim_mat))
})

# norm_cpm ----------------------
testthat::test_that('norm_cpm()', {
	mat_cpm <- norm_cpm(sim_mat)

	testthat::expect_true(is.matrix(mat_cpm))
	testthat::expect_identical(dim(mat_cpm), dim(sim_mat))
})

# norm_cpm_top ----------------------
testthat::test_that('norm_cpm_top()', {
	mat_cpm_top <- norm_cpm_top(sim_mat, 20)

	testthat::expect_true(is.matrix(mat_cpm_top))
	testthat::expect_identical(dim(mat_cpm_top), dim(sim_mat))
})

testthat::test_that('norm_cpm_top() error', {
	norm_cpm_top(sim_mat, nrow(sim_mat))

	testthat::expect_error(
		norm_cpm_top(sim_mat, nrow(sim_mat) + 1),
		'two few feature for CPM top k normalization'
	)
})

# norm_cpm_rm -----------------
testthat::test_that('norm_cpm_rm()', {
	mat_cpm_rm <- norm_cpm_rm(sim_mat, c('miRNA', 'piRNA'))

	testthat::expect_true(is.matrix(mat_cpm_rm))
	testthat::expect_identical(dim(mat_cpm_rm), dim(sim_mat))
})


testthat::test_that('norm_cpm_rm()', {
	norm_cpm_rm(sim_mat, 'lncRNA')

	testthat::expect_error(
		norm_cpm_rm(sim_mat, 'noRNA'),
		'unknown transcript type to remove for CPM normalization'
	)
})

# norm_cpm_refer -----------------
testthat::test_that('norm_cpm_refer()', {
	mat_cpm_refer <- norm_cpm_refer(sim_mat, suggest_refer$id)

	testthat::expect_true(is.matrix(mat_cpm_refer))
	testthat::expect_identical(dim(mat_cpm_refer), dim(sim_mat))
})


testthat::test_that('norm_cpm_refer()', {
	# test only one reference gene
	norm_cpm_refer(sim_mat, rownames(sim_mat)[1])

	testthat::expect_error(
		norm_cpm_refer(sim_mat, 'non-exist-RNA'),
		'can\'t find any reference transcript in the matrix for CPM normalization'
	)
})

