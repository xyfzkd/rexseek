testthat::context('Testing matrix process')
if (basename(getwd()) == 'testthat') setwd('../..')

mat_raw <- read_mat('inst/extdata/scirep_sequential_qc.txt', n_max = 1000)

testthat::test_that('read_mat()', {
	testthat::expect_true(is.matrix(mat_raw))
	testthat::expect_identical(typeof(mat_raw), 'integer')
	testthat::expect_identical(dim(mat_raw), c(1000L, 191L))
});



testthat::test_that('filter_low()', {
	testthat::expect_true(is.matrix(filter_low(mat_raw)))
	testthat::expect_identical(typeof(filter_low(mat_raw)), 'integer')

	testthat::expect_identical(dim(filter_low(mat_raw)), c(100L, 191L))
	testthat::expect_identical(
		dim(filter_low(mat_raw, min_count = 10)),
		c(24L, 191L)
	)
	testthat::expect_identical(
		dim(filter_low(mat_raw, min_sample_per_gene = 20)),
		c(58L, 191L)
	)

	dim(filter_low(mat_raw, min_sample_per_gene = 20))
});



testthat::test_that('plot_PCA()', {
	testthat::expect_true(T);

	as_SingleCellExperiment(mat1000) %>% plot_PCA()
	as_SingleCellExperiment(mat1000, sample_class) %>% plot_PCA(shape = 'label')
	as_SingleCellExperiment(mat1000, sample_class) %>% plot_PCA(color = 'label')
});




testthat::test_that('plot_TSNE()', {
	testthat::expect_true(T);

	as_SingleCellExperiment(mat1000) %>% plot_TSNE()
	as_SingleCellExperiment(mat1000, sample_class) %>% plot_TSNE(shape = 'label')
	as_SingleCellExperiment(mat1000, sample_class) %>% plot_TSNE(color = 'label')
});

