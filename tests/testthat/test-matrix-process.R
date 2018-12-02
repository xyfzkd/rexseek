testthat::context('Testing matrix process')
if (basename(getwd()) == 'testthat') setwd('../..')

# read_mat --------------------
testthat::test_that('read_mat()', {
	temp_file <- tempfile()
	sim_mat %>% tibble::as_tibble(rownames = 'transcript') %>%
		readr::write_tsv(temp_file);
	mat <- read_mat(temp_file)

	testthat::expect_true(is.matrix(mat))
	testthat::expect_identical(typeof(mat), 'integer')
	testthat::expect_identical(dim(mat), c(1000L, 50L))
});


# filter_low --------------------
testthat::test_that('filter_low()', {
	testthat::expect_true(is.matrix(filter_low(sim_mat)))
	testthat::expect_identical(typeof(filter_low(sim_mat)), 'integer')

	testthat::expect_identical(dim(filter_low(sim_mat/4)), c(120L, 50L))
	testthat::expect_identical(
		dim(filter_low(sim_mat, min_count = 7)),
		c(672L, 50L)
	)
	testthat::expect_identical(
		dim(filter_low(sim_mat, min_sample_per_gene = 40)),
		c(912L, 50L)
	)
});


# plot_highest_exprs --------------------
testthat::test_that('plot_highest_exprs()', {
	testthat::expect_true(T);

	plot_highest_exprs(sim_mat)
	plot_highest_exprs(sim_mat[1:10,])  # test `top_n` > `nrow(mat)`
});

# plot_PCA --------------------
testthat::test_that('plot_PCA()', {
	testthat::expect_true(T);

	as_SingleCellExperiment(sim_mat) %>% plot_PCA()
	as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
	as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
});



# plot_TSNE --------------------
testthat::test_that('plot_TSNE()', {
	testthat::expect_true(T);

	as_SingleCellExperiment(sim_mat) %>% plot_TSNE()
	as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_TSNE(shape = 'label')
	as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_TSNE(color = 'label')
});




testthat::test_that('cv_fun()', {
	set.seed(0)

	testthat::expect_equal(cv_fun(rnorm(5)), 0.9238771)
	testthat::expect_equal(cv_fun(rpois(5, 4)), 0.5832981)
});


# plot_cv_density ----------------------
testthat::test_that('plot_cv_density()', {
	plot_cv_density(sim_mat, suggest_refer$id) # only show one which exsits in matrix
	plot_cv_density(sim_mat, suggest_refer$id, suggest_refer$name)
	plot_cv_density(sim_mat, rownames(sim_mat)[1:6])
	plot_cv_density(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))

	#warnings
	testthat::expect_warning(
		plot_cv_density(sim_mat, 'non-exist-RNA'),
		"None refer gene found in the count matrix"
	)
	testthat::expect_warning(
		plot_cv_density(sim_mat, rownames(sim_mat)[1:6], letters[1:3]),
		"Ignoring refer_gene_name, since it isn't the same length as refer_gene_id"
	)
});


# plot_refer_violin ----------------------
testthat::test_that('plot_refer_violin()', {
	plot_refer_violin(sim_mat, suggest_refer$id) # only show one which exsits in matrix
	plot_refer_violin(sim_mat, suggest_refer$id, suggest_refer$name)
	plot_refer_violin(sim_mat, rownames(sim_mat)[1:6])
	plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))

	#warnings
	testthat::expect_warning(
		plot_refer_violin(sim_mat, 'non-exist-RNA'),
		"None refer gene found in the count matrix"
	)
	testthat::expect_warning(
		plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], letters[1:3]),
		"Ignoring refer_gene_name, since it isn't the same length as refer_gene_id"
	)
});



