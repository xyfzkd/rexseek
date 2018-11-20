#matrix_path = 'inst/extdata/scirep_sequential_qc.txt'
#classinfo_path = 'scirep_classes.txt'

#' @title read counts matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()]
#'
#' @return integer matrix
#' @export
#'
#' @examples NULL
read_mat <- function(path, ...) {
	path %>% readr::read_tsv(T, readr::cols(.default = 'c'), ...) %>%
		dplyr::mutate_at(-1, readr::parse_integer) %>%
		dplyr::rename('transcript' = 1) %>% as.data.frame() %>%
		tibble::column_to_rownames('transcript') %>% as.matrix()
}

#' @title filter genes with low expression values
#'
#' @param mat integer matrix.
#' @param min_count, min_sample_per_gene integer scalar. For each gene, it must contain at least `min_count` reads in at least `min_sample_per_gene` samples. Otherwise, it would be dropped.
#'
#' @return integer matrix.
#' @export
#'
#' @examples NULL
filter_low <- function(mat, min_count = 2, min_sample_per_gene = 5) {
	low_per_row <- rowSums(mat > min_count)
	keeped_row <- low_per_row > min_sample_per_gene
	mat[keeped_row, ]
}



#' @export
plot_highest_exprs <- function(sce) {
	sce %>% {suppressMessages(scater::calculateQCMetrics(.))} %>%
		scater::plotHighestExprs(n = 20)
}

#' @title plot PCA, TSNE
#'
#' @param sce A SingleCellExperiment object.
#' @param shape, color string. specify a column in `col_data` of [as_SingleCellExperiment()] to shape/color by
#'
#' @export
plot_PCA <- function(sce, shape = NULL, color = NULL) {
	sce %>% scater::plotPCA(
		shape_by = shape, colour_by = color,
    	run_args = list(exprs_values = 'counts')
	)
}



#' @describeIn plot_PCA
#'
#' @export
plot_TSNE <- function(sce, shape = NULL, color = NULL) {
	sce %>% scater::plotTSNE(
	    shape_by = shape, colour_by = color,
	    run_args = list(exprs_values = "counts")
	)
}





