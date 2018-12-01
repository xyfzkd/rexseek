#matrix_path = 'data-raw/external/scirep_sequential_qc.txt'
#classinfo_path = 'scirep_classes.txt'

#' @title Read counts matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()].
#'
#' @return integer matrix. counts.
#'
#' @details
#' In any case, first part (separated by `|`) of row names must be
#' Ensembl transcript id
#'
#' @export

# path = 'data-raw/external/scirep_sequential_qc.txt'
read_mat <- function(path, ...) {
	path %>% readr::read_tsv(T, readr::cols(.default = 'c'), ...) %>%
		dplyr::mutate_at(-1, readr::parse_integer) %>%
		dplyr::rename('transcript' = 1) %>% as.data.frame() %>%
		tibble::column_to_rownames('transcript') %>% as.matrix()
}



#' @title filter genes with low expression values
#'
#' @param mat integer matrix. counts.
#' @param min_count, min_sample_per_gene integer scalar. For each gene, it must
#'   contain at least `min_count` reads in at least `min_sample_per_gene`
#'   samples. Otherwise, it would be dropped.
#'
#' @return integer matrix. counts.
#'
#' @examples
#' filter_low(sim_mat)
#'
#' @export
filter_low <- function(mat, min_count = 2, min_sample_per_gene = 5) {
	low_per_row <- rowSums(mat > min_count)
	keeped_row <- low_per_row > min_sample_per_gene
	mat[keeped_row, ]
}


#' @title Plot highest expressed genes
#'
#' @details Each row represents a gene. In a row, each bar represents a sample,
#' the x axis tells the percentage of counts accounted for across the whole
#' dataset. (counts of that gene in that sample / total counts of all genes in
#' all samples)
#'
#' Genes are sorted by average expression.
#'
#' @param mat numeric matrix. counts.
#' @param top_n integer scalar. How many genes to show. If greater that total
#'   gene number, show all genes.
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_highest_exprs()
#'
#' @export
plot_highest_exprs <- function(mat, top_n = 20) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressMessages(scater::calculateQCMetrics(.))} %>%
		scater::plotHighestExprs(n = top_n)
}


# plot_group --------------

#' @title workhorse of `plot_PCA/TSNE()`
#'
#' @details `plot_PCA(sce, shape, color) -> `plot_group_impl(sce, shape, color, scater::plotPCA)`
#'
#' @keywords internal
plot_group_impl <- function(sce, shape = NULL, color = NULL, plot_fun) {
 	plot_fun(
 		sce,
		shape_by = shape, colour_by = color,
    	run_args = list(exprs_values = 'counts')
	)
}

#' @title plot PCA, TSNE
#'
#' @param sce SingleCellExperiment object.
#' @param shape, color string. Specify a column in `col_data` of
#'   [as_SingleCellExperiment()] to shape/color by.
#'
#' @return ggplot object.
#'
#' @name plot_group
NULL


#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#'
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_PCA <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotPCA)
}



#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#'
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_TSNE <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotTSNE)
}

# plot CV -------------------------

#' @title calculate coefficient of variance
#'
#' @param x numeric.
#'
#' @return numeric scalar.
cv_fun <- function(x) {
	sd(x, na.rm = T) / mean(x, na.rm = T)
}


#' @title density plot of coefficient of variation
#'
#' @param mat integer matrix. counts
#' @param refer_gene_id character. Ensembl transcript id, add a vertical line for each gene to mark the corresponding CV (on x axis). Only genes in counts matrix would be shown. Usually these genes should be the reference genes you want to use for normalization.
#' @param refer_gene_name character. Transcript name
#'
#' @return [ggplot2::ggplot()] object
#'
#' @examples
#'
#'
#' @export

# mat = sim_mat
# refer_gene_id = suggest_refer$id
# refer_gene_name = suggest_refer$name
plot_cv_density <- function(mat, refer_gene_id = '', refer_gene_name = refer_gene_id) {
	cv <- mat %>% apply(1, cv_fun) %>%
		{tibble::tibble(id = names(.), value = .)} %>%
		dplyr::mutate(id = stringr::str_extract(id, '[^|]+'))
	plot <- ggplot2::ggplot(cv) +
		ggplot2::geom_density(ggplot2::aes(value), color = 'blue') +
		ggplot2::labs(x = 'coefficient of variation')

	if (length(refer_gene_id) != length(refer_gene_name)) {
		warning("Ignoring refer_gene_name, since it isn't the same length as refer_gene_id")
		refer_gene_name = refer_gene_id
	}
	cv_refer <- tibble::tibble(id = refer_gene_id, name = refer_gene_name) %>%
		dplyr::inner_join(cv, by = 'id')
	if (nrow(cv_refer) == 0L) return(plot)

	plot +
		ggplot2::geom_vline(xintercept = cv_refer$value, color = 'green') +
		ggplot2::geom_point(
			ggplot2::aes(x = value, y = seq_along(value)),
			data = cv_refer, size = 2, shape = 1
		) +
		ggrepel::geom_label_repel(
			ggplot2::aes(x = value, y = seq_along(value), label = name),
			data = cv_refer, hjust = 0.5
		)
}






