norm_mat_check_arg <- function(norm_methods, top_n, rm_gene_type, refer_gene_id) {
	if ('CPM_top' %in% norm_methods && is.null(top_n))
		stop('top_n must be specified for CPM_top method')

	if ('CPM_rm' %in% norm_methods && is.null(rm_gene_type))
		stop('rm_gene_type must be specified for CPM_rm method')

	if ('CPM_refer' %in% norm_methods && is.null(refer_gene_id))
		stop('refer_gene_id must be specified for CPM_refer method')

}

norm_mat_rmd <- function() {

}

#' @example
#' \donotrun{
#'     norm_mat(
#'         '/path/to/matrix', c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'),
#'         top_n = 20, rm_gene_type = c('miRNA', 'piRNA'), refer_gene_id = c("ENST00000408438.1", "ENST00000385271.1", "ENST00000607334.3", "ENST00000385059.1", "ENST00000362134.1", "ENST00000385245.1", "ENST00000385045.1", "ENST00000362117.1", "ENST00000384832.1", "ENST00000579846.3")
#'     )
#' }
#' @export
norm_mat <- function(
	counts_mat_path,
	norm_methods = c('SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer'),
	top_n = NULL, rm_gene_type = NULL, refer_gene_id = NULL,
	min_count = 2, min_sample_per_gene = 5,
	highest_exprs_top_n = 20,
	sample_class_path = NULL, PCA_label_by = NULL, PAC_color_by = NULL,
	refer_gene_name = refer_gene_id,
	output_dir = '.', output_file = 'norm'
) {
	norm_mat_check_arg(norm_methods, top_n, rm_gene_type, refer_gene_id)
	mat <- read_mat(counts_mat_path) %>% filter_low(min_count, min_sample_per_gene)

	if ('SCnorm' %in% norm_methods)    mat_SCnorm <- norm_SCnorm(mat)
	if ('TMM' %in% norm_methods)       mat_tmm <- norm_tmm(mat)
	if ('RLE' %in% norm_methods)       mat_rle <- norm_rle(mat)
	if ('CPM' %in% norm_methods)       mat_cpm <- norm_cpm(mat)
	if ('CPM_top' %in% norm_methods)   mat_cpm_top <- norm_cpm_top(mat, top_n)
	if ('CPM_rm' %in% norm_methods)    mat_cpm_rm <- norm_cpm_rm(mat, rm_gene_type)
	if ('CPM_refer' %in% norm_methods) mat_cpm_refer <- norm_cpm_refer(mat, refer_gene_id)

	mat_names <- ls(pattern = '^mat_')

	for(mat_name in mat_names) {
		readr::write_tsv(get(mat_name), paste0(output_dir, '/', output_file, stringr::str_extract(mat_name, '_\\w+', '.tsv')))

	}

}













#' @export
norm_SCnorm <- function(mat, ...) {
	mat %>% {suppressMessages(SCnorm::SCnorm(Conditions = rep(1, ncol(.)), Data = ., ...))} %>%
		SingleCellExperiment::normcounts()
}

#' @title TMM/RLE normalization
#'
#' @param mat integer matrix. counts
#' @name norm_scater
NULL


#' @describeIn norm_scater TMM normalization
#'
#' @export
norm_tmm <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
}

#' @describeIn norm_scater RLE normalization
#' @export
norm_rle <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
}


#' @param mat integer matrix. counts
#' @param row integer or logical. Use which rows (genes) to perfrom normalization
norm_cpm_impl <- function(mat, row) {
	t(t(mat*1e6) / colSums(mat[row, , drop = F], na.rm = T))
}

# norm_cpm ------------------
#' @title CPM normalization
#'
#' @param mat integer matrix. counts.
#'
#' @details some functions may throw errors
#'
#' @name norm_cpm
NULL


#' @describeIn norm_cpm normalize by all genes
#' @export
norm_cpm <- function(mat) {
	row_all <- nrow(mat) %>% seq_len()

	norm_cpm_impl(mat, row_all)
}

#' @describeIn norm_cpm normalize by top k genes sorted by expression level
#' @export
norm_cpm_top <- function(mat, top_n) {
	if (nrow(mat) < top_n)
		stop('two few feature for CPM top k normalization')

	row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
		 {.$ix[seq_len(top_n)]}

	norm_cpm_impl(mat, row_top)
}


#' @describeIn norm_cpm normalize by removing genes of given type (such as piRNAs)
#' @param gene_type character.
#'
#' @export
norm_cpm_rm <- function(mat, gene_type) {
	if (!all(gene_type %in% rexseek::rna_type$gene_type))
		stop('unknown transcript type to remove for CPM normalization')

	id_rm <- tibble::tibble(gene_type) %>%
		dplyr::left_join(rexseek::rna_type, by = 'gene_type') %>% .$transcript_id
	row_rm <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% id_rm}

	norm_cpm_impl(mat, !row_rm)
}



#' @describeIn norm_cpm normalize by of given reference genes
#'
#' @param refer_gene_id character.
#'
#' @examples
#' norm_cpm_refer(sim_mat, suggest_refer$id)
#'
#' @export

# mat = sim_mat
# refer_gene_id = suggest_refer$id
norm_cpm_refer <- function(mat, refer_gene_id) {
	row_refer <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% refer_gene_id}
	if (!any(row_refer))
		stop('can\'t find any reference transcript in the matrix for CPM normalization')

	norm_cpm_impl(mat, row_refer)
}
