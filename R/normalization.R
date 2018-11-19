
#' @export
norm_SCnorm <- function(mat) {
	mat %>% {SCnorm::SCnorm(., Conditions = rep(1, ncol(.)))} %>%
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
	t(t(mat*1e6) / colSums(mat[row, ], na.rm = T))
}

# norm_cpm ------------------
#' @title CPM normalization
#'
#' @param mat integer matrix. counts. first part (separated by `|`) of row names must be transcript id
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
#' @param transcript_type character.
#'
#' @export
norm_cpm_rm <- function(mat, transcript_type) {
	if (!all(transcript_type %in% rna_type$transcript_type))
		stop('unknown transcript type to remove for CPM normalization')

	id_rm <- tibble::tibble(transcript_type) %>%
		dplyr::left_join(rna_type, by = 'transcript_type') %>% .$transcript_id
	row_rm <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% id_rm}

	norm_cpm_impl(mat, !row_rm)
}



#' @describeIn norm_cpm normalize by of given reference genes
#' @param reference_transcript_id character.
#'
#' @export
norm_cpm_refer <- function(mat, reference_transcript_id) {
	row_refer <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
		{. %in% reference_transcript_id}
	if (length(row_refer) == 0L)
		stop('can\'t find any reference transcript in the matrix for CPM normalization')

	norm_cpm_impl(mat, row_refer)
}
#
# row_refer <- mat %>% rownames() %>% {. %in% reference_transcript_id}
# mat_cpm_refer <- t(t(mat*1e6) / colSums(mat[row_refer, ], na.rm = T))
#
# mat_cpm_refer %T>% SCnorm::plotCountDepth() %>% colMeans() %>% head()




# mi_pi_id <- annotation %>%
# 	dplyr::filter(transcript_type %in% c('miRNA', 'piRNA')) %>% .$transcript_id
