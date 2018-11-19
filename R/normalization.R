
#' @export
norm_SCnorm <- function(mat) {
	mat %>% {SCnorm::SCnorm(., Conditions = rep(1, ncol(.)))} %>%
		SingleCellExperiment::normcounts()
}

#' @export
norm_tmm <- function(mat) {
	mat_tmm <- mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
	mode(mat_tmm) <- 'integer'
}

#' @export
norm_rle <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()
}

#norm_tmm(mat_raw)


