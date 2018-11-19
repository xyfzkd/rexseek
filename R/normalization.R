
#' @export
norm_SCnorm <- function(mat) {
	mat %>% {SCnorm::SCnorm(., Conditions = rep(1, ncol(.)))} %>%
		SingleCellExperiment::normcounts()
}

#' @export
norm_tmm <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts() %T>%
		{mode(.) <- 'integer'}
}

#' @export
norm_rle <- function(mat) {
	mat %>% as_SingleCellExperiment() %>%
		{suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
		scater::normalise() %>% SingleCellExperiment::normcounts()%T>%
		{mode(.) <- 'integer'}
}

#norm_tmm(mat_raw)


