#' @export
as_SingleCellExperiment <- function(mat, col_data = NULL) {
	assays = list(counts = mat)
	if (is.null(col_data))
		SingleCellExperiment::SingleCellExperiment(assays = assays)
	else
		SingleCellExperiment::SingleCellExperiment(assays = assays, colData = col_data)
}


#' @title calculate coefficient of variance
#'
#' @param x numeric.
#'
#' @return numeric scalar.
cv_fun <- function(x) {
	sd(x, na.rm = T) / mean(x, na.rm = T)
}





