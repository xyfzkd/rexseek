#' @export
as_SingleCellExperiment <- function(mat, col_data = NULL) {
	assays = list(counts = mat)
	if (is.null(col_data))
		SingleCellExperiment::SingleCellExperiment(assays = assays)
	else
		SingleCellExperiment::SingleCellExperiment(assays = assays, colData = col_data)
}