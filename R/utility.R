as_SingleCellExperiment <- function(mat, info = NULL) {
	assays = list(counts = mat)
	if (is.null(info))
		SingleCellExperiment::SingleCellExperiment(assays = assays)
	else
		SingleCellExperiment::SingleCellExperiment(assays = assays, colData = info)
}