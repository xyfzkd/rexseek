% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix-process.R
\name{plot_group}
\alias{plot_group}
\alias{plot_PCA}
\alias{plot_TSNE}
\title{plot PCA, TSNE}
\usage{
plot_PCA(sce, shape = NULL, color = NULL)

plot_TSNE(sce, shape = NULL, color = NULL)
}
\arguments{
\item{sce}{SingleCellExperiment object.}

\item{shape, }{color string. Specify a column in \code{col_data} of
\code{\link[=as_SingleCellExperiment]{as_SingleCellExperiment()}} to shape/color by.}
}
\value{
ggplot object.
}
\description{
plot PCA, TSNE
}
\examples{
as_SingleCellExperiment(sim_mat) \%>\% plot_PCA()

as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA()
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(shape = 'label')
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(color = 'label')
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(shape = 'label', color = 'label')

as_SingleCellExperiment(sim_mat) \%>\% plot_PCA()

as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA()
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(shape = 'label')
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(color = 'label')
as_SingleCellExperiment(sim_mat, sim_sample_class) \%>\% plot_PCA(shape = 'label', color = 'label')

}
