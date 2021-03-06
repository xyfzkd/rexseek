% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix-process.R
\name{plot_highest_exprs}
\alias{plot_highest_exprs}
\title{Plot highest expressed genes}
\usage{
plot_highest_exprs(mat, top_n = 20)
}
\arguments{
\item{mat}{numeric matrix. counts.}

\item{top_n}{integer scalar. How many genes to show. If greater that total
gene number, show all genes.}
}
\description{
Plot highest expressed genes
}
\details{
Each row represents a gene. In a row, each bar represents a sample,
the x axis tells the percentage of counts accounted for across the whole
dataset. (counts of that gene in that sample / total counts of all genes in
all samples)

Genes are sorted by average expression.
}
\examples{
as_SingleCellExperiment(sim_mat) \%>\% plot_highest_exprs()

}
