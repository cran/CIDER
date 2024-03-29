% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{hdbscan.seurat}
\alias{hdbscan.seurat}
\title{Initial clustering for evaluating integration}
\usage{
hdbscan.seurat(seu, reduction = "pca", dims = seq_len(15), minPts = 25)
}
\arguments{
\item{seu}{a Seurat object containing integrated or batch corrected
PCA.}

\item{reduction}{Character. Name of the dimension reduction after
integration or batch correction. (Default: PCA)}

\item{dims}{Numeric vector. Dimensions used for initial clustering.
(Default: 1:15)}

\item{minPts}{Interger. Minimum size of clusters. Will be passed
to the `hdbscan` function. (Default: 25)}
}
\value{
A Seurat object having two additional columns in its
meta.data: dbscan_cluster and initial_cluster.
}
\description{
This function applies HDBSCAN, a density-based
clustering method, on the corrected dimension reduction.
}
\seealso{
Usage of this function should be followed by
\code{\link{getIDEr}} and \code{\link{estimateProb}}.
}
