% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDistMat.R
\name{plotHeatmap}
\alias{plotHeatmap}
\title{Plot Heatmap for the IDER-Based Similarity Matrix}
\usage{
plotHeatmap(seu, ider, batch.var = "Batch")
}
\arguments{
\item{seu}{A Seurat object.}

\item{ider}{The output list from the \code{getIDEr} function.}

\item{batch.var}{Character string specifying the metadata column that contains batch information. Default is "Batch".}
}
\value{
A heatmap displaying the similarity between shared groups across batches.
}
\description{
This function generates a heatmap that visualises the similarity between shared 
groups across batches, as computed by \code{getIDEr}.
}
\seealso{
\code{\link{getIDEr}}
}
