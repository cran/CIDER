% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNetwork.R
\name{plotNetwork}
\alias{plotNetwork}
\title{Plot Network Graph}
\usage{
plotNetwork(
  seu,
  ider,
  batch.var = "Batch",
  colour.by = NULL,
  weight.factor = 6.5,
  col.vector = NULL,
  vertex.size = 1
)
}
\arguments{
\item{seu}{Seurat S4 object after the step of \code{getIDER}, containing 
\code{initial_cluster} and \code{Batch} in its meta.data. Required.}

\item{ider}{A list. Output of `getIDER`. Required.}

\item{batch.var}{Character. Metadata colname containing batch information.
(Default: \code{Batch})}

\item{colour.by}{Character. It should be one of the colnames of Seurat 
object meta.data.It is used to colour the vertex of the network graph. 
(Default: \code{NULL})}

\item{weight.factor}{Numerical. Adjust the thickness of the edges. 
(Default: 6.5)}

\item{col.vector}{A vector of Hex colour codes. If no value is given 
(default), a vector of 74 colours will be used.}

\item{vertex.size}{Numerical. Adjsut the size of vertexes. (Default: 1)}
}
\value{
An igraph object
}
\description{
Visualise the network based on an IDER-based similarity matrix. 
The vertexes are initial clusters, and
the edge width denotes the similarity between two initial clusters.
}
\seealso{
\code{\link{getIDEr}}
}
