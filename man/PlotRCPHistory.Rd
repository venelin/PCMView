% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotRCPHistory.R
\name{PlotRCPHistory}
\alias{PlotRCPHistory}
\title{Plot the search path from a recursive clade partition search}
\usage{
PlotRCPHistory(fit, sizeGreyNodepoints = 2.2,
  sizeColorNodepoints = 2.2, sizeBlackAllowedModelTypes = 1.4,
  sizeColorAllowedModelTypes = 1.4, sizeRankInQueue = 1.4,
  vjustBlackAllowedModelTypes = -1.6,
  vjustColorAllowedModelTypes = -1.6, ...)
}
\arguments{
\item{fit}{an object of S3 class `PCMFitModelMappings` returned by a call
to `PCMFitMixed`.}

\item{sizeGreyNodepoints, sizeColorNodepoints, sizeBlackAllowedModelTypes, sizeColorAllowedModelTypes, sizeRankInQueue, vjustBlackAllowedModelTypes, vjustColorAllowedModelTypes}{graphical parameters (see function description).}

\item{...}{additional parameters passed to \code{\link{PCMTreePlot}}.}
}
\value{
a list of annotated ggtree plots. Some of the entries in this list
can be NULL to indicate that no score improvement has been achieved at the
corresponding iteration. The example below shows how to filter these out.
}
\description{
This function produces a list of tree-plots representing the
search-path during a recrusive clade partition search. Each plot has a label
on top showing:
\describe{
 \item{}{a number in parentheses (i) describes iteration i of the main loop.}
 \item{}{the score, the log-likelihood and the number of parameters of the model.}
 \item{}{a coloured node with a number i is the partition root for the iteration.}
 \item{}{Nodes in grey represent the potential shift points - these are
 descendants from the partition root, which have not been "cut out" by a shift
 and have at least q descendants, themselves.}
 \item{}{Letters in braces denote the candidate model-types for each shift-node.}
}
}
\examples{
lstPlots <- PlotRCPHistory(
  fitBMWithShiftsToMammalData, layout = "fan")
cowplot::plot_grid(plotlist = lstPlots[!sapply(lstPlots, is.null)])

}
