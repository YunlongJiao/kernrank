% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kendall_top.R
\name{kendall_top}
\alias{kendall_top}
\title{Kendall kernel for top-k rankings}
\usage{
kendall_top(x, y)
}
\arguments{
\item{x, y}{Vector. 
If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
NAs replace unobserved values.}
}
\value{
Kendall kernel for top-k rankings,
where ties (supposed few) are broken by adopting a convolution kernel averaging compatible rankings without ties.
}
\description{
Calculates Kendall kernel between top-k rankings in time \code{O(klogk)},
where ties (supposed few) are broken by adopting a convolution kernel
averaging compatible rankings without ties.
}
\examples{
x <- c(1.5, 0.1, NA, -4, NA)
y <- c(NA, NA, 0, 3, NA)
kendall_top(x, y)

}
\references{
Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
}
\author{
Yunlong Jiao
}
\keyword{Kendall}
\keyword{Kernel}
\keyword{TopkRanking}
