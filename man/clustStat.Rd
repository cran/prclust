\name{clustStat}
\alias{clustStat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
External Evaluation of Cluster Results
}
\description{
clusterStat returns Rand, adjusted Rand, Jaccard index, measuring the agreement between estimated cluster and the truth with a higher value indicating a higher agreement.
}
\usage{
clustStat(trueGroup, group)
}
\arguments{
  \item{trueGroup}{The true cluster results.}
  \item{group}{The estimated cluster results, not neccessary calculating by PRclust.}
}
\references{
Rand, William M. "Objective criteria for the evaluation of clustering methods." \emph{Journal of the American Statistical association} 66.336 (1971): 846-850.

Hubert, Lawrence, and Phipps Arabie. "Comparing partitions." \emph{Journal of classification} 2.1 (1985): 193-218.
}
\value{
The return value is a "clusterStat" class, providing the following information.
 \item{Rand}{Rand Index}
 \item{AdjustedRand}{Adjusted Rand Index}
 \item{Jaccard}{Jaccard Index}
}

\author{
Chong Wu
}
\examples{
	a <- rep(1:3,3)
	a
	b <- rep(c(4:6),3)
	b
	clustStat(a,b)
}

\keyword{clustStat}
