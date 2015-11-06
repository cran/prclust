\name{PRclust}
\alias{PRclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the Solution of Penalized Regression-Based Clustering.
}
\description{
	Clustering is unsupervised and exploratory in nature. Yet, it can be performed through penalized regression with grouping pursuit. Prclust helps us peform penalized regression-based clustering with various loss functions and grouping penalities via two algorithm (DC-ADMM and quadratic penalty).
}
\usage{
PRclust(data, lambda1, lambda2, tau, 
	loss.method = c("quadratic","lasso"), 
	group.method = c("gtlp","lasso","SCAD","MCP"), 
	algorithm = c("ADMM","Quadratic"), epsilon=0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	 input matrix, of dimension nvars x nobs; each column is an observation vector.
}
   \item{lambda1}{
     Tuning parameter or step size: lambda1, typically set at 1 for quadratic penalty based algorithm; 0.4 for revised ADMM.
  }
  \item{lambda2}{
  	 Tuning parameter: lambda2, the magnitude of grouping penalty.
}
  \item{tau}{
  	 Tuning parameter: tau, related to grouping penalty.
}

\item{loss.method }{
      The loss method. "lasso" stands for \eqn{L_1} loss function, while "quadratic" stands for the quadratic loss function.  
}
  \item{group.method}{
    character: may be abbreviated. "gtlp" means generalized group lasso is used for grouping penalty. "lasso" means lasso is used for grouping penalty. "SCAD" and "MCP" are two other non-convex penalty.
}
  \item{algorithm}{
	character: may be abbreviated. The algorithm to use for finding the solution. The default algorithm is "ADMM", which stands for the new algorithm we developed.
}
\item{epsilon}{The stopping critetion parameter. The default is 0.001.
}
}
\details{
Clustering analysis has been widely used in many fields. In the absence of a class label, clustering analysis is also called unsupervised learning. However, penalized regression-based clustering adopts a novel framework for clustering analysis by viewing it as a regression problem. In this method, a novel non-convex penalty for grouping pursuit was proposed which data-adaptively encourages the equality among some unknown subsets of parameter estimates. This new method can deal with some complex clustering situation, for example, in the presence of non-convex cluster, in which the K-means fails to work, PRclust might perform much better.
}
\value{
The return value is a list. In this list, it contains the following matrix.
 \item{mu}{The centroid of the each observations.}
 \item{theta}{The theta value for the data set, not very useful.}
 \item{group}{The group for each points.}
 \item{count}{The iteration times.}
}
\references{
Pan Wei, Xiaotong Shen, and Binghui Liu. "Cluster Analysis: Unsupervised Learning via Supervised Learning with a Non-convex Penalty." \emph{The Journal of Machine Learning Research} 14.1 (2013):1865-1889.

Chong Wu, Sunghoon Kwon, Xiaotong Shen and Wei Pan. "A new Algorithm and Theory for Penalized Regression-based Clustering", submitted. 
}
\author{
Chong Wu, Wei Pan
}
\note{
Choosing tunning parameter is kind of time consuming job. It is always based on "trials and errors".
}


\examples{
library("prclust")
# To let you have a better understanding about the power and strength
# of PRclust method, 6 examples in original prclust paper were provided.
################################################
### case 1
################################################
## generate the data
data = matrix(NA,2,100)
data[1,1:50] = rnorm(50,0,0.33)
data[2,1:50] = rnorm(50,0,0.33)
data[1,51:100] = rnorm(50,1,0.33)
data[2,51:100] = rnorm(50,1,0.33)
## set the tunning parameter
lambda1 =1
lambda2 = 3
tau = 0.5
a =PRclust(data,lambda1,lambda2,tau)
a

## quadratic penalty
lambda1 =1
lambda2 = 1
tau = 0.5
a =PRclust(data ,lambda1,lambda2,tau, algorithm ="Quadratic")
a

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PRclust}
