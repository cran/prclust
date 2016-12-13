\name{GCV}
\alias{GCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the Generalized Cross-Validation Statistic (GCV)}
\description{
Calculate the generalized cross-validation statistic with generalized degrees of freedom.
}
\usage{
GCV(data,lambda1,lambda2,tau,sigma,B=100,
	loss.method = c("quadratic","lasso"),
	grouping.penalty = c("gtlp","L1","SCAD","MCP"), 
	algorithm = c("ADMM","Quadratic"), epsilon =0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{
	 Numeric data matrix .
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
  \item{sigma}{
  	 The perturbation size.
  }
  \item{B}{
  	 The Monte Carlo time. The defualt value is 100.
}
\item{loss.method }{
     character may be abbreviated. "lasso" stands for \eqn{L_1} loss function, while "quadratic" stands for the quadratic loss function.  
}
  \item{grouping.penalty}{
  	character: may be abbreviated. "gtlp" means generalized group lasso is used for grouping penalty. "lasso" means lasso is used for grouping penalty. "SCAD" and "MCP" are two other non-convex penalty.
}
  \item{algorithm}{
	character: may be abbreviated. The algorithm will use for finding the solution. The default algorithm is "ADMM", which stands for the DC-ADMM.
}
\item{epsilon}{The stopping critetion parameter. The default is 0.001.
}
}
\details{
A bonus with the regression approach to clustering is the potential application of many existing model selection methods for regression or supervised learning to clustering. We propose using generalized cross-validation (GCV). GCV can be regarded as an approximation to leave-one-out cross-validation (CV). Hence, GCV provides an approximately unbiased estimate of the prediction error.

We use the generalized degrees of freedom (GDF) to consider the data-adaptive nature in estimating the centroids of the observations.

The chosen tuning parameters are the one giving the smallest GCV error.
}
\value{
Return value: the Generalized cross-validation statistic (GCV)
}
\references{
Pan, W., Shen, X., & Liu, B. (2013). Cluster analysis: unsupervised learning via supervised learning with a non-convex penalty. \emph{Journal of Machine Learning Research}, 14(1), 1865-1889.
}
\author{
Chong Wu, Wei Pan
}
\examples{
set.seed(1)
library("prclust")
data = matrix(NA,2,50)
data[1,1:25] = rnorm(25,0,0.33)
data[2,1:25] = rnorm(25,0,0.33)
data[1,26:50] = rnorm(25,1,0.33)
data[2,26:50] = rnorm(25,1,0.33)

#case 1
gcv1 = GCV(data,lambda1=1,lambda2=1,tau=0.5,sigma=0.25,B =10)
gcv1

#case 2
gcv2 = GCV(data,lambda1=1,lambda2=0.7,tau=0.3,sigma=0.25,B = 10)
gcv2

# Note that the combination of tuning parameters in case 1 are better than 
# the combination of tuning parameters in case 2 since the value of GCV in case 1 is
# less than the value in case 2.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GCV}
