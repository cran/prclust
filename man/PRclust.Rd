\name{PRclust}
\alias{PRclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the Solution of Penalized Regression-Based Clustering.
}
\description{
	Clustering analysis is usually regarded as unsupervised learning in absence of a class label, as opposed to supervised learning. It can be performed through penalized regression with grouping pursuit. Prclust helps us peform penalized regression-based clustering with various loss functions and grouping penalities via two algorithm (DC-ADMM and quadratic penalty).
}
\usage{
PRclust(data, rho, lambda, tau, 
  loss.function = c("quadratic","L1"),
  grouping.penalty = c("gtlp","tlp"), 
  algorithm = c("DCADMM","Quadratic"),
  abs_res = 0.5,rel_res = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	 Input matrix. Each column is an observation vector.
}
   \item{rho}{
     Tuning parameter or step size: rho, typically set at 1 for quadratic penalty based algorithm; 0.4 for DC-ADMM. (Note that rho is the lambda1 in quadratic penalty based algorithm.)
  }
  \item{lambda}{
  	 Tuning parameter: lambda, the magnitude of grouping penalty.
}
  \item{tau}{
  	 Tuning parameter: tau, a nonnegative tuning parameter controll ing the trade-off between the model fit and the number of clusters.
}

\item{loss.function }{
      The loss function. "L1" stands for \eqn{L_1} loss function, while "quadratic" stands for the quadratic loss function.  
}
  \item{grouping.penalty}{
    Grouping penalty: "gtlp" is the group truncated lasso peanlty, while "tlp" is the truncated lasso. Note that when tau is large, "gtlp" and "tlp" are the $L_2$ penalty and lasso penalty, respectively.
}
  \item{algorithm}{
  Two algorithms for PRclust. "DC-ADMM" and "Quadratic" stand for the DC-ADMM and quadratic penalty based criterion respectively. "DC-ADMM" is much faster than "Quadratic" and thus recommend it here. 
}
\item{abs_res}{Primal residual, the stopping critetion parameter corresponding to DC-ADMM. The default is 0.5.
}
\item{rel_res}{Dual residual, the stopping critetion parameter corresponding to DC-ADMM. The default is 0.5.
}
}
\details{
Clustering analysis has been widely used in many fields. In the absence of a class label, clustering analysis is also called unsupervised learning. However, penalized regression-based clustering adopts a novel framework for clustering analysis by viewing it as a regression problem. In this method, a novel non-convex penalty for grouping pursuit was proposed which data-adaptively encourages the equality among some unknown subsets of parameter estimates. This new method can deal with some complex clustering situation, for example, in the presence of non-convex cluster, in which the K-means fails to work, PRclust might perform much better.
}
\value{
The return value is a list, which contains the following matrix.
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
Chong Wu
}

\note{
Choosing tunning parameter is kind of time consuming job. It is always based on "trials and errors".
}

\examples{
# Other examples are available upon request. 
## generate the data
data = matrix(NA,2,100)
data[1,1:50] = rnorm(50,0,0.33)
data[2,1:50] = rnorm(50,0,0.33)
data[1,51:100] = rnorm(50,1,0.33)
data[2,51:100] = rnorm(50,1,0.33)
## set the tunning parameter
a =PRclust(data,rho=1,lambda=3,tau=0.5)
a
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PRclust}
