\name{stability}
\alias{stability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the stability based statistics}
\description{
Calculate the  the stability based statistics. We try with various tuning parameter values, obtaining their corresponding statbility based statistics average prediction strengths, then choose the set of the tuning parameters with the maximum average prediction stength.
}
\usage{
stability(data,rho,lambda,tau,
    loss.function = c("quadratic","L1","MCP","SCAD"),
    grouping.penalty = c("gtlp","tlp"), 
    algorithm = c("DCADMM","Quadratic"),
    epsilon = 0.001,n.times = 10)
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
  \item{grouping.penalty}{Grouping penalty. Character: may be abbreviated. "gtlp" means generalized group lasso is used for grouping penalty. "lasso" means lasso is used for grouping penalty. "SCAD" and "MCP" are two other non-convex penalty.
}
  \item{algorithm}{
  Two algorithms for PRclust. "DC-ADMM" and "Quadratic" stand for the DC-ADMM and quadratic penalty based criterion respectively. "DC-ADMM" is much faster than "Quadratic" and thus recommend it here. 
}
\item{epsilon}{The stopping critetion parameter corresponding to DC-ADMM. The default is 0.001.
}
\item{n.times}{Repeat times. Based on our limited simulations, we find 10 is usually good enough.}
}

\details{
A generalized degrees of freedom (GDF) together with generalized cross validation (GCV) was proposed for selection of tuning parameters for clustering (Pan et al., 2013). This method, while yielding good performance, requires extensive computation and specification of a hyper-parameter perturbation size. Here, we provide an alternative by modifying a stability-based criterion (Tibshirani and Walther, 2005; Liu et al., 2016) for determining the tuning parameters.

The main idea of the method is based on cross-validation. That is, (1) randomly partition the entire data set into a training set and a test set with an almost equal size; (2) cluster the training and test sets separately via PRclust with the same tuning parameters; (3) measure how well the training set clusters predict the test clusters.

We try with various tuning parameter values, obtaining their corresponding statbility based statistics average prediction strengths, then choose the set of the tuning parameters with the maximum average prediction stength.
}
\value{
Return value: the average prediction score.
}

\references{
Wu, C., Kwon, S., Shen, X., & Pan, W. (2016). A New Algorithm and Theory for Penalized Regression-based Clustering. \emph{Journal of Machine Learning Research}, 17(188), 1-25.
}
\author{
Chong Wu
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
stab1 = stability(data,rho=1,lambda=1,tau=0.5,n.times = 2)
stab1

#case 2
stab2 = stability(data,rho=1,lambda=0.7,tau=0.3,n.times = 2)
stab2
# Note that the combination of tuning parameters in case 1 are better than 
# the combination of tuning parameters in case 2 since the value of GCV in case 1 is
# less than the value in case 2.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{stability}
