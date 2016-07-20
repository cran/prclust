\name{prclust-package}
\docType{package}
\title{
Penalized Regression Based Cluster Method
}
\description{
Clustering analysis is widely used in many fields. Traditionally clustering is regarded as unsupervised learning for its lack of a class label or a quantitative response variable, which in contrast is present in supervised learning such as classification and regression. Here we formulate clustering as penalized regression with grouping pursuit. In addition to the novel use of a non-convex group penalty and its associated unique operating characteristics in the proposed clustering method, a main advantage of this formulation is its allowing borrowing some well established results in classification and regression, such as model selection criteria to select the number of clusters, a difficult problem in clustering analysis. In particular, we propose using the generalized cross-validation (GCV) based on generalized degrees of freedom (GDF) to select the number of clusters. we further develop this method by developing a more efficient algorithm for scalable computation as well as a new theory for PRclust. This algorithm, called DC-ADMM, combines difference of convex programming with the alternating direction method of multipliers (ADMM). This method is more efficient than the quadratic penalty algorithm used in Pan et al. (2013) due to the availability of closed-form updating formulas.
}
\details{
\tabular{ll}{
Package: \tab prclust\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2016-7-18\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Chong Wu, Wei Pan\cr
Maintainer: Chong Wu <wuxx0845@umn.edu>
}
\references{
Pan Wei, Xiaotong Shen, and Binghui Liu. "Cluster Analysis: Unsupervised Learning via Supervised Learning with a Non-convex Penalty." \emph{The Journal of Machine Learning Research} 14.1 (2013):1865-1889.

Chong Wu, Sunghoon Kwon, Xiaotong Shen and Wei Pan. "A new Algorithm and Theory for Penalized Regression-based Clustering", submitted. 
}
\examples{
## In default, we use DC-ADMM, a faster algorithm to solve 
## the objective function and get the clustering result.

## generate the data
data = matrix(NA,2,100)
data[1,1:50] = rnorm(50,0,0.33)
data[2,1:50] = rnorm(50,0,0.33)
data[1,51:100] = rnorm(50,1,0.33)
data[2,51:100] = rnorm(50,1,0.33)

# clustering via PRclsut 
a =PRclust(data,rho=0.4,lambda=1,tau=0.5)
a$mu
a$group
}
