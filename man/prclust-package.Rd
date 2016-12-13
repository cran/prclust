\name{prclust-package}
\alias{prclust-package}
\alias{prclust}
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
Version: \tab 1.3\cr
Date: \tab 2016-12-12\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Chong Wu, Wei Pan\cr
Maintainer: Chong Wu <wuxx0845@umn.edu>
}
\references{
Pan, W., Shen, X., & Liu, B. (2013). Cluster analysis: unsupervised learning via supervised learning with a non-convex penalty. \emph{Journal of Machine Learning Research}, 14(1), 1865-1889.

Wu, C., Kwon, S., Shen, X., & Pan, W. (2016). A New Algorithm and Theory for Penalized Regression-based Clustering. \emph{Journal of Machine Learning Research}, 17(188), 1-25.
}
\keyword{ prclust-package }
\examples{
## In default, we use DC-ADMM, a faster algorithm to solve 
## the objective function and get the clustering result.
library("prclust")
## generate the data
data = matrix(NA,2,100)
data[1,1:50] = rnorm(50,0,0.33)
data[2,1:50] = rnorm(50,0,0.33)
data[1,51:100] = rnorm(50,1,0.33)
data[2,51:100] = rnorm(50,1,0.33)

# clustering via PRclsut 
a =PRclust(data,lambda1=0.4,lambda2=1,tau=0.5)
a$mu
a$group
}
