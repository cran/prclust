# prclust
Penalzied regression based clustering 

## Install the package
We test it on R 3.2.1 in Linux server and 3.2.2 in Windows and Mac. *Note that for windows user, we need install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. For Mac user, we need install [gfortran](https://cran.r-project.org/bin/macosx/tools/).*

```
library(devtools)
install_github("ChongWu-Biostat/prclust") # install the prclust packages
```

## PRclust
Clustering is unsupervised and exploratory in nature. Yet, it can be performed through penalized regression with grouping pursuit
. Prclust helps us peform penalized regression-based clustering with various loss functions and grouping penalities via two algorithm (DC-ADMM and quadratic penalty).

```
library("prclust")
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
# More examples can be found in the manual
```

## Generalized Cross Validation

A bonus with the regression approach to clustering is the potential application of many existing model selection methods for regression or supervised learning to clustering. We propose using generalized cross-validation (GCV). GCV can be regarded as an approximation to leave-one-out cross-validation (CV). Hence, GCV provides an approximately unbiased estimate of the prediction error.

We use the generalized degrees of freedom (GDF) to consider the data-adaptive nature in estimating the centroids of the observations.

The chosen tuning parameters are the one giving the smallest GCV error.

```
library("prclust")
set.seed(1)
data = matrix(NA,2,50)
data[1,1:25] = rnorm(25,0,0.33)
data[2,1:25] = rnorm(25,0,0.33)
data[1,26:50] = rnorm(25,1,0.33)
data[2,26:50] = rnorm(25,1,0.33)

#case 1
gcv1 = GCV(data,lambda1=1,lambda2=1,tau=0.5,sigma=0.25)
gcv1

#case 2
gcv2 = GCV(data,lambda1=1,lambda2=0.7,tau=0.3,sigma=0.25)
gcv2

# Note that the combination of tuning parameters in case 1 are better than 
# the combination of tuning parameters in case 2 since the value of GCV in case 1 is
# less than the value in case 2.
```

## Evaluating the clustering results
Suppose we know the true cluster results beforehand. clusterStat provides Rand, adjusted Rand, Jaccard index to measure the quality of a cluster results.

```
a <- rep(1:3,3)
a
b <- rep(c(4:6),3)
b
clusterStat(a,b)
```

## Manual
If you like our pakcage *prclust*, please give us a star. You can download the *prclust* package manual [here](https://cutpi.com/upimages/1446325172.pdf). 




