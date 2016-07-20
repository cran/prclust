###################################
## general degree of freedom and GCV
###################################

GDF.Step23 <- function(seed,data,lambda1,lambda2,tau,mumethods, methods,sigma,algorithm,abs_res,rel_res)
{
    set.seed(seed) ## Set the random seed for this simulation
    # peturbation of data
    deltaB = matrix(rnorm(dim(data)[1]*dim(data)[2],0,sigma),dim(data)[1],dim(data)[2])
    data1 = data + deltaB
    if (algorithm == 1){
        a = .Call('prclust_DCADMM', PACKAGE = 'prclust', data1, lambda1, lambda2, tau, abs_res, rel_res, mumethods, methods)
    } else {
        a =    .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data1, lambda1, lambda2, tau,mumethods, methods)

    }
    out = list()
    out[[1]] = deltaB
    out[[2]] = a$mu
    out[[3]] = a$group
    out
}

##########################################

##########################################

GCV <- function(data,rho,lambda,tau,sigma,B=100,loss.function = c("quadratic","L1"),grouping.penalty = c("gtlp","tlp"), algorithm = c("DCADMM","Quadratic"),abs_res = 0.5,rel_res = 0.5)
{
    ## judge for different situation
    mumethods = switch(match.arg(loss.function), `quadratic` = 0,L1 = 1)
    methods = switch(match.arg(grouping.penalty), `gtlp` = 0,tlp = 1)
    nalgorithm = switch(match.arg(algorithm), `DCADMM` = 1,Quadratic = 2)
    lambda1 = rho
    lambda2 = lambda
    if(is.character(lambda1))
    stop("lambda1 must be a number")
    if(is.character(sigma))
    stop("sigma must be a postive number")
    if(is.character(B))
    stop("B must be a postive number")
    if(is.character(lambda2))
    stop("lambda must be a number")
    if(is.character(tau))
    stop("tau must be a number")
    
    if(lambda1<0 | is.na(lambda1))
    stop("rho must be a postive number.")
    if(lambda<0 | is.na(lambda2))
    stop("lambda must be a postive number.")
    if(tau<0 | is.na(tau))
    stop("tau must be a postive number.")
    if(sigma<0 | is.na(sigma))
    stop("sigma must be a postive number.")
    if(B<0 | is.na(B))
    stop("B must be a postive integer.")
    
    B = as.integer(B)
    data = as.matrix(data)
    if(sum(is.na(data)))
    stop("Clustering data contains NA or character value. The current version does not support missing data.")
    
    ##require("multicore")
    if( nalgorithm ==2) {
        if (mumethods!= 0) {
            stop("Quadtraic penalty based algorithm cannot deal with the selected objective function. You can try ADMM instead.")
        }
    }

    res = mclapply(1:B,GDF.Step23,data = data,lambda1 = lambda1,lambda2 = lambda2,
    tau = tau, mumethods = mumethods, methods = methods,sigma = sigma,algorithm = nalgorithm,abs_res = abs_res,rel_res = rel_res)
    nrows = dim(data)[1]
    ncols = dim(data)[2]
    num = nrows * ncols
    slope = matrix(NA,1,num)
    for(ii in 1:num)
    {
        x = matrix(NA,B,1)
        y = matrix(NA,B,1)
        row = floor((ii-1)/ncols)+1
        col = ii %% ncols
        if (col == 0)
        {
            col = ncols
        }
        for(i in 1:B){
            
            x[i] = res[[i]][[1]][row,col]
            y[i] = res[[i]][[2]][row,col]
        }
        slope[ii] = coefficients(lm(y~x))[2]
    }
    
    GDF = sum(slope)
    ## calculate the original mu
    if (nalgorithm == 1){
        rho = lambda1
        a = .Call('prclust_DCADMM', PACKAGE = 'prclust', data, rho, lambda2, tau, abs_res, rel_res, mumethods, methods)
  } else {
        a =    .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data, lambda1, lambda2, tau, mumethods, methods)
        
    }
    
    GCV = sum((data-a$mu)^2)/(nrows*ncols- GDF)^2
    groupNum = length(unique(a$group))
    
    ## esitmated sigma
    estSigma = sum((data-a$mu)^2)/(nrows*ncols- GDF)
    #  sum((data[1:2,]-a$mu[1:2,])^2)/(2 * ncols- GDF)^2
    out = t(as.matrix(c(GDF,GCV,groupNum,estSigma)))
    colnames(out) = c("GDF","GCV","groupNum","estSigmaSquare")
    out
}
