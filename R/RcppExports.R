clustStat <- function(trueGroup, group) {
    x = as.vector(trueGroup)
    y = as.vector(group)
    if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
    tab <- table(x, y)
    if (all(dim(tab) == c(1, 1)))
    ARI <- 1
    a <- sum(choose(tab, 2))
    b <- sum(choose(rowSums(tab), 2)) - a
    c <- sum(choose(colSums(tab), 2)) - a
    d <- choose(sum(tab), 2) - a - b - c
    ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b +
    a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
    
    tempres <- .Call('prclust_clusterStat', PACKAGE = 'prclust', x, y)
    RAND <- (tempres$a + tempres$b) /(tempres$a + tempres$b + tempres$c + tempres$d)
    Jaccard <- (tempres$a) /(tempres$a + tempres$c + tempres$d)
    
    
    out = list()
    #out["TrueGroup"] <- trueGroup
    #out["EstimatedGroup"] <- group
    out["Rand"] <- RAND
    out["AdjustedRand"] <- ARI
    out["Jaccard"] <- Jaccard
    class(out) <- "clustStat"
    out
}

PRclust <- function(data, rho, lambda, tau, loss.function = c("quadratic","L1"),grouping.penalty = c("gtlp","tlp"), algorithm = c("DCADMM","Quadratic"),abs_res = 0.5,rel_res = 0.5) {
    
    ## judge for different situation
    mumethod = switch(match.arg(loss.function), `quadratic` = 0,L1 = 1)
    methods = switch(match.arg(grouping.penalty), `gtlp` = 0,tlp = 1)
    nalgorithm = switch(match.arg(algorithm), `DCADMM` = 1,Quadratic = 2)
    lambda1 = rho
    if(is.character(lambda1))
        stop("rho must be a postive number, you can use either GCV or stability based method to choose good tunning parameters.")
    if(is.character(lambda))
        stop("lambda must be a postive number,you can use either GCV or stability based method to choose good tunning parameters.")
    if(is.character(tau))
        stop("tau must be a postive number,you can use either GCV or stability based method to choose good tunning parameters.")
        
    if(lambda1<0 | is.na(lambda1))
        stop("rho must be a postive number, you can use GCV to choose the 'best' tunning parameter.")
    if(lambda<0 | is.na(lambda))
        stop("lambda must be a postive number, you can use GCV to choose the 'best' tunning parameter.")
    if(tau<0 | is.na(tau))
        stop("tau must be a postive number,you can use either GCV or stability based method to choose good tunning parameters.")
    
    data = as.matrix(data)
    if(sum(is.na(data)))
        stop("Clustering data contains NA or character value. The current version does not support missing data.")
    
    if( nalgorithm ==1){
        res = .Call('prclust_DCADMM', PACKAGE = 'prclust', data, rho, lambda, tau, abs_res , rel_res,mumethod, methods )
        final.count = sum(res$count2)
    } else {
        if (mumethod!= 0) {
            stop("Quadtraic penalty based algorithm cannot deal with the selected objective function. You can try ADMM instead.")
        }
        
        res = .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data, lambda1, lambda, tau, mumethod,methods)
        final.count = res$count
    }
    
    res = list(mu = res$mu,count = final.count,group = res$group,
    theta = res$theta,rho = lambda1, lambda = lambda,tau = tau, method = methods, algorithm = nalgorithm)
    class(res) = "prclust"
    res
}

print.clustStat <-function(x, ...) {
    cat("External evaluation of cluster results:\n")
    cat(paste("The Rand index: ",x$Rand,"\n",sep = ""))
    cat(paste("The adjusted rand index: ",x$AdjustedRand,"\n",sep = ""))
    cat(paste("The Jaccard index: ", x$Jaccard,"\n",sep = ""))
}

print.prclust <- function(x, ...) {
    temp.group = x$group
    max.groupnum = max(temp.group)
    cat(paste("Penalized regression-based clustering (prclust) with ",max.groupnum," clusters.\n",sep = ""))
    cat(paste("The iteration time is ",x$count,".\n",sep = ""))
    
    cat("\nThe centroids of observations:\n")
    print(x$mu)
    
    cat("\nClustering vector:\n")
    print(x$group)
    invisible(x)
}
