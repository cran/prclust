stability<- function(data, rho,lambda, tau,loss.function = c("quadratic","L1","MCP","SCAD"),grouping.penalty = c("gtlp","tlp"), algorithm = c("DCADMM","Quadratic"), epsilon = 0.001,n.times = 10) {
    
    
    
    ## judge for different situation
    mumethod = switch(match.arg(loss.function), `quadratic` = 0,L1 = 1)
    methods = switch(match.arg(grouping.penalty), `gtlp` = 0,L1 = 1, MCP = 2, SCAD = 3)
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
    
    
    
    #n.times = 100
    final.res = matrix(NA,n.times,1)
    for(iteration in 1:n.times) {
        n= dim(data)[2]
        index = sample(1:n,n)
        index.tr = index[1:(n/2)]
        index.test = index[(n/2 + 1):n]
        data.tr = data[,index.tr]
        data.test = data[,index.test]
        
        data.tmp = t(data[,index])
        dist.tmp = as.matrix(dist(data.tmp))
        
        dist.tmp = dist.tmp[(n/2 + 1):n,]
        dist.tmp = dist.tmp[,1:(n/2)]
        
        testtr.index = sapply(seq(nrow(dist.tmp)), function(i) {
            j <- order(dist.tmp[i,])[2]
        })
        
        if( nalgorithm ==1){
            a = .Call('prclust_PRclustADMM', PACKAGE = 'prclust', data.tr, rho, lambda, tau,mumethod, methods,epsilon)
            #a = .Call('prclust_DCADMM', PACKAGE = 'prclust', data.tr, rho, lambda, tau, abs_res , rel_res,mumethod, methods )
            tr.res  = as.matrix(a$group)
            
            tmp = tr.res[testtr.index,]
            a = .Call('prclust_PRclustADMM', PACKAGE = 'prclust', data.test, rho, lambda, tau,mumethod, methods,epsilon)
            #a = .Call('prclust_DCADMM', PACKAGE = 'prclust', data.test, rho, lambda, tau, abs_res , rel_res,mumethod, methods )
        } else {
            if (mumethod!= 0) {
                stop("Quadtraic penalty based algorithm cannot deal with the selected objective function. You can try ADMM instead.")
            }
            a = .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data.tr, lambda1, lambda, tau, mumethod,methods)
            
            tr.res  = as.matrix(a$group)
            tmp = tr.res[testtr.index,]
            a = .Call('prclust_PRclustOriginal', PACKAGE = 'prclust', data.test, lambda1, lambda, tau, mumethod,methods)
            
        }
        
        # calculate aRand
        x = as.vector(tmp)
        y = as.vector(a$group)
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
    
        final.res[iteration,1] = ARI
    }
    
    mean(final.res)
}

