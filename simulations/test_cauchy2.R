library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)
library(ggplot2)

multivarGaussian_cauchy = function(n, d, out_perc, out_mag, independent = TRUE, cov_scale = 1){
    # Mean of the cluster
    maxl <- runif(1,7,14)
    msd <- sqrt(maxl)
    mu = runif(d,1,100)
    if (independent){
        ## Diagonal
        sigma = diag(runif(d, 1, maxl))
    } else{
        sigma = matrix(runif(d*d, -maxl/2, maxl/2), d, d)
        sigma = (sigma + t(sigma))/2
        eigs = eigen(sigma)$values
        if (min(eigs) <= 0.5) {sigma = sigma - (min(eigs) - 0.5) * diag(d)}
    }
    sigma = cov_scale * sigma
    if (out_perc == 0 ) {
        gauss = mvrnorm(n, mu = mu, Sigma = sigma)
    }
    else {
        n_inlier <- round(n*(1-out_perc))
        n_outlier <- n - n_inlier
        
        # Else, pull from an MVN with the specified covariance matrix
        gauss1 = mvrnorm(n_inlier, mu = mu, Sigma = sigma)
        
        # Then pull the outliers from a uniform dist.
        outliers <- c()
        count <- 0
        while(count!=n_outlier) {
            out_coef <- rmvc(1, mu, sigma)
            mhd <- mahalanobis(out_coef, mu, sigma)
            if (mhd > 21 & mhd < 300) {
                outliers <- rbind(outliers,out_coef)
                count = count + 1
            } 
        }
        
        # Combine for one dataset
        gauss = rbind(gauss1, outliers)
    }
    
    mvGauss = list(mu = mu, sigma = sigma, gauss = gauss, outlier = outliers, n_inlier = n_inlier)
    return(mvGauss)
    
}

n = 100
d = 2
out_perc = 0.1
out_mag = 1

a <- multivarGaussian_cauchy(n, d, out_perc, out_mag, independent = TRUE, cov_scale = 1)
ggplot() + geom_point(aes(x = a$gauss[,1], y = a$gauss[,2]), col = 'darkgrey') + 
           geom_point(aes(x = a$outlier[,1], y = a$outlier[,2]), pch = 1, cex = 2)


simMultGauss_cauchy = function(n, d, cluster, out_perc, out_mag, cov_scale = 1){
    samples_simMultGauss = replicate(cluster, multivarGaussian_cauchy(n = n, d = d,
                                                                      out_perc = out_perc, out_mag = out_mag, cov_scale = 1))
    
    sample_inline = do.call(rbind, samples_simMultGauss[5,])
    loop = 0
    idx = c()
    for (i in c(1:cluster)){
        temp = loop*n+c(1:sample_inline[i,])
        idx=c(idx,temp)
        loop = loop+1
    }
    sampleMu = do.call(rbind, samples_simMultGauss[1,])
    sampleSigma = lapply(samples_simMultGauss[2,], function(y) as.matrix(y))
    simSamp = do.call(rbind, samples_simMultGauss[3,])
    simIn = simSamp[idx,]
    simOut = simSamp[-idx,]
    Label <- rep(1:cluster, each=n)
    Labelo <- Label
    Labelo[-idx] <- cluster+1
    return(list(mus = sampleMu, sigmas = sampleSigma, simdata = simSamp, outliers= simOut, inliers = simIn, Label = Label, LabelO = Labelo))
}
samples_simMultGauss = replicate(cluster, multivarGaussian_cauchy(n = n, d = d,
                                                                  out_perc = out_perc, 
                                                                  out_mag = out_mag, 
                                                                  cov_scale = 1))


cluster = 7
d = 2
seperation = 29
n = 50
out_perc <- 0.1

sim_cauchy <- function(cluster, d, seperation, n, out_perc) {
    n_inlier <- round(n*(1-out_perc))
    n_outlier <- n - n_inlier
    ## Simulate mu that are 
    mu <- runif(d, 0, 100)
    count = 0
    while (count < cluster) {
        new_toss <- runif(2,1,100)
        temp <- rbind(mu, new_toss)
        if (min(dist(temp)) < seperation) {
            print('re-toss')
            next
        }
        mu <- temp
        count = count + 1
    }
    sigma <- c()
    for(i in c(1:cluster)) {
        sigma[[i]] = diag(runif(d, 7, 16))
    }
    gauss <- c()
    label <- rep(c(1:cluster),each = n_inlier)
    for(i in c(1:cluster)) {
        gauss <- rbind(gauss, LaplacesDemon::rmvt(n_inlier, mu = mu[i,], S = sigma[[i]], df = 7))
        # gauss <- rbind(gauss, mvrnorm(40, mu = mu[i,], Sigma = sigma[[i]]))
    }
    
    
    ## add outliers
    outliers <- c()
    for (i in c(1:cluster)) {
        count <- 0
        while(count!=n_outlier) {
            out_coef <- rmvc(1, mu[i,], sigma[[i]])
            mhd <- mahalanobis(out_coef, mu[i,], sigma[[i]])
            if (mhd > 16 & mhd < 70) {
                outliers <- rbind(outliers,out_coef)
                count = count + 1
            } 
        } 
    }
    ggplot() + geom_point(aes(x = gauss[,1], y = gauss[,2]), col = label) +geom_point(aes(x = outliers[,1], y = outliers[,2]), pch = 1)
}

sim_cauchy(cluster, d, seperation, n, out_perc)
