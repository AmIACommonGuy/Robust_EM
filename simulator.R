library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)
multivarGaussian_unif = function(n, d, out_perc, out_mag, independent = TRUE, cov_scale = 1){
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
        for (i in c(1:n_outlier)) {

            out_coef <- runif(d,-3.5,3.5)
            while (sum(out_coef^2) < 2.8) {
                out_coef <- runif(d,-3.5,3.5)
            }
            a <- out_coef*out_mag*sqrt(diag(sigma)) + mu
            outliers = rbind(outliers , a)
        }
        
        # Combine for one dataset
        gauss = rbind(gauss1, outliers)
    }
    
    mvGauss = list(mu = mu, sigma = sigma, gauss = gauss,outlier = outliers, n_inlier = n_inlier)
    return(mvGauss)
    
}

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
        while(length(outliers)!=n_outlier) {
            out_coef <- rmvc(1, mu, sigma)
            
            mhd <- (out_coef-mu)%*%solve(sigma)%*%t(out_coef-mu)
            if (mhd > 21 & mhd < 300) {
                outliers <- rbind(outliers,out_coef)
            } 
        }
        
        # Combine for one dataset
        gauss = rbind(gauss1, outliers)
    }
    
    mvGauss = list(mu = mu, sigma = sigma, gauss = gauss, outlier = outliers, n_inlier = n_inlier)
    return(mvGauss)
    
}




simMultGauss_unif = function(n, d, cluster, out_perc, out_mag, cov_scale = 1){
    samples_simMultGauss = replicate(cluster, multivarGaussian_unif(n = n, d = d,
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
    while (min(dist(sampleMu)) < 21 | max(dist(sampleMu)) > 63) {
        samples_simMultGauss = c()
        samples_simMultGauss = replicate(cluster, 
                                         multivarGaussian_unif(n = n, d = d,
                                                               out_perc = out_perc, out_mag = out_mag, cov_scale = 1))
        sampleMu = do.call(rbind, samples_simMultGauss[1,])
        sampleSigma = lapply(samples_simMultGauss[2,], function(y) as.matrix(y))
        }
    simSamp = do.call(rbind, samples_simMultGauss[3,])
    simIn = simSamp[idx,]
    simOut = simSamp[-idx,]
    Label <- rep(1:cluster, each=n)
    Labelo <- Label
    Labelo[-idx] <- cluster+1
    return(list(mus = sampleMu, sigmas = sampleSigma, simdata = simSamp, outliers= simOut, inliers = simIn, Label = Label, LabelO = Labelo))
}

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















