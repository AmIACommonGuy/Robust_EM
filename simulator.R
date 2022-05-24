library(MASS)
library(stats)
multivarGaussian_unif = function(n, d, out_perc, out_mag, independent = TRUE, cov_scale = 1){
    # Mean of the cluster
    maxl <- runif(1,3,4)
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
        gauss = mvrnorm(n, mu = mu, Sigma = sigma) #, tol = 1)
    }
    else {
        n_inlier <- round(n*(1-out_perc))
        n_outlier <- n - n_inlier
        
        # Else, pull from an MVN with the specified covariance matrix
        gauss1 = mvrnorm(n_inlier, mu = mu, Sigma = sigma) #, tol = 1)
        
        # Then pull the outliers from a uniform dist.
        outliers <- c()
        for (i in c(1:n_outlier)) {
            a <- runif(d,-3*msd,3*msd)
            outliers = rbind(outliers ,mu + (sign(a)*runif(1,-3*msd,3*msd) + a)*cov_scale * out_mag)
        }
        
        # Combine for one dataset
        gauss = rbind(gauss1, outliers)
    }
    
    mvGauss = list(mu = mu, sigma = sigma, gauss = gauss,outlier = outliers, n_inlier = n_inlier)
    return(mvGauss)
    
}

simMultGauss_unif = function(n, d, cluster, out_perc, out_mag, cov_scale = 1){
    samples_simMultGauss = replicate(cluster, multivarGaussian_unif(n = n, d = d,
                                                                    out_perc = out_perc, out_mag = out_mag, cov_scale))
    sample_inline = do.call(rbind, samples_simMultGauss[5,])
    loop = 0
    idx = c()
    for (i in c(1:cluster)){
        print(i)
        temp = loop*n+c(1:sample_inline[i,])
        idx=c(idx,temp)
        loop = loop+1
    }
    sampleMu = do.call(rbind, samples_simMultGauss[1,])
    sampleSigma = lapply(samples_simMultGauss[2,], function(y) as.matrix(y))
    simSamp = do.call(rbind, samples_simMultGauss[3,])
    simOut = do.call(rbind, samples_simMultGauss[4,])
    return(list(mus = sampleMu, sigmas = sampleSigma, simdata = simSamp, outliers= simOut, idx))
}



multivarGaussian_nm = function(n, d, out_perc, out_mag, independent = TRUE, cov_scale = 1){
    mu = runif(d,1,50)
    
    if (independent){
        sigma = diag(runif(d, 1, 3))
    } else{
        sigma = matrix(runif(d*d, -3, 3), d, d)
        sigma = (sigma + t(sigma))/2
        eigs = eigen(sigma)$values
        if (min(eigs) <= 0.5) {sigma = sigma - (min(eigs) - 0.5) * diag(d)}
    }
    
    sigma = cov_scale * sigma * runif(1,1,4)
    sigma_out = cov_scale * sigma * out_mag
    if (out_perc == 0 ) {
        gauss = mvrnorm(n, mu = mu, Sigma = sigma) #, tol = 1)
    }
    else {
        n_inlier <- round(n*(1-out_perc))
        n_outlier <- n - n_inlier
        
        # Else, pull from an MVN with the specified covariance matrix
        gauss1 = mvrnorm(n_inlier, mu = mu, Sigma = sigma) #, tol = 1)
        
        # Then pull the outliers from an MVN with the scaled up covariance matrix
        gauss2 = mvrnorm(n_outlier, mu = mu, Sigma = sigma_out) #, tol= 1)
        
        # Combine for one dataset
        gauss = rbind(gauss1, gauss2)
    }
    
    mvGauss = list(mu = mu, sigma = sigma, gauss = gauss)
    return(mvGauss)
}

simMultGauss_nm = function(n, d, cluster, out_perc, out_mag, cov_scale = 1){
    samples_simMultGauss = replicate(cluster, multivarGaussian_nm(n = n, d = d,
                                                                  out_perc = out_perc, out_mag = out_mag, cov_scale))
    
    sampleMu = do.call(rbind, samples_simMultGauss[1,])
    sampleSigma = lapply(samples_simMultGauss[2,], function(y) as.matrix(y))
    simSamp = do.call(rbind, samples_simMultGauss[3,])
    
    return(list(mus = sampleMu, sigmas = sampleSigma, simdata = simSamp))
}
