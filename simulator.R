library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)


sim_mixture <- function(cluster, d, seperation, n, out_perc, out_mag = 1, cov_scale = 1, df = Inf, outl_mode = 'cauchy') {
    ## Number of outliers determined here
    n_inlier <- round(n*(1-out_perc))
    n_outlier <- n - n_inlier
    ## Simulate mu
    mu <- runif(d, 0, 100)
    count = 0
    while (count < cluster) {
        new_mu <- runif(d, 0, 100)
        temp <- rbind(mu, new_mu)
        if (min(dist(temp)) < seperation) {
            next
        }
        mu <- temp
        count = count + 1
    }
    ## Simulate sigma
    sigma <- c()
    for(i in c(1:cluster)) {
        sigma[[i]] = diag(runif(d, 4, 25))*cov_scale
    }
    gauss <- c()
    label <- rep(c(1:cluster),each = n_inlier)
    for(i in c(1:cluster)) {
        gauss <- rbind(gauss, LaplacesDemon::rmvt(n_inlier, mu = mu[i,], S = sigma[[i]], df = df))
    }
    
    ## add outliers
    outliers <- c()
    for (i in c(1:cluster)) {
        count <- 0
        if (outl_mode == 'cauchy') {
            while(count!=n_outlier) {
                ## Draw outliers from multivariate Cauchy
                out_coef <- rmvc(1, mu[i,], sigma[[i]])
                mhd <- mahalanobis(out_coef, mu[i,], sigma[[i]]) * out_mag
                if (mhd > (5-out_mag) & mhd < (45+out_mag)) {
                    outliers <- rbind(outliers,out_coef)
                    count = count + 1
                } 
            } 
        } else if (outl_mode == 'unif') {
            while(count!=n_outlier) {
                ## Draw outliers from multivariate Cauchy
                out_coef <- runif(d,-3.5,3.5)
                out_coef <- out_coef + sign(out_coef)* diag(sigma[[i]]) * 2
                mhd <- mahalanobis(out_coef, mu[i,], sigma[[i]]) * out_mag
                if (mhd > (5-out_mag) & mhd < (45+out_mag)) {
                    outliers <- rbind(outliers,out_coef)
                    count = count + 1
                } 
            } 
            
        } else {
            stop('invalid mode')
        }

    
    }
    labelo <- c(label, rep(0, dim(outliers)[1]))
    return(list(x = rbind(gauss, outliers), label=labelo, label_in = label, gauss = gauss, outlier=outliers, mu = mu, sigma = sigma))
}













