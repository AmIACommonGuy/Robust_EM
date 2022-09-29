library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)

# function: sim_mixture
# Description: function that generates simulation data with outliers.
# Input:  
#         cluster       =       number of clusters
#         d             =       The number of dimensions of the data
#         separation    =       minimum distance between centers of clusters
#         n             =       number of points in each cluster
#         out_perc      =       percentage of outliers in the dataset 
#         cov_scale     =       scaling coefficient for covariance
#         df            =       degree of freedom, Inf when normal
#         outl_mode     =       Distributuin of outliers. There are two 
#                               options -'cauchy' or 'uniform'
# Return:
#         A list of simulated observations.
#         x             =       all observations.
#         labelo        =       label with outliers (0)
#         labelin       =       label without outliers
#         gauss         =       non-outliers
#         outliers      =       outliers
#         mu            =       list of averaged cluster locations
#         sigma         =       list of cluster covariance
#
sim_mixture <- function(cluster, d, separation, n, out_perc, cov_scale = 1, df = Inf, outl_mode = 'cauchy') {
    ## Number of outliers determined here
    n_inlier <- round(n*(1-out_perc))
    n_outlier <- n - n_inlier
    ## Simulate mu
    mu <- runif(d, 0, 100)
    count = 0
    while (count < cluster) {
        new_mu <- runif(d, 0, 100)
        temp <- rbind(mu, new_mu)
        if (min(dist(temp)) < separation) {
            next
        }
        mu <- temp
        count = count + 1
    }
    ## Simulate sigma
    sigma <- c()
    for(i in c(1:cluster)) {
        sigma[[i]] = diag(runif(d, 4, 25)) * cov_scale
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
                mhd <- mahalanobis(out_coef, mu[i,], sigma[[i]])
                if (mhd > 10 & mhd < 45) {
                    outliers <- rbind(outliers,out_coef)
                    count = count + 1
                } 
            } 
        } else if (outl_mode == 'unif') {
            while(count!=n_outlier) {
                ## Draw outliers from multivariate Cauchy
                out_coef <- runif(d,-3.5,3.5)
                out_coef <- out_coef + sign(out_coef)* diag(sigma[[i]]) + mu[i,]
                outliers <- rbind(outliers,out_coef)
                count = count + 1
            } 
            
        } else {
            stop('invalid mode')
        }

    
    }
    labelo <- c(label, rep(0, dim(outliers)[1]))
    return(list(x = rbind(gauss, outliers), label=labelo, label_in = label, gauss = gauss, outlier=outliers, mu = mu, sigma = sigma))
}













