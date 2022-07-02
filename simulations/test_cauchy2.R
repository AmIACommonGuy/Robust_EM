library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)
library(ggplot2)

sim_cauchy <- function(cluster, d, seperation, n, out_perc, cov_scale, df) {
    n_inlier <- round(n*(1-out_perc))
    n_outlier <- n - n_inlier
    ## Simulate mu that are 
    mu <- runif(d, 0, 100)
    count = 0
    while (count < cluster) {
        new_toss <- runif(2,1,100)
        temp <- rbind(mu, new_toss)
        if (min(dist(temp)) < seperation) {
            next
        }
        mu <- temp
        count = count + 1
    }
    sigma <- c()
    for(i in c(1:cluster)) {
        sigma[[i]] = diag(runif(d, 10, 25))*cov_scale
    }
    gauss <- c()
    label <- rep(c(1:cluster),each = n_inlier)
    for(i in c(1:cluster)) {
        gauss <- rbind(gauss, LaplacesDemon::rmvt(n_inlier, mu = mu[i,], S = sigma[[i]], df = df))
        # gauss <- rbind(gauss, mvrnorm(40, mu = mu[i,], Sigma = sigma[[i]]))
    }
    
    
    ## add outliers
    outliers <- c()
    for (i in c(1:cluster)) {
        count <- 0
        while(count!=n_outlier) {
            out_coef <- rmvc(1, mu[i,], sigma[[i]])
            mhd <- mahalanobis(out_coef, mu[i,], sigma[[i]])
            if (mhd > 7 & mhd < 40) {
                outliers <- rbind(outliers,out_coef)
                count = count + 1
            } 
        } 
    }
    
    labelo <- c(label, rep(0, dim(outliers)[1]))
    return(list(x = rbind(gauss, outliers), label=labelo, label_in = label, gauss= gauss, outlier=outliers, mu = mu, sigma = sigma))
}

# a <- sim_cauchy(cluster, d, seperation, n, out_perc, cov_scale = cov_scale)
# ggplot() + geom_point(aes(x = a$x[,1], y = a$x[,2]), col = a$label+1)
# 
# NoiseInit_by_guess = sample(c(TRUE,FALSE), size = n*cluster, 
#                             replace = TRUE, prob = c(1,5)/6)
# result_mclust = Mclust(a$x, verbose=F, G=cluster, modelNames = "VVV",
#                        initialization = list(noise = NoiseInit_by_guess),
#                        control = emControl(tol=c(1.e-6, 1.e-7), itmax=15),)
# plot_mc(result_mclust)
# rand.index(result_mclust$classification, a$label)
# adjustedRandIndex(result_mclust$classification, a$label)
# result_mclust$parameters$mean
# result_mclust$parameters$variance
# result_mclust$classification
# classError(result_mclust$classification, a$label)

