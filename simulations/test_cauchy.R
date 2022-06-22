library(MASS)
library(stats)
library(mvtnorm)
library(LaplacesDemon)

cluster = 2
lam = 2.5
n = 100
d = 2
maxl <- runif(1,7,14)
(mu = runif(d,1,100))
out_perc = 0.1
sigma = diag(runif(d, 1, maxl))
n_inlier <- round(n*(1-out_perc))
n_outlier <- n - n_inlier
gauss1 = mvrnorm(n_inlier, mu = mu, Sigma = sigma)
ggplot() + geom_point(aes(x = gauss1[,1], y = gauss1[,2]), col = 'darkgrey')+ geom_point(aes(x = mu[1], y =mu[2]), col = 'red')

outliers <- c()
while(dim(outliers)[1]!=n_outlier) {
    out_coef <- rmvc(1, mu, sigma)
    
    mhd <- (out_coef-mu)%*%solve(sigma)%*%t(out_coef-mu)
    if (mhd > 14 & mhd < 140) {
        outliers <- rbind(outliers,out_coef)
    } 
}

a <- multivarGaussian_cauchy(n, d, out_perc = 0.2, out_mag=1, independent = TRUE, cov_scale = 1)
a$outlier
ggplot() + geom_point(aes(x = a$gauss[,1], y = a$gauss[,2]), col = 'darkgrey')+ 
    geom_point(aes(x = a$mu[1], y =a$mu[2]), col = 'red')+
    geom_point(aes(x = a$outlier[,1], y =a$outlier[,2]), col = 'blue')
