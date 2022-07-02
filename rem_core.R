# Libraries
library(MASS)
library(Matrix)
library(mvtnorm)
library(car)
library(gmp)
library(mclust)
library(dplyr)

# function: EM_robust
# Input:  c           =       number of clusters
#         sampleMat   =       n*d matrix of data points (no cluster specified)
#         lambda      =       final lambda value - will go lambda^10 to lambda
#                             lambda will automatically equal 1, if not specified
#                             This is equivalent to a standard EM algorithm
#         d           =       The number of dimensions of the data
#         sigma_str   =       The covariance structure for each cluster - options
#                             are unstr, diag, or const
#         inits       =       Initial condition for EM. It contains cluster proportions,
#                             cluster average, and covariance matrixes for clusters. 
EM_robust = function(sampleMat, c, lambda = Inf, d, sigma_str, inits) {
    n=nrow(sampleMat) # n is the number of observations
    ## Observations
    x = sampleMat
    ## Initial values
    tau = matrix(inits[[3]], c, 1)
    mu = inits[[1]]
    sigma = inits[[2]]
    ## Posterior likelihood
    T_mat = matrix(0,c,n)
    max_it = 15
    if (lambda == Inf) {
        ll <- c(1)
    } else {
        ll <-unique(c(Inf, lambda^(5:1))) 
    }
    for (l in ll) {
        old_mu = mu
        for (it in 1:max_it) {
            old_T_mat = T_mat
            # E STEP: Construct the vector of the denom for T_mat for each i
            for (j in 1:c) {
                if (det(sigma[[j]]) < 1e-7) {
                    sigma[[j]] = diag(d)*1e-1
                }
                T_mat[j,] = log(tau[j]) + logdmvnorm(x,mu[j,],sigma = sigma[[j]])
            }
            for (i in 1:n) {
                scale = max(T_mat[,i])
                T_mat[,i] = exp(T_mat[,i] - scale - log(sum(exp(T_mat[,i] - scale))))
                ## if outlier detected, its better to assign it to the c
                if (any(is.na(T_mat[,i]))) {
                    T_mat[,i] = old_T_mat[,i]
                }
            }
            # Robust M STEP: Update tau, mu and sigma with l0 penalty
            for (j in 1:c) {
                e = matrix(0, n, d)
                
                # MHdist: Square of Mahalanobis distance
                x1 = x - rep(mu[j, ], each = n)
                sigma_inv = solve(as.matrix(sigma[[j]]))
                x1_sigma_inv = x1 %*% sigma_inv
                MHdist = matrix(rowSums(x1_sigma_inv * x1), n, 1)
                if (l == 1) {
                    e = matrix(0, n, d)
                } else {
                    for (i in 1:n) {
                        if (MHdist[i,] < l^2) {
                            e[i,] = 0 
                        }
                        else { e[i,] = (x[i,] - mu[j,]) }
                    }
                }

                # Update sigma
                num = matrix(0,d,d)
                denom = 0
                indices = which(e[,1]==0)
                # Update tau
                tau[j,] = (1/length(indices))*sum(T_mat[j,indices])
                # Update mu
                mu[j,] = colSums(T_mat[j,indices]%*%(x[indices,]-e[indices,]))/colSums(as.matrix(T_mat[j,indices]))
                if (length(indices) == 0) {
                    stop('lambda is too small, choose a bigger lambda to avoid no inliers case!')
                }
                if (sigma_str == "unstr") {
                    denom = sum(T_mat[j,indices])
                    if (l == 1) {
                        x_prime = x[indices,]-matrix(mu[j,], ncol=d, nrow=length(indices), byrow=T)
                    } else {
                        x_prime = x[indices,]-matrix(mu[j,], ncol=d, nrow=length(indices), byrow=T) 
                    }
                    num = t(x_prime)%*%diag(T_mat[j,indices], nrow=length(indices))%*%x_prime
                    sigma[[j]] = matrix(num/denom, d, d)
                }
                else if (sigma_str == "const") {
                    denom = sum(T_mat[j,indices])
                    if (l == 1) {
                        x_prime = x[indices,]-matrix(mu[j,], ncol=d, nrow=length(indices), byrow=T)
                    } else {
                        x_prime = x[indices,]-matrix(mu[j,], ncol=d, nrow=length(indices), byrow=T) 
                    }
                    num = sum(T_mat[j,indices]*x_prime^2)/d * diag(d)
                    denom = sum(T_mat[j,indices])
                    sigma[[j]] =  matrix(num/denom, d, d)
                }
                # Diagonal structure
                else {
                    for (i in 1:n){
                        num = num + T_mat[j,i]*(matrix(x[i,]-e[i,]-mu[j,])%*%t(matrix(x[i,]-e[i,]-mu[j,])))*diag(d)
                        denom = denom + T_mat[j,i]
                    }
                    sigma[[j]] = (num/denom)
                }
                # Eliminate uninformative cluster
                if (sum(T_mat[j,]!=0) <= 3) {
                    mu[j,] = matrix(0, 1, d)
                    sigma[[j]] = diag(d) * 1e-2
                }
            }
            if (any(is.na(max(abs(old_T_mat - T_mat))))) {stop('The lambda value is too small')}
            if (max(abs(old_T_mat - T_mat)) < 1e-7) {
                break
            }
        }
    }
    inlier_list <- c()
    ## Predict outliers
    for (j in 1:c) {
        e = matrix(0, n, d)
        # MHdist: Square of Mahalanobis distance
        x1 = x - rep(mu[j, ], each = n)
        sigma_inv = solve(as.matrix(sigma[[j]]))
        x1_sigma_inv = x1 %*% sigma_inv
        MHdist = matrix(rowSums(x1_sigma_inv * x1), n, 1)
        for (i in 1:n) {
            ## try out dif qchisq 099 is perferred
            if (MHdist[i,] < lam^2) {
                inlier_list <- c(inlier_list, i)
            }
        }
    }
    
    inlier_list <- unique(inlier_list)
    # Assign the column names to the specified cluster
    Soft_assign = t(T_mat)
    colnames(Soft_assign) = paste(1:c)
    hard_assign = matrix(paste(apply(Soft_assign, 1, which.max)), n, 1)
    hard_assign_out <- hard_assign
    hard_assign_out[-inlier_list, 1] <- 0
    # Construct return list
    returnList = list(mu, sigma, T_mat, tau, c, d, n, lambda, hard_assign, hard_assign_out, inlier_list)
    names(returnList) = c("mu", "sigma", "T_mat", "tau", "c", "d", "n", "lambda", "hard_assign","hard_assign_out", "inliers")
    return(returnList)
}
