library(MASS)
library(Matrix)
library(mvtnorm)
library(mixtools)
library(dplyr)
library(stats)

## This function is used to create the initial guess for the EM algorithm
initial_hier = function(sampleMat, cluster) {
    ## Same as mclust might stuck in the local maximum
    hclass1 = hclass(hc(sampleMat), G=cluster) 

    mat_withlabels <- as.data.frame(sampleMat) %>% mutate(class = hclass1)
    initial_means = mat_withlabels %>% group_by(class) %>%
        summarise_all(mean) %>% select(-class) %>% as.matrix()
    initial_cov = lapply(1:cluster, function(x) {sampleMat[which(hclass1 == x), ] %>% cov()})
    initial_tau = sapply(1:cluster, function(x) {sum(hclass1 == x) / length(hclass1)})
    return(list(initial_means, initial_cov, initial_tau))
}








