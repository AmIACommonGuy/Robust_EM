library(MASS)
library(Matrix)
library(mvtnorm)
library(mixtools)
library(dplyr)
library(stats)

# function: initial_hier
# Description:        This function generate the initial condition needed for the robust EM method.
# Input:
#         sampleMat   =       n*d matrix of data points. (no cluster specified, n stands for number of points,
#                             d stands for dimension)
#         cluster     =       number of clusters.
#         mode        =       There are two methods implemented to generate the initial condition. 'hc' stands for
#                             hierarchical agglomeration, which is the standarded method used in Mclust. 'denoise'
#                             will return the denoised mclust result as the initial condition.
# Return:
#         A list of 3 elements. initial_means is a cluster*d array of cluster averages. initial_cov is a list whose
#         length is cluster. Each element of initial_cov is a d by d covaraince matrix. initial_tau is a array of 
#         length cluster whose components are mixing proportions.
#
initial_hier = function(sampleMat, cluster, mode = 'hc') {
    if (mode == 'hc') {
        hclass1 = hclass(hc(sampleMat), G=cluster) 
        mat_withlabels <- as.data.frame(sampleMat) %>% mutate(class = hclass1)
        initial_means = mat_withlabels %>% group_by(class) %>%
            summarise_all(mean) %>% select(-class) %>% as.matrix()
        initial_cov = lapply(1:cluster, function(x) {sampleMat[which(hclass1 == x), ] %>% cov()})
        initial_tau = sapply(1:cluster, function(x) {sum(hclass1 == x) / length(hclass1)})
    } else if (mode == 'denoise') {
        NoiseInit_by_guess = sample(c(TRUE,FALSE), size = dim(sim_info$x)[1],
                                    replace = TRUE, prob = c(1,5)/6)
        result_mclust = Mclust(sim_info$x, verbose=F, G=cluster, modelNames = "VVV",
                               initialization = list(noise = NoiseInit_by_guess),
                               control = emControl(tol=c(1.e-6, 1.e-7), itmax=15))
        hclass1 = result_mclust$classification
        mat_withlabels <- as.data.frame(sim_info$x) %>% mutate(class = hclass1)
        initial_means = mat_withlabels %>% group_by(class) %>%
            summarise_all(mean) %>% select(-class) %>% as.matrix()
        initial_means = initial_means[-1,]
        initial_cov = lapply(1:cluster, function(x) {sim_info$x[which(hclass1 == x), ] %>% cov()})
        initial_tau = sapply(1:cluster, function(x) {sum(hclass1 == x) / length(hclass1)})
    } else {
        stop('invalid mode')
    }
    return(list(initial_means, initial_cov, initial_tau))
}







