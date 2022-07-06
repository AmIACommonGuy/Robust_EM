library(MASS)
library(Matrix)
library(mvtnorm)
library(mixtools)
library(dplyr)
library(stats)

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







