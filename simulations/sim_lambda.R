library(MASS)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(ggpubr)
library(mixtools)
library(RColorBrewer)
library(fossil)
library(dplyr)
library(mclust)

num = 2231
set.seed(num)

## parameters
d = 4
n = 100
c = 3
rand_outliers = NULL
num_sim = 50

for (i in 1:4) {
    # Let's do 50 different sets per k
    
    lam <- 2+(i)*0.2
    
    for (j in 1:num_sim) {
        # for each k, we want the repeated samples to be the same
        seed_num = num + j
        set.seed(seed_num)
        print(seed_num)
        
        ################ Create the simulated data #################
        # Create the simulation data
        samples = replicate(c, multivarGaussian_unif(n=n, d=d, out_perc=0.15, out_mag=1)) # changed from 5 to 10
        
        # Combine the c samples so that all samples are in one matrix
        sampleMat = samples[,1]$gauss
        
        if (c > 1) {
            for (k in 2:c) {sampleMat = rbind(sampleMat, as.matrix(samples[,k]$gauss))}
        }
        # Combine the c mu's into one matrix
        sampleMu = t(sapply(1:c, function(x) samples[,x]$mu))
        
        # Combine the c sigmas into one list
        sampleSigma = lapply(1:c, function(x) as.matrix(samples[2,][[x]]))
        
        # Label the sample's clusters
        sampleMat1 = as.data.frame(cbind(sampleMat, rep(1:c, each=n)))
        colnames(sampleMat1) = c("X", "Y", "Cluster")
        sampleMat1$Cluster = as.character(sampleMat1$Cluster)
        init_cond <- initial_hier(sampleMat, cluster = c)
        result_rem = EM_robust(sampleMat =  sampleMat, c = c, d=d, lambda = lam, sigma_str = 'unstr', inits = init_cond)
        result_standard = EM_robust(sampleMat =  sampleMat, c = c, d=d, lambda = 1, sigma_str = 'unstr', inits = init_cond)
        
        true = as.numeric(sampleMat1$Cluster) ## The true label
        
        acc_std = rand.index(true, as.numeric(result_standard$hard_assign))
        acc_rem = rand.index(true, as.numeric(result_rem$hard_assign))
        
        curr_stand = c(lam, seed_num, "standard", acc_std)
        curr_rob = c(lam, seed_num, "robust", acc_rem)
        
        rand_outliers = rbind(rand_outliers, curr_stand, curr_rob)
    }
}

rand_outliers = as.data.frame(rand_outliers)
colnames(rand_outliers) = c("lambda", "Seed_number", "Type_EM", "Rand_index")
rand_outliers = rand_outliers %>% mutate(lambda = as.numeric(as.character(lambda)),
                                         Seed_number = as.numeric(as.character(Seed_number)),
                                         Rand_index = as.numeric(as.character(Rand_index)))

rand_outliers_mean = rand_outliers %>% group_by(lambda, Type_EM) %>%
    summarise(mean=mean(Rand_index),
              lower = mean - 1.96*sd(Rand_index)/num_sim,
              upper = min(mean + 1.96*sd(Rand_index)/num_sim,1))


## Plot the accuracy over the percent outliers
ggplot(rand_outliers_mean, aes(x = lambda, y = mean, colour = Type_EM)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, colour = NULL, group = Type_EM),
                  width=.4) +
    geom_line(stat = "identity") +
    geom_point(size=2) +
    labs(y = "Rand Index", x = "lambda", colour = "Type of EM\nalgorithm") +
    ggtitle("Clustering accuracy over lambda") +
    scale_color_brewer(palette="Dark2") +
    theme(plot.title = element_text(hjust = 0.5))

