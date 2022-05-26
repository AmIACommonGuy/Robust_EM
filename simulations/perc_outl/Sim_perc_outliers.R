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

num = 261
set.seed(num)

## parameters
d = 3
n = 100
c = 3
rand_outliers = NULL
num_sim = 40
lam = 2.5

for (i in 1:5) { # for 1:5 we range from 5% to 25%
    # Let's do 50 different sets per k
    
    # Set the percent outliers
    perc_out = i*0.05
    print(i)
    
    for (j in 1:num_sim) {
        # for each k, we want the repeated samples to be the same
        seed_num = num + j
        set.seed(seed_num)
        print(seed_num)
        
        ################ Create the simulated data #################
        # Create the simulation data
        samples = replicate(c, multivarGaussian_unif(n=n, d=d, out_perc=perc_out, out_mag=1)) # changed from 5 to 10
        
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
        result_mclust = Mclust(sampleMat, verbose=F, G=c, modelNames = "VVV",
                               initialization = list(hcPairs = hc(sampleMat)),
                               control = emControl(tol=c(1.e-5, 1.e-6), itmax=15))
        
        true = as.numeric(sampleMat1$Cluster) ## The true label
        
        acc_mcl = rand.index(true, result_mclust$classification)
        acc_std = rand.index(true, as.numeric(result_standard$hard_assign))
        acc_rem = rand.index(true, as.numeric(result_rem$hard_assign))
        
        curr_mclust = c(perc_out*100, seed_num, "mclust", acc_mcl)
        curr_stand = c(perc_out*100, seed_num, "standard", acc_std)
        curr_rob = c(perc_out*100, seed_num, "robust", acc_rem)
        
        rand_outliers = rbind(rand_outliers, curr_stand, curr_rob, curr_mclust)
    }
}

rand_outliers = as.data.frame(rand_outliers)
colnames(rand_outliers) = c("Perc_outliers", "Seed_number", "Type_EM", "Rand_index")
rand_outliers = rand_outliers %>% mutate(Perc_outliers = as.numeric(as.character(Perc_outliers)),
                                         Seed_number = as.numeric(as.character(Seed_number)),
                                         Rand_index = as.numeric(as.character(Rand_index)))

## fix
rand_outliers_mean = rand_outliers %>% group_by(Perc_outliers, Type_EM) %>%
    summarise(mean=mean(Rand_index),
              lower = mean - 1.96*sd(Rand_index)/sqrt(num_sim),
              upper = min(mean + 1.96*sd(Rand_index)/sqrt(num_sim),1))


## Plot the accuracy over the percent outliers
ggplot(rand_outliers_mean, aes(x = Perc_outliers, y = mean, colour = Type_EM)) +
    geom_errorbar(aes(ymin=lower, ymax=upper, colour = NULL, group = Type_EM),
                  width=.4, position = position_dodge(0.6)) +
    geom_line(position = position_dodge(0.6), stat = "identity") +
    geom_point(position = position_dodge(0.6), size=2) +
    labs(y = "Rand Index", x = "Percent Outliers (%)", colour = "Type of EM\nalgorithm") +
    ggtitle("Clustering accuracy over percent outliers") +
    scale_color_brewer(palette="Dark2") +
    theme(plot.title = element_text(hjust = 0.5))

