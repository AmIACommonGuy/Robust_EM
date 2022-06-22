cluster = 2
lam = 2.5
n = 500
am <- c()
sm <- c()
rm <- c()
for (i in c(1:100)) {
    sim_info <- simMultGauss_unif(n = n, 
                                  d = 2,
                                  cluster = cluster,
                                  out_perc = 0.1,
                                  out_mag = 1.4)
    sampleMat = as.data.frame(cbind(sim_info$simdata, rep(1:cluster, each=n)))
    colnames(sampleMat) = c("X", "Y", "Cluster")
    sampleMat$Cluster = as.character(sampleMat$Cluster)
    true = sim_info$Label
    true_outl <- sim_info$LabelO
    init_cond <- initial_hier(sim_info$simdata, cluster = cluster)
    NoiseInit_by_guess = sample(c(TRUE,FALSE), size = n*cluster, 
                                replace = TRUE, prob = c(1,5)/6)
    result_mclust = Mclust(sim_info[["simdata"]], verbose=F, G=cluster, modelNames = "VVV",
                           initialization = list(hc(sampleMat), noise = NoiseInit_by_guess),
                           control = emControl(tol=c(1.e-6, 1.e-7), itmax=15))
    ## Standard EM
    result_std <- EM_robust(sampleMat =  sim_info$simdata, c = cluster, d=d, lambda = Inf, sigma_str = 'unstr', inits = init_cond)
    ## Robust EM
    result_robust <- EM_robust(sampleMat =  sim_info$simdata, c = cluster, d=d, lambda = lam, sigma_str = 'unstr', inits = init_cond)
    ## 1. check adjusted randindex
    (acc_mcl = rand.index(true_outl, result_mclust$classification))
    am <- c(am,acc_mcl)
    (acc_std = rand.index(true_outl, as.numeric(result_std$hard_assign_out)))
    sm <- c(sm, acc_std)
    (acc_rem = rand.index(true_outl, as.numeric(result_robust$hard_assign_out)))
    rm <- c(rm, acc_rem)
}
boxplot(am,rm,sm,names = c('mclust_noise','standard','robust'), ylab = 'Rand Index')
