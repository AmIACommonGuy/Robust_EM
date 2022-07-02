library(MASS)

plot_mc <- function(x) {
    mu = as.data.frame(t(x$parameters$mean))
    sigma = x$parameters$variance$sigma
    c = x$G
    d = x$d
    n = x$n
    # Assign the cluster number to whichever one has the highest probability
    sampleMat = as.data.frame(x$data)
    colnames(sampleMat) <- c("Axis1", "Axis2")
    sampleMat$Cluster = as.character(x$classification)
    inliers_ind <- which(x$classification !=0)
    plot.new()
    ellipse0 = lapply(1:nrow(mu), function(j) mixtools::ellipse(mu=x$parameters$mean[,j],
                                                                sigma=x$parameters$variance$sigma[,,j], npoints = 300, newplot=F))
    ellipse_list <- c()
    covar_plot_f_EM <- function(j){
        boundary <- data.frame(ellipse0[j])
        colnames(boundary) <- c("Axis1", "Axis2")
        return(geom_path(data = boundary, aes(x=Axis1, y=Axis2), size=0.5, color = 'red'))
    }
    covar_EM_plots <- lapply(1:nrow(mu), covar_plot_f_EM)
    dev.off()
    ggplot() +
        geom_point(aes(x = sampleMat[inliers_ind,1], y = sampleMat[inliers_ind,2]),color= sampleMat$Cluster[inliers_ind]) +
        geom_point(aes(x = sampleMat[-inliers_ind,1], y = sampleMat[-inliers_ind,2]),pch = 4) +
        geom_point(data = mu, aes(x=V1, y=V2), shape=2, size = 1, stroke = 2, color = 'black') +
        covar_EM_plots
}

plot_rem = function(sim_EMfit, sampleMat){
    
    swap <- function(vec, from, to) {
        tmp <- to[ match(vec, from) ]
        tmp[is.na(tmp)] <- vec[is.na(tmp)]
        return(tmp)
    }
    
    # Assign values from returned list to mu, sigma, T_mat, tau, e
    mu = as.data.frame(sim_EMfit$mu)
    sigma = sim_EMfit$sigma
    T_mat = sim_EMfit$T_mat
    tau = sim_EMfit$tau
    #error = sim_EMfit$e
    c = sim_EMfit$c
    d = sim_EMfit$d
    n = sim_EMfit$n
    lambda = sim_EMfit$lambda
    inliers <- which(sim_EMfit$hard_assign_out != 0)
    
    # Think about assigning the points to the clusters
    Soft_assign = t(T_mat)
    
    colnames(Soft_assign) = paste(1:c)
    
    # Assign the cluster number to whichever one has the highest probability
    hard_assign = matrix(paste(apply(Soft_assign, 1, which.max)), n, 1)
    #hard_assign = swap(hard_assign, c(1,2), c(2,1))
    
    sampleMat2 = as.data.frame(sampleMat)
    sampleMat2$Cluster =sim_EMfit$hard_assign_out
    
    plot.new()
    ellipse0 = lapply(1:c, function(j) mixtools::ellipse(mu=sim_EMfit$mu[j,],
                                                         sigma=sim_EMfit$sigma[[j]], npoints = 300, newplot=F))
    
    ellipse_list <- c()
    covar_plot_f_EM <- function(j){
        boundary <- data.frame(ellipse0[j])
        colnames(boundary) <- c("Axis1", "Axis2")
        return(geom_path(data = boundary, aes(x=Axis1, y=Axis2), size=0.5, color = 'red'))
    }
    covar_EM_plots <- lapply(1:nrow(mu), covar_plot_f_EM)
    dev.off()
    g <- ggplot() + geom_point(aes(x = sampleMat[inliers, 1], y = sampleMat[inliers, 2]), color = sampleMat2$Cluster[inliers]) +
        geom_point(aes(x = sampleMat[-inliers, 1], y = sampleMat[-inliers, 2]), pch = 4) +
        geom_point(data = mu, aes(x=V1, y=V2), shape=2, size = 1, stroke = 2, color = 'black') +
        covar_EM_plots
    return(g)
}














