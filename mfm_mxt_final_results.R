########################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack",
              "mvtnorm","dplyr","armspp","ggcorrplot","reshape2","ggpubr",
              "cowplot","doParallel","foreach","MixMatrix","abind",
              "mclust", "mcclust","fossil","kernlab", "mixtools",
              "npmr","tibble","doMC","clusterGeneration","ClusterR","Rcpp","sna")
# https://cran.r-project.org/web/packages/MixMatrix/vignettes/matrixnormal.html
ipak(packages)
##########################
source("MFM_MN.R")
#
M_list_10by6 <- readRDS("C:/Users/Asus/Desktop/Matrix-t-clustering/bayesian_clustering_supplement/M_10by6_list.rds")
#
MCMC.total <- 100
MCMC.burnin <- ceiling(MCMC.total*0.30)
# MCMC.thin <- 10
#registerDoParallel(cores=8)
# helper function AR(1) correlation matrix
# ar1_cor <- function(n, rho) {
#   exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
#                     (1:n - 1))
#   rho^exponent
# }
###
# http://sherrytowers.com/2013/10/24/k-means-clustering/
# AIC for kmeans
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}
# BIC for kmeans
kmeansBIC <- function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}
###
simulation_settings <- data.frame(expand.grid(signal_strength = 1, 
                                              noise_factor = c(1.5, 1, 0.5), 
                                              total_sample_size = c(100, 200, 300),
                                              rho_factor = c(0.9,0.3),
                                              cluster_num = c(3,6))) %>% rowid_to_column("setting_id")
################
# iii = as.numeric(readline(prompt = "Enter the setting id: "))
iii = 1

#GAMMA = 1
# GAMMA = as.numeric(readline(prompt = "Enter the value of concentration parameter for Dirichlet distribution: "))
################
signal_strength <- simulation_settings$signal_strength[iii]
noise_factor <- simulation_settings$noise_factor[iii]
rho_factor <- simulation_settings$rho_factor[iii]
cluster_num <- simulation_settings$cluster_num[iii]
total_sample_size <- simulation_settings$total_sample_size[iii]
setting_id <- simulation_settings$setting_id[iii]
############################
for(i in 1:cluster_num){
  assign(paste0("M",i), M_list_10by6[[i]])
}


################
### values of params####
# data = dat_sim
# alpha = 5.5
# beta = 5.5
# RHO=1
# a_gamma = 6
# b_gamma = 0.1
# a_nu = 5
# b_nu = 0.1
# M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2
# Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ) 
# Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ) 
# initNClusters =3
# 
# ###True params
# gamma <- rgamma(1, shape = a_gamma, scale = b_gamma)
# Sigma1 <- rWishart(1,2*alpha,(1/gamma)*solve(PSI1)) ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
# Sigma1 <- matrix(Sigma1,nrow = nrow(Sigma1) ,ncol = ncol(Sigma1))
# nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be >1 bcz nu+p-1>dim( )
# W <- MCMCpack::riwish(nu+p-1,Sigma1)
# Sigma2 <- MCMCpack::riwish(2*beta, gamma*PSI2)

# set.seed(0)
num_replicates <- 1
# 
if(cluster_num == 3){
  cluster_proportion <- c(0.15,0.35,0.5)
} else if(cluster_num == 6){
  cluster_proportion <- c(0.35, 0.20, 0.15, 0.1, 0.1, 0.1)
  # rep(1/cluster_num, cluster_num)
}
#registerDoMC(1)
# sim_10by6 <- vector(mode="list", length = length(total_sample_size))
start_time <- Sys.time()
###
sim_10by6 <- foreach(j = 1:num_replicates, .errorhandling = "pass")%dopar%{
  # 
  set.seed(j+123456*setting_id) 
  # W <- clusterGeneration::rcorrmatrix(d = nrow(M1))
  ###########################
  #### true values of params####
  alpha = 5.5
  beta = 5.5
  RHO = 1
  a_gamma = 5
  b_gamma = 0.1
  a_nu = 5
  b_nu = 0.1
  p = 10
  q = 6
  #### generated params #### 
  gamma <- rgamma(1, shape = a_gamma, scale = b_gamma)
  
  PSI1_cluster1= (0.5^2)*diag(1,p)
  PSI1_cluster2= diag(1,p)
  PSI1_cluster3= (1.5^2)*diag(1,p)
  ##
  PSI2_cluster1= (0.5^2)*diag(1,q)
  PSI2_cluster2= diag(1,q)
  PSI2_cluster3= (1.5^2)*diag(1,q)
  ##  
  Sigma1_cluster1 <- matrix(rWishart(1,2*alpha,(1/gamma)*solve(PSI1_cluster1)),nrow = p ,ncol = p) ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
  Sigma1_cluster2 <- matrix(rWishart(1,2*alpha,(1/gamma)*solve(PSI1_cluster2)),nrow = p ,ncol = p)
  Sigma1_cluster3 <- matrix(rWishart(1,2*alpha,(1/gamma)*solve(PSI1_cluster3)),nrow = p ,ncol = p)
  ##
  Sigma2_cluster1 <- MCMCpack::riwish(2*beta, gamma*PSI2_cluster1)
  Sigma2_cluster2 <- MCMCpack::riwish(2*beta, gamma*PSI2_cluster2)
  Sigma2_cluster3 <- MCMCpack::riwish(2*beta, gamma*PSI2_cluster3)
  
  
  nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be >1 bcz nu+p-1>dim( )
  # W <- MCMCpack::riwish(nu+p-1,Sigma1)

  if(cluster_num == 3){
    dat_sim_cluster1 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[1], df = nu,
                                             mean = M1*signal_strength, 
                                             U = Sigma1_cluster1 , 
                                             V = Sigma2_cluster1)  
    dat_sim_cluster2 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[2], df = nu,
                                             mean = M2*signal_strength, 
                                             U = Sigma1_cluster2 , 
                                             V = Sigma2_cluster2)
    dat_sim_cluster3 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[3], df = nu,
                                             mean = M3*signal_strength, 
                                             U = Sigma1_cluster3 , 
                                             V = Sigma2_cluster3)
    # 
    dat_sim <- abind(dat_sim_cluster1, dat_sim_cluster2, dat_sim_cluster3,
                     along = 3)
    # 
    true_cluster_membership <- rep(c(1,2,3), c(total_sample_size*cluster_proportion[1],
                                               total_sample_size*cluster_proportion[2],
                                               total_sample_size*cluster_proportion[3]))
    
  } else if(cluster_num == 6){
    dat_sim_cluster1 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[1], df = nu,
                                             mean = M1*signal_strength, 
                                             U = Sigma1_cluster1 , 
                                             V = Sigma2_cluster1)  
    dat_sim_cluster2 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[2], df = nu,
                                             mean = M2*signal_strength, 
                                             U = Sigma1_cluster2 , 
                                             V = Sigma2_cluster2)
    dat_sim_cluster3 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[3], df = nu,
                                             mean = M3*signal_strength, 
                                             U = Sigma1_cluster3 , 
                                             V = Sigma2_cluster3)
    dat_sim_cluster4 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[4], df = nu,
                                             mean = M4*signal_strength, 
                                             U = Sigma1_cluster1 , 
                                             V = Sigma2_cluster2)
    dat_sim_cluster5 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[5], df = nu,
                                             mean = M5*signal_strength, 
                                             U = Sigma1_cluster3 , 
                                             V = Sigma2_cluster1)
    dat_sim_cluster6 <-  MixMatrix::rmatrixt(total_sample_size*cluster_proportion[6], df = nu,
                                             mean = M6*signal_strength, 
                                             U = Sigma1_cluster2 , 
                                             V = Sigma2_cluster3)
    
    dat_sim <- abind(dat_sim_cluster1, dat_sim_cluster2, dat_sim_cluster3,  dat_sim_cluster4,
                     dat_sim_cluster5, dat_sim_cluster6,
                     along = 3)
    
    true_cluster_membership <- rep(c(1,2,3,4,5,6), c(total_sample_size*cluster_proportion[1],
                                                     total_sample_size*cluster_proportion[2],
                                                     total_sample_size*cluster_proportion[3],
                                                     total_sample_size*cluster_proportion[4],
                                                     total_sample_size*cluster_proportion[5],
                                                     total_sample_size*cluster_proportion[6]))
  }
  # combine them
  dat_sim_vector <- t(apply(dat_sim, c(3), c))
  
  ###params###
  M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2
  Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 )
  Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 )

  # pairwise distance between matrices
  # nuclear norm, spectral norm
  dat_sim_matrix_dist_spectral <- matrix(0, nrow = total_sample_size, ncol = total_sample_size)
  dat_sim_matrix_dist_nuclear <- matrix(0, nrow = total_sample_size, ncol = total_sample_size)
  for(ii in 1:total_sample_size){
    for(jj in 1:total_sample_size){
      if(ii != jj){
        dat_sim_matrix_dist_spectral[ii,jj] <- norm(dat_sim[,,ii] - dat_sim[,,jj], type = "2")
        dat_sim_matrix_dist_nuclear[ii,jj] <- npmr::nuclear(dat_sim[,,ii] - dat_sim[,,jj])
      }
    }
  }
  # dat_sim <- array(NA, dim = c(nrow(M1),ncol(M1),total_sample_size))
  # run MFM analysis
  temp_time1 <- Sys.time()
  
  data = dat_sim
  initNClusters = 3
  niterations = 20
  
  temp_MFM_rlt <- MFM_Mxt_equal_cov(data, niterations, alpha, beta, RHO, a_gamma, b_gamma, a_nu, b_nu, M0, Sigma0, Omega0, initNClusters, MLE.initial=FALSE)

  # 
  ii = length(temp_MFM_rlt$Iterates) 
  RI_MFM_trace <- sapply(1:MCMC.total, function(x) fossil::rand.index(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout))
  ARI_MFM_trace <- sapply(1:MCMC.total, function(x) mclust::adjustedRandIndex(true_cluster_membership,temp_MFM_rlt$Iterates[[ii]]$zout))
  
  ### 
  temp_time2 <- Sys.time()
  temp_DP_rlt <- DP_Mxt_equal_cov(data, niterations, alpha, beta, RHO,a_gamma,b_gamma,a_nu, b_nu, M0, Sigma0, Omega0, initNClusters, MLE.initial=FALSE)

  #
  RI_DP_trace <- sapply(1:MCMC.total, function(x) fossil::rand.index(true_cluster_membership,temp_DP_rlt$Iterates[[ii]]$zout))
  ARI_DP_trace <- sapply(1:MCMC.total, function(x) mclust::adjustedRandIndex(true_cluster_membership,temp_DP_rlt$Iterates[[ii]]$zout))
  #
  temp_time3 <- Sys.time()
  # do.call("c", lapply( ,function(x) fossil::rand.index(true_cluster_membership,x$zout) ) )
  temp_MFM_Dahl_rlt <- getDahl(temp_MFM_rlt,burn=MCMC.burnin)
  temp_DP_Dahl_rlt <- getDahl(temp_DP_rlt,burn=MCMC.burnin)
  ## k-means
  k_vec <- 1 : (total_sample_size - 1)
  temp_kmeans_list <- vector(mode="list", length=length(k_vec))
  temp_kmeans_aic <- rep(NA, length=length(k_vec))
  temp_kmeans_bic <- rep(NA, length=length(k_vec))
  ##
  iters <- temp_MFM_rlt$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)

  for(k in 1:length(k_vec)){
    # kmeans each row is a data point
    temp_kmeans_list[[k]] <- kmeans(dat_sim_vector, 
                                    centers = k_vec[k])
    temp_kmeans_aic[k] <- kmeansAIC(temp_kmeans_list[[k]])
    temp_kmeans_bic[k] <- kmeansBIC(temp_kmeans_list[[k]])
  }
  # 
  temp_kmeans_rlt_aic <- temp_kmeans_list[which(temp_kmeans_aic==min(temp_kmeans_aic))]
  temp_kmeans_rlt_bic <- temp_kmeans_list[which(temp_kmeans_bic==min(temp_kmeans_bic))]
  temp_kmeans_rlt_oracle <- temp_kmeans_list[which(k_vec==cluster_num)]
  #
  if(max(temp_MFM_Dahl_rlt$zout) == total_sample_size){
    temp_kmeans_rlt_MFM <- vector(mode = "list", length = 1)
    temp_kmeans_rlt_MFM[[1]]$cluster <- 1:total_sample_size
  } else{
    temp_kmeans_rlt_MFM <- temp_kmeans_list[which(k_vec==max(temp_MFM_Dahl_rlt$zout))] 
  }
  # temp_kmeans_rlt_MFM <- temp_kmeans_list[which(k_vec==max(temp_MFM_Dahl_rlt$zout))]
  #
  temp_specc_rlt_oracle <- kernlab::specc(dat_sim_vector, centers=cluster_num)
  temp_specc_rlt_MFM <- try(kernlab::specc(dat_sim_vector, centers=max(temp_MFM_Dahl_rlt$zout)))
  # temp_time <- Sys.time()
  
  # finite mixtures of gaussian
  temp_mvnormalmixEM_rlt_oracle <- try(mvnormalmixEM(dat_sim_vector,
                                                     k = cluster_num))
  temp_mvnormalmixEM_rlt_MFM <- try(mvnormalmixEM(dat_sim_vector,
                                                  k = max(temp_MFM_Dahl_rlt$zout)))
  
  # matrix norm k-centroid clustering
  temp_kcentroid_rlt_nuclear_oracle <- Cluster_Medoids(dat_sim_matrix_dist_nuclear,
                                                       clusters = cluster_num)
  temp_kcentroid_rlt_spectral_oracle <- Cluster_Medoids(dat_sim_matrix_dist_spectral,
                                                        clusters = cluster_num)
  #
  temp_kcentroid_rlt_nuclear_MFM <- Cluster_Medoids(dat_sim_matrix_dist_nuclear,
                                                    clusters = cluster_num)
  temp_kcentroid_rlt_spectral_MFM <- Cluster_Medoids(dat_sim_matrix_dist_spectral,
                                                     clusters = cluster_num)
  #########################
  #########################
  # clustering performance
  ARI_MFM <-  mclust::adjustedRandIndex(true_cluster_membership,
                                        temp_MFM_Dahl_rlt$zout)
  ARI_DP <-  mclust::adjustedRandIndex(true_cluster_membership,
                                       temp_DP_Dahl_rlt$zout)
  #
  ARI_kmeans_aic <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_aic[[1]]$cluster)
  ARI_kmeans_bic <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_bic[[1]]$cluster)
  ARI_kmeans_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                 temp_kmeans_rlt_oracle[[1]]$cluster)
  
  ARI_kmeans_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                              temp_kmeans_rlt_MFM[[1]]$cluster)
  
  ARI_specc_MFM <- ifelse(class(temp_specc_rlt_MFM) == "try-error", 
                          NA,
                          mclust::adjustedRandIndex(true_cluster_membership,
                                                    temp_specc_rlt_MFM@.Data))
  
  ARI_specc_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                temp_specc_rlt_oracle@.Data)
  #
  ARI_mvnormalmixEM_oracle <- ifelse(class(temp_mvnormalmixEM_rlt_oracle) == "try-error",
                                     NA,mclust::adjustedRandIndex(true_cluster_membership,
                                                                  apply(temp_mvnormalmixEM_rlt_oracle$posterior,
                                                                        1, function(x) which(x == max(x)))) )
  #  
  ARI_mvnormalmixEM_MFM <- ifelse(class(temp_mvnormalmixEM_rlt_MFM) == "try-error",
                                  NA,mclust::adjustedRandIndex(true_cluster_membership,
                                                               apply(temp_mvnormalmixEM_rlt_MFM$posterior,
                                                                     1, function(x) which(x == max(x)))) )
  # 
  ARI_kcentroid_nuclear_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                            temp_kcentroid_rlt_nuclear_oracle$clusters)
  
  ARI_kcentroid_nuclear_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                                         temp_kcentroid_rlt_nuclear_MFM$clusters)
  
  ARI_kcentroid_spectral_oracle <- mclust::adjustedRandIndex(true_cluster_membership,
                                                             temp_kcentroid_rlt_spectral_oracle$clusters)
  
  ARI_kcentroid_spectral_MFM <- mclust::adjustedRandIndex(true_cluster_membership,
                                                          temp_kcentroid_rlt_spectral_MFM$clusters)
  # 
  Khat_MFM <- max(temp_MFM_Dahl_rlt$zout)
  Khat_DP <- max(temp_DP_Dahl_rlt$zout)
  Khat_kmeans_aic <- max(temp_kmeans_rlt_aic[[1]]$cluster)
  Khat_kmeans_bic <- max(temp_kmeans_rlt_bic[[1]]$cluster)
  Khat_kmeans_oracle <- cluster_num
  Khat_kmeans_MFM <- Khat_MFM
  Khat_specc_oracle <- cluster_num
  Khat_specc_MFM <- Khat_MFM 
  
  Sigma2Sigma1_RMSE_MFM <- sqrt(mean(((temp_MFM_Dahl_rlt$Sigma2out %x% temp_MFM_Dahl_rlt$Sigma1out) - (Sigma2 %x% Sigma1))^2))
  Sigma2Sigma1_RMSE_DP <- sqrt(mean(((temp_DP_Dahl_rlt$Sigma2out %x% temp_DP_Dahl_rlt$Sigma1out) - (Sigma2%x% Sigma1))^2))
  # 
  M_RMSE_MFM <- sqrt(mean((temp_MFM_Dahl_rlt$Mout[,,temp_MFM_Dahl_rlt$zout] - dat_sim)^2))
  M_RMSE_DP <- sqrt(mean((temp_DP_Dahl_rlt$Mout[,,temp_DP_Dahl_rlt$zout] - dat_sim)^2))
  
  temp_rlt <- list(MFM_Dahl_rlt = temp_MFM_Dahl_rlt,
                   DP_Dahl_rlt = temp_DP_Dahl_rlt,
                   kmeans_rlt_aic = temp_kmeans_rlt_aic[[1]],
                   kmeans_rlt_bic = temp_kmeans_rlt_bic[[1]],
                   kmeans_rlt_oracle = temp_kmeans_rlt_oracle[[1]],
                   MFM_DP_time = difftime(temp_time3, temp_time2, units = "mins"),
                   MFM_Dahl_time = difftime(temp_time2,temp_time1,units = "mins"),
                   Khat_MFM=Khat_MFM,
                   Khat_DP=Khat_DP,
                   Khat_kmeans_aic=Khat_kmeans_aic,
                   Khat_kmeans_bic=Khat_kmeans_bic,
                   Khat_kmeans_oracle=Khat_kmeans_oracle,
                   Khat_kmeans_MFM=Khat_kmeans_MFM,
                   Khat_specc_oracle=Khat_specc_oracle,
                   Khat_specc_MFM=Khat_specc_MFM,
                   ARI_MFM=ARI_MFM,
                   ARI_DP=ARI_DP,
                   ARI_kmeans_aic=ARI_kmeans_aic,
                   ARI_kmeans_bic=ARI_kmeans_bic,
                   ARI_kmeans_oracle=ARI_kmeans_oracle,
                   ARI_kmeans_MFM=ARI_kmeans_MFM,
                   ARI_specc_oracle=ARI_specc_oracle,
                   ARI_specc_MFM=ARI_specc_MFM,
                   ARI_mvnormalmixEM_oracle = ARI_mvnormalmixEM_oracle,
                   ARI_mvnormalmixEM_MFM = ARI_mvnormalmixEM_MFM,
                   ARI_kcentroid_nuclear_oracle = ARI_kcentroid_nuclear_oracle,
                   ARI_kcentroid_nuclear_MFM = ARI_kcentroid_nuclear_MFM,
                   ARI_kcentroid_spectral_oracle = ARI_kcentroid_spectral_oracle,
                   ARI_kcentroid_spectral_MFM = ARI_kcentroid_spectral_MFM,
                   RI_MFM_trace=RI_MFM_trace,
                   ARI_MFM_trace=ARI_MFM_trace,
                   VU_RMSE_MFM=VU_RMSE_MFM,
                   VU_RMSE_DP=VU_RMSE_DP,
                   M_RMSE_MFM=M_RMSE_MFM,
                   M_RMSE_DP=M_RMSE_DP)
  #
  temp_rlt
}
#####################################
end_time <- Sys.time()
#####################################
# saveRDS(sim_10by6,
#         paste0("sim_10by6_",setting_id,".rds"))
#####################################
#####################################

