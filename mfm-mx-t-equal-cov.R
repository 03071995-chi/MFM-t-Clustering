
#####values of the params################
# data = dat_sim
# alpha = (1+dim(dat_sim)[1])/2
# beta = (1+dim(dat_sim)[1])/2
# RHO=1
# a_gamma = 3
# b_gamma = 0.1
# a_nu = 2
# b_nu = 0.1
# #nu = 5
# M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2
# Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ) 
# Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ) 
# initNClusters =3
rm(list=ls()) ### clears thr R workspace
##########################################
#########################################
MFM_Mxt_equal_cov <- function(data, niterations, alpha, beta, RHO, a_gamma, b_gamma, a_nu, b_nu, M0, Sigma0, Omega0, initNClusters, MLE.initial=FALSE)
{
  ## Model: Y_{i}|z,theta \sim MN(theta_{z_i}) ##
  ##        M_{r} \sim MN(M0, Sigma0, Omega0), r = 1,...,k ##
  ##        U ~ IW(2alpha, (2beta)^{-1})
  ##        V ~ IW(2psi, (2rho)^{-1})
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  
  ################################################################
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters (shape, rate) for the prior on elements in lambda in Gamma distribution ##
  ##        RHO = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = possion parameters, a k by 1 vector ##
  
  #################################################################
  # n = length(data)
  n = dim(data)[3]
  #precomputation for prespecified coefficient VN
  N=n ## n is the number of oberservations
  # dimension of the each data observation
  p <- nrow(data[,,1])
  q <- ncol(data[,,1])
  solve_Omega0 <- solve(Omega0)
  solve_Sigma0 <- solve(Sigma0)

  PSI1 <- diag(1,p) ## I_p
  PSI2 <- diag(1,q) ## I_q
  
  lambda = 1 
  # N = 400
  VN<-0
  tmax = 400+10
  for (t in 1:tmax)
  {
    r = log(0)#-inf
    for (k in t:500)
    {
      b = sum(log((k-t+1):k))-sum(log((k*RHO):(k*RHO+N-1))) + dpois(k-1, lambda, log = TRUE)
      m = max(b,r)
      r = log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] = r
  }  
  
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  # need to initialize M, W, Sigma2
  if(MLE.initial){
    MLE_data <- MixMatrix::MLmatrixt(data)
    M <- array(0, dim=c(p,q,initNClusters))
    for(r in 1:initNClusters){
      M[,,r] <- MLE_data$mean
    }
    # the last dimension corresponds to the sample size, initNClusters
    W <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
    Sigma2 <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
  }
  else{
    M <- array(mniw::rMNorm(initNClusters,
                            Lambda = M0,
                            SigmaR = Sigma0, SigmaC = Omega0),
               dim=c(p,q,initNClusters)) # three dimensions
    gamma <- rgamma(1, shape = a_gamma, scale = b_gamma)
    Sigma1 <- rWishart(1,2*alpha,(1/gamma)*solve(PSI1)) ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
    Sigma1 <- matrix(Sigma1,nrow = nrow(Sigma1) ,ncol = ncol(Sigma1))
    nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be >1 bcz nu+p-1>dim( )
    W <- MCMCpack::riwish(nu+p-1,Sigma1)
    Sigma2 <- MCMCpack::riwish(2*beta, gamma*PSI2) ## we multiply 2 with beta to ensure that nu parameter is GREATER than dim(PSI2)
  }
  # U <- MCMCpack::riwish(2*alpha, 2*beta)
  # V <- MCMCpack::riwish(2*psi, 2*rho) 
  solve_Sigma2 = solve(Sigma2)
  solve_W = solve(W)
  # print(M)
  
  # MCMCpack::InvWishart
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    print(iter)
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    # avoid repeated calculations in marginal likelihood part, pre-compute these common numbers in advance
    solve_Sigma_tilde <- solve_Sigma2 %x% solve_W + solve_Omega0 %x% solve_Sigma0
    Sigma_tilde <- solve( solve_Sigma_tilde )
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 ###??
        ###
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(RHO + c.counts.noi[x]) + MixMatrix::dmatrixt(x = data[,,i],df = nu,
                                                    mean = M[,,x], 
                                                    U = Sigma1, 
                                                    V = Sigma2, log=TRUE)
        })                                                       ###log=TRUE gives the logarithm of Mxt density
        clusterAssign_1 = clusterAssign   ##??
        clusterAssign_1[i] = nClusters+1  ##??
        
        Mu_tilde <- Sigma_tilde %*% ( (solve_Sigma2 %x% solve_W ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        
        log_numerator <- log(MCMCpack::diwish(W, nu+p-1, Sigma1)) + (p*q/2) * log(det(Sigma_tilde))  + (1/2) * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - (1/2) * ( t(c(data[,,i])) %*% (solve_Sigma2 %x% solve_W) %*% c(data[,,i]) - (1/2) * t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + (p/2) * log(det(Sigma2)) + (q/2) * log(det(W)) + (p/2) * log(det(Omega0)) + (q/2) * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(RHO) + (VN[nClusters+1]-VN[nClusters]) + (log_numerator - log_denominator)
            
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ) )##??--->logSum part is for scaling?
        print(exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
        clusterAssign[i] <- cluster.i
        #
        if (cluster.i > nClusters) ### when a new cluster arrives
        {
          Mnew = array(0, dim=c(p,q,nClusters+1))
          Mnew[,,1:nClusters] = M
          Mnew[,,nClusters+1] = mniw::rMNorm(1, Lambda = M0, SigmaR = Sigma0, SigmaC = Omega0) # simulate from prior
          M = Mnew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {# phi = phi 
            M = M
            clusterSizes <- table (as.factor(clusterAssign))
            nClusters <- length(clusterSizes)}
      } else {
        if(length(dim(M))==2){
          M = array(M , dim=c(p,q,1))
        }
        # a singleton, have |C| choices 
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - RHO # can offset the gamma adding later
        #finding the probs for sampling process
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(RHO+c.counts.noi[x]) + MixMatrix::dmatrixt(x = data[,,i],df = nu,
                                                         mean = M[,,x], 
                                                         U = Sigma1, 
                                                         V = Sigma2, log=TRUE)
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        #
        Mu_tilde <- Sigma_tilde %*% ( (solve_Sigma2 %x% solve_W ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- log(MCMCpack::diwish(W, nu+p-1, Sigma1)) + (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_Sigma2 %x% solve_W) %*% c(data[,,i]) -1/2* (t(c(M0))) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(Sigma2)) + q/2 * log(det(W)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(RHO) + (VN[nClusters+1]-VN[nClusters]) + (log_numerator - log_denominator)
        # log_clusterProbs[nClusters+1]<-log(GAMMA) + VN[nClusters+1]-VN[nClusters] + dMNorm_marginal(Y=data[,,i], U=U, 
        #                                                                                     V=V, M0=M0, 
        #                                                                                     Sigma0=Sigma0, 
        #                                                                                     Omega0=Omega0,
        #                                                                                     log=TRUE)
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
        print(exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          # phi = phi[-cur.cluster.i]}
          M = M[,,-cur.cluster.i] }
      }
    }
    
    n_Cluster_size <- table( factor(clusterAssign, levels=1:max(clusterAssign)) )
    # ensure M is a 3-D array
    M <- array(M, dim = c(p,q,nClusters))
    for (r in 1:nClusters){
      temp_zY_sum <- matrix(0, nrow=p, ncol=q)
      for(ii in 1:n){
        if(clusterAssign[ii] == r){
          temp_zY_sum <- temp_zY_sum + data[,,ii]
        }
      }
      temp_zeta <- n_Cluster_size[r] * (solve_Sigma2 %x% solve_W) + (solve_Omega0 %x% solve_Sigma0)
      temp_xi <- c( solve_Sigma0 %*% M0 %*% solve_Omega0 + solve_W %*% temp_zY_sum %*% solve_Sigma2 )
      # first generate a vector
      solve_temp_zeta <- solve(temp_zeta)
      # print(solve_temp_zeta %*% temp_xi)
      temp_M_vec <- mvtnorm::rmvnorm(1, mean = solve_temp_zeta %*% temp_xi, sigma = solve_temp_zeta)

      M[,,r] <- matrix(temp_M_vec, nrow = p, ncol = q) # convert a vector to matrix by column
      # clusterAssign is a vector, convert it to Z_ij
      # rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
   
     # update W
    W_bar <- matrix(Sigma1, nrow = nrow(Sigma1) , ncol = ncol(Sigma1))
    #
    for(ii in 1:n){
      W_bar = W_bar + 1 * (data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_Sigma2 %*% t(data[,,ii] - M[,,clusterAssign[ii]])

    }
    W <- MCMCpack::riwish((nu+p+q-1), W_bar)
    # fix the trace of   to be p
    W <- p * W / sum(diag(W))
    solve_W = solve(W) 
    
    #Update sigma1  
    Sigma1_scale <- solve(n*solve(W) + gamma*PSI1)
    df = 2*alpha + n*(nu+p-1)
    Sigma1 <- rWishart(1, df,  Sigma1_scale)
    Sigma1 <- matrix(Sigma1, nrow = nrow(Sigma1), ncol = ncol(Sigma1))
    solve_Sigma1  = solve(Sigma1)
    
    # Update sigma2
    Sigma2_scale <- matrix(0, nrow=q, ncol=q) + gamma*PSI2
    for(ii in 1:n){
      Sigma2_scale <- Sigma2_scale + 1 * t(data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_W %*% (data[,,ii] - M[,,clusterAssign[ii]])
    }
   # Sigma2 <- MCMCpack::riwish(beta+n*p*nClusters, Sigma2_scale)

    Sigma2 <- MCMCpack::riwish(2*beta+n*p, Sigma2_scale)
    # fix the trace of sigma2 to be q
    Sigma2 <- q * Sigma2 / sum(diag(Sigma2))
    solve_Sigma2 = solve(Sigma2)
    
    # Update gamma parameter
    gamma_shape = a_gamma + 1/2 *(2*p*alpha + 2*q*beta)
    gamma_scale = 1/ ((1/b_gamma) + 1/2 * (sum(diag(PSI1 %*% Sigma1)) + sum(diag(PSI2 %*% solve_Sigma2))))
    
    gamma = rgamma(1 , gamma_shape , gamma_scale)
    
    # Update nu parameter 
    # constant = (mvgamma(1/2 * (nu+p+q-1))/ mvgamma(1/2 * (nu+p-1)))**n
    # shape = a_nu
    # scale = 1/b_nu 
    # for(ii in 1:n){
    #   scale <- scale + 1/2 * 1 * log(det(diag(p)+ solve_Sigma1 %*% (data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_Sigma2 %*% t(data[,,ii] - M[,,clusterAssign[ii]])))
    # }
    # scale = 1/scale
    # 
    # #nu = constant * rtgamma(1, shape, scale, min = 1.0001, max = 1e+09)
    # nu = rtgamma(1, shape, scale, min = 1.0001, max = 1e+09)
    
    # print(W)
    # print(Sigma2)
    # print(Sigma1)
    # History[[iter]] <- list(zout = clusterAssign, phiout = phi)
    History[[iter]] <- list(zout = clusterAssign, Mout = M, Wout = W, Sigma2out = Sigma2, Sigma1out = Sigma1)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}
###########################################

# mbm <- microbenchmark(matrixsampling:::rwishart(1, nu = 10, Sigma = diag(4))
# ,Sigma1 <- rWishart(1, 10, diag(4)),Sigma2<-MCMCpack::rwish(10, diag(4)))



