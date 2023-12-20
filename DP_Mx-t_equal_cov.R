#### MFM for matrix normal
#################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack","mvtnorm","sna","MixMatrix","LaplacesDemon","cascsim","CholWishart")
ipak(packages)

library(LaplacesDemon)


#################################
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  MFMfit = temp_MFM_rlt
  burn = MCMC.burnin
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

#####values of the params################
data = dat_sim
alpha = (1+dim(dat_sim)[1])/2
beta = (1+dim(dat_sim)[1])/2
RHO=1
a_gamma = 3
b_gamma = 0.1
a_nu = 2
b_nu = 0.1
#nu = 5
M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2
Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ) 
Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ) 
initNClusters =3
##########################################
################  
################
## function for Collapsed sampler
DP_Mxt_equal_cov <- function(data, niterations, alpha, beta, RHO,a_gamma,b_gamma,a_nu, b_nu, M0, Sigma0, Omega0, initNClusters, MLE.initial=FALSE)
{             
  ## Model: Y_{i}|M_{z_i},W_{i},Sigma2 \sim MN(M_{z_i},W_{i},Sigma2) ##
  ##        M_{r} \sim MN(M0, Sigma0, Omega0), r = 1,...,k ##
  ##        W_{i} ~ IW_{p}(nu+p-1,Sigma1) , i=1,2...n
  ##        P(z_i = j) = \pi_j, j = 1,...,k 
  ##        Sigma1 ~ W_{p}(alpha, (gamma*PSI_1)^{-1})
  ##        Sigma2 ~ IW_{q}(beta, gamma*PSI_2)
  ##        \pi ~ Dirichlet_{k}(RHO,...,RHO) 
  ##        gamma ~ Gamma(a_{gamma},b_{gamma})
  ##        nu ~ Truncated_gamma(a_{nu},b_{nu})
  ##        k-1 ~ poisson(1) ##
  
  # # define a helper function to calculate the marginal probability Y|W_i,Sigma2 
  # dMNorm_marginal <- function(Y, W, Sigma2, M0, Sigma0, Omega0, log=FALSE){
  #   p <- nrow(Y)
  #   q <- ncol(Y)
  #   solve_Sigma2 = solve(Sigma2)
  #   solve_W = solve(W)
  #   solve_Sigma0 = solve(Sigma0)
  #   solve_Omega0 = solve(Omega0)
  #   Sigma_tilde <- solve( solve_Sigma2 %x% solve_W + solve_Omega0 %x% solve_Sigma0 )
  #   Mu_tilde <- Sigma_tilde %*% ( (solve_Sigma2 %x% solve_W ) %*% c(Y) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
  #   log_numerator <- (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve(Sigma_tilde) %*% Mu_tilde - 1/2 * ( t(c(Y)) %*% (solve_Sigma2 %x% solve_W) %*% c(Y) + t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
  #   log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(Sigma2)) + q/2 * log(det(W)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
  #   if(log == TRUE){
  #     return(log_numerator-log_denominator)
  #   }
  #   else{
  #     return( exp(log_numerator - log_denominator) )
  #   }
  # }
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters (shape, rate) for the prior on elements in lambda in Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distribution ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = poisson parameters, a k by 1 vector ##
  
  #################################################################
  # n = length(data)
  n = dim(data)[3]
  #precomputation for prespecified coefficient VN
  # lambda <- LAMBDA
  #gamma <- GAMMA
  N=n ## n is the number of oberservations
  
  # dimension of the each data observation
  p <- nrow(data[,,1])
  q <- ncol(data[,,1])
  solve_Omega0 <- solve(Omega0)
  solve_Sigma0 <- solve(Sigma0)
  
  PSI1 <- diag(1,nrow = p)
  PSI2 <- diag(1,nrow = q)
  #
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
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        ###
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(c.counts.noi[x]) + MixMatrix::dmatrixt(x = data[,,i],df = nu,
                                                     mean = M[,,x], 
                                                     U = Sigma1, 
                                                     V = Sigma2, log=TRUE)  # log==TRUE gives the logarithm of Mxt density
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        #
        Mu_tilde <- Sigma_tilde %*% ( (solve_Sigma2 %x% solve_W ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- log(MCMCpack::diwish(W, nu+p-1, Sigma1)) + (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_Sigma2 %x% solve_W) %*% c(data[,,i]) -1/2 * t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(Sigma2)) + q/2 * log(det(W)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(RHO) + (log_numerator - log_denominator)
        # + (VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ) )
        clusterAssign[i] <- cluster.i
        #
        
        if (cluster.i > nClusters)
        {
          # phinew = rep(0,nClusters+1)
          # phinew[1:nClusters] = phi
          # phinew[nClusters+1] = rgamma(1, shape = alpha, rate = beta)
          # phi = phinew
          Mnew = array(0, dim=c(p,q,nClusters+1))
          Mnew[,,1:nClusters] = M
          Mnew[,,nClusters+1] = mniw::rMNorm(1, Lambda = M0, SigmaR = Sigma0, SigmaC = Omega0) # simulate from prior
          M = Mnew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {# phi = phi
            M = M
            clusterSizes <- table(as.factor(clusterAssign))
            nClusters <- length(clusterSizes)}
      } else {
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 # can offset the gamma adding later
        #finding the probs for sampling process
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(c.counts.noi[x]) + MixMatrix::dmatrixt(x = data[,,i],df = nu,
                                                     mean = M[,,x], 
                                                     U = Sigma1, 
                                                     V = Sigma2, log=TRUE)
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        # 
        Mu_tilde <- Sigma_tilde %*% ( (solve_Sigma2 %x% solve_W ) %*% c(data[,,i]) + (solve_Omega0 %x% solve_Sigma0) %*% c(M0) )
        log_numerator <- log(MCMCpack::diwish(W, nu+p-1, Sigma1)) + (p*q/2) * log(det(Sigma_tilde))  + 1/2 * t(Mu_tilde) %*% solve_Sigma_tilde %*% Mu_tilde - 1/2 * ( t(c(data[,,i])) %*% (solve_Sigma2 %x% solve_W) %*% c(data[,,i]) -1/2* t(c(M0)) %*% (solve_Omega0 %x% solve_Sigma0) %*% c(M0) ) 
        log_denominator <- (p*q/2) * log(2*pi) + p/2 * log(det(Sigma2)) + q/2 * log(det(W)) + p/2 * log(det(Omega0)) + q/2 * log(det(Sigma0))
        
        log_clusterProbs[nClusters+1] <- log(RHO) + (log_numerator - log_denominator)
        # (VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
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
        # temp_zY_sum <- temp_zY_sum + Z[ii,r] * data[,,i]
      }
      # print(temp_zY_sum)
      temp_zeta <- n_Cluster_size[r] * (solve_Sigma2 %x% solve_W) + (solve_Omega0 %x% solve_Sigma0)
      temp_xi <- c( solve_Sigma0 %*% M0 %*% solve_Omega0 + solve_W %*% temp_zY_sum %*% solve_Sigma2 )
      # first generate a vector
      solve_temp_zeta <- solve(temp_zeta)
      # print(solve_temp_zeta %*% temp_xi)
      temp_M_vec <- mvtnorm::rmvnorm(1, mean = solve_temp_zeta %*% temp_xi, sigma = solve_temp_zeta)
      # browser()
      #print(temp_M_vec)
      #print(r)
      #print(nClusters)
      #print(dim(M))
      #print(p)
      #print(q)
      M[,,r] <- matrix(temp_M_vec, nrow = p, ncol = q) # convert a vector to matrix by column
      # clusterAssign is a vector, convert it to Z_ij
      # rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
    # update W
    W_bar <- matrix(0, nrow=p, ncol=p) + Sigma1
    #
    for(ii in 1:n){
      W_bar <- W_bar + 1 * (data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_Sigma2 %*% t(data[,,ii] - M[,,clusterAssign[ii]])
    }
    W <- MCMCpack::riwish(nu+p+q-1, W_bar)
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
    # History[[iter]] <- list(zout = clusterAssign, phiout = phi)
    History[[iter]] <- list(zout = clusterAssign, Mout = M, Wout = W, Sigma2out = Sigma2, Sigma1out = Sigma1)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}
#########################################
#########################################


####don't need dnorm marginal function####



