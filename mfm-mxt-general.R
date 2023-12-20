#### MFM for matrix t
#################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("matrixNormal","mniw","MCMCpack","mvtnorm","sna","MixMatrix")
ipak(packages)
#################################
#####values of the params################
data = dat_sim
alpha = (1+dim(dat_sim)[1])/2
beta = (1+dim(dat_sim)[1])/2
RHO=1
a_gamma = 6
b_gamma = 0.1
a_nu = 5
b_nu = 0.1
#nu = 5
M0=(apply(dat_sim, c(1,2), max) + apply(dat_sim, c(1,2), min))/2
Sigma0=diag( ((apply(dat_sim, c(1), max)-apply(dat_sim, c(1), min))/4)^2 ) 
Omega0=diag( ((apply(dat_sim, c(2), max)-apply(dat_sim, c(2), min))/4)^2 ) 
initNClusters =3
#########################################
### general case: each component has its own 
### covariance matrices
#########################################
MFM_Mxt_general <- function(data, niterations, alpha, beta, RHO, a_gamma,b_gamma,a_nu,b_nu, M0, Sigma0, Omega0, initNClusters, 
                            VN, m = 5, MLE.initial=FALSE)
{
  ## Model: Y_{i}|z,theta \sim MN(theta_{z_i}) ##
  ##        M_{r} \sim MN(M0, Sigma0, Omega0), r = 1,...,k 
  ##        U_r ~ IW(2alpha, (2beta)^{-1})
  ##        V_r ~ IW(2psi, (2rho)^{-1})
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-MxN ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        M0 = the parameter for matrix normal distribution ##
  ##        initNClusters = the initial number of clusters ##
  ##        m = additional components needed at each step when updating cluster specific parameters ##
  ##
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = matrix normal parameters, a k by 1 vector ##
  
  #################################################################
  # n = length(data)
  n = dim(data)[3]
  N = n ## n is the number of oberservations
  m=5
  ##################################################################
  ## Calculate VN(t)
  # gamma = 1; 
  lambda = 1 
  
  # N = 400
  VN<-0
  tmax = 400+10
  for (t in 1:tmax)
  {
    r = log(0)
    for (k in t:500)
    {
      b = sum(log((k-t+1):k))-sum(log((k*RHO):(k*RHO+N-1))) + dpois(k-1, lambda, log = TRUE)
      m = max(b,r)
      r = log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] = r
  }
  ##################################################################
  # dimension of the each data observation
  p <- nrow(data[,,1])
  q <- ncol(data[,,1])
  solve_Omega0 <- solve(Omega0)
  solve_Sigma0 <- solve(Sigma0)
  
  PSI1 <- diag(1,p)
  PSI2 <- diag(1,q)
  
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE)) # make sure no component is empty
  #
  # print(clusterAssign)
  
  # need to intialize M, Wi(each Wi for each Yi), Sigma2
  # M: p*q*initNClusters
  # Wi: p*p*n
  W_list <- array(NA, dim = c(p,p,n))
  solve_W_list <- array(NA, dim = c(p,p,n))
  
  #
  if(MLE.initial){
    # MLE for all the data
    MLE_data <- MixMatrix::MLmatrixt(data)
    # initialize M
    M <- array(0, dim=c(p,q,initNClusters))
    #
    W_list <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
    Sigma2 <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V))) ### U corresponds to W, V corresponds to Sigma2. MLE_data(MXnormal) has two cov matrices $U and $V. Check MLE_data$U and MLE_data$V
    for(r in 1:initNClusters){
      # calculate cluster specific MLE
      temp_MLE_data <- try(MixMatrix::MLmatrixt(data[,,which(clusterAssign == r)]))
      # 
      if(class(temp_MLE_data) == "try-error"){
        # if error, use the MLE from general results
        M[,,r] <- MLE_data$mean
        W_list[,,r] <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
        Sigma2 <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
        
      } else{
        # if not, use the MLE from cluster specific results
        M[,,r] <- temp_MLE_data$mean
        W_list[,,r] <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
        Sigma2 <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
      }
      
      # M[,,r] <- MLE_data$mean
    }
    # the last dimension corresponds to the sample size, initNClusters
    # U <- (sum(diag(MLE_data$V))/nrow(MLE_data$V))*MLE_data$U
    # V <- nrow(MLE_data$V)*((MLE_data$V)/sum(diag(MLE_data$V)))
    # 
  }
  else{
    M <- array(mniw::rMNorm(initNClusters,
                            Lambda = M0,
                            SigmaR = Sigma0, 
                            SigmaC = Omega0),
               dim=c(p,q,initNClusters)) # three dimensions, random initializations
    gamma <- rgamma(1, shape = a_gamma, scale = b_gamma)
    Sigma2 <- MCMCpack::riwish(2*beta, gamma*PSI2)
    for(i in 1:n){
      set.seed(i)
      # W
      Sigma1 <- rwishart(2*alpha,(1/gamma)*solve(PSI1)) ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
      nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be >1 bcz nu+p-1>dim(sigma1)
      W_list[,,i] <- MCMCpack::riwish(nu+p-1, Sigma1)
      solve_W_list[,,i] <- solve(W_list[,,i])
    }

  }
  # initialization end
  #
  # browser()
  # solve_V = solve(V)
  # solve_U = solve(U)
  # 
  # 
  # print(M)
  # MCMCpack::InvWishart
  History <- vector("list", niterations)
  ##
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    print(iter)
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    
    # avoid repeated calculations in marginal likelihood part, pre-compute these common numbers in advance
    # solve_Sigma_tilde <- solve_V %x% solve_U + solve_Omega0 %x% solve_Sigma0
    # Sigma_tilde <- solve( solve_Sigma_tilde )
    
    ## 
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      print(i)
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
        ##
        # print(dim(M))
        log_clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          log(RHO+c.counts.noi[x]) + MixMatrix::dmatrixt(x = data[,,i],df = nu,
                                                         mean = M[,,x], 
                                                         U = Sigma1, 
                                                         V = Sigma2, log=TRUE) # dpois(data[i],phi[x])
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        
        # new cluster prob
        # auxiliary variable based methods
        # temp_M_list <- vector(mode="list",length=m)
        # temp_W_list <- vector(mode="list",length=m)
        # temp_Sigma2_list <- vector(mode="list",length=m)
        
        temp_M_list <- list()
        temp_Sigma1_list <-list()
        temp_W_list <- list()
        temp_Sigma2_list <- list()
        #
        for(jj in 1:m){
          # sample from the priors
          temp_M <- mniw::rMNorm(1,
                                 Lambda = M0,
                                 SigmaR = Sigma0, SigmaC = Omega0)
          # 
          temp_Sigma1 <- rwishart(2*alpha,(1/gamma)*solve(PSI1))  ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
          nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be > 1 because nu+p-1 > dim(sigma1)
          temp_W<- MCMCpack::riwish(nu+p-1, temp_Sigma1)
          #
          temp_Sigma2 <- MCMCpack::riwish(2*beta, gamma*PSI2) 
          
          # matrix normal likelihood?
          log_clusterProbs[nClusters+jj] <- log(RHO) + (VN[nClusters+1]-VN[nClusters]) - log(m) + 
            MixMatrix::dmatrixt(x = data[,,i],df = nu, mean = temp_M, 
                                U = temp_Sigma1, 
                                V = temp_Sigma2, log=TRUE)
                                                                                                                              
          #
          temp_M_list[[jj]] <- temp_M
          temp_Sigma1_list[[jj]] <- temp_Sigma1
          temp_W_list[[jj]] <- temp_W
          temp_Sigma2_list[[jj]] <- temp_Sigma2
        }
        
             
        #choose the cluster number for ith observation
        # use logSum to avoid underflow/overflow
        cluster.i <- sample(1:(nClusters+m), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ) ) 
        clusterAssign[i] <- ifelse(cluster.i <= nClusters, cluster.i, nClusters + 1)
        # 
        # 
        if (cluster.i > nClusters)
        {
          # update M
          Mnew = array(0, dim=c(p,q,nClusters+1))
          Mnew[,,1:nClusters] = M
          Mnew[,,nClusters+1] = temp_M_list[[cluster.i-nClusters]]
          M = Mnew
          
          # update Sigma1
          Sigma1new = array(0, dim=c(p,p,nClusters+1))
          Sigma1new[,,1:nClusters] = Sigma1
          Sigma1new[,,nClusters+1] = temp_Sigma1_list[[cluster.i-nClusters]]
          Sigma1 = Sigma1new
          
          # update W_list
          Wnew = array(0, dim=c(p,p,nClusters+1))
          Wnew[,,1:nClusters] = W_list
          Wnew[,,nClusters+1] = temp_W_list[[cluster.i-nClusters]]
          W_list = Wnew
          
          # update Sigma2
          Sigma2new = array(0, dim=c(q,q,nClusters+1))
          Sigma2new[,,1:nClusters] = Sigma2
          Sigma2new[,,nClusters+1] = temp_Sigma2_list[[cluster.i-nClusters]]
          Sigma2 = Sigma2new
          
          # 
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {# phi = phi
            M = M
            W_list = W_list
            Sigma2 = Sigma2
            clusterSizes <- table(as.factor(clusterAssign))
            nClusters <- length(clusterSizes)}
      } else {
        #
        if(length(dim(M))==2){
          M = array(M, dim=c(p,q,1))
        }
        if(length(dim(Sigma1[[1]]))==2){
          Sigma1[[1]] = array(Sigma1[[1]], dim=c(p,p,1))
        }
        if(length(dim(W_list[[1]]))==2){
          W_list[[1]] = array(W_list[[1]], dim=c(p,p,1))
        }
        if(length(dim(Sigma2[[1]]))==2){
          Sigma2[[1]] = array(Sigma2[[1]], dim=c(q,q,1))
        }
        #
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - RHO # can offset the gamma adding later
        #finding the probs for sampling process
        #
        # print(dim(M))
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
        
        # new cluster prob
        # auxiliary variable based methods
        # temp_M_list <- vector(mode="list",length=m)
        # temp_W_list <- vector(mode="list",length=m)
        # temp_Sigma2_list <- vector(mode="list",length=m)
        temp_M_list <- list()
        temp_Sigma1_list <-list()
        temp_W_list <- list()
        temp_Sigma2_list <- list()
        #
        m=5
        #
        for(jj in 1:m){
          # sample from the priors
          temp_M <- mniw::rMNorm(1,
                                 Lambda = M0,
                                 SigmaR = Sigma0, 
                                 SigmaC = Omega0)
          # 
          temp_Sigma1 <- rwishart(2*alpha,(1/gamma)*solve(PSI1)) ## we multiply 2 with alpha to ensure that nu parameter is GREATER than dim(PSI1)
          nu <- rtgamma(1, a_nu, b_nu, min = 1.0001, max = 1e+09) ## nu must be >1 bcz nu+p-1>dim(sigma1)
          temp_W<- MCMCpack::riwish(nu+p-1, Sigma1)
          #
          temp_Sigma2 <- MCMCpack::riwish(2*beta, gamma*PSI2)  
          # matrix normal likelihood?
          log_clusterProbs[nClusters+jj] <- log(RHO) + (VN[nClusters+1]-VN[nClusters]) - log(m) + 
            MixMatrix::dmatrixt(x = data[,,i],df = nu, mean = temp_M, 
                                U = temp_Sigma1, 
                                V = temp_Sigma2, log=TRUE)
          #
          temp_M_list[[jj]] <- temp_M
          temp_Sigma1_list[[jj]] <- temp_Sigma1
          temp_W_list[[jj]] <- temp_W
          temp_Sigma2_list[[jj]] <- temp_Sigma2
        }
        
        #choose the cluster number for ith observation3
        cluster.i <- sample(1:(nClusters+m), size = 1,
                            prob = exp(log_clusterProbs - sna::logSum(log_clusterProbs) ))
        # 
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleton one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels

        } else {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          # phi = phi[-cur.cluster.i]}
          M = M[,,-cur.cluster.i]
          Sigma1 = Sigma1[,,-cur.cluster.i]
          W_list = W_list[,,-cur.cluster.i]
          Sigma2 = Sigma2[,,-cur.cluster.i] 
          
          if(length(dim(M))==2){
            M = array(M, dim=c(p,q,1))
          }
          if(length(dim(Sigma1[[1]]))==2){
            Sigma1 = array(Sigma1, dim=c(p,p,1))
          }
          if(length(dim(W_list[[1]]))==2){
            W = array(W, dim=c(p,p,1))
          }  
          if(length(dim(Sigma2[[1]]))==2){
            Sigma2 = array(Sigma2, dim=c(q,q,1))
          } }
      }
    } # end for loop over subjects i
    

    n_Cluster_size <- table( factor(clusterAssign, levels=1:nClusters) )
    
    ##  solve_U, solve_V calculation to avoid redundant calculations ##
    solve_Sigma1 <- array(NA, dim=c(p,p,nClusters))
    solve_W_list <- array(NA, dim=c(p,p,nClusters))
    solve_Sigma2 <- array(NA, dim=c(q,q,nClusters))
    ## 
    for(r in 1:nClusters){
      solve_Sigma1[,,r] <- solve(Sigma1[,,r])
      solve_W_list[,,r] <- solve(W_list[,,r])
      solve_Sigma2[,,r] <- solve(Sigma2[,,r])
    }
    
    # ensure M is a 3-D array
    # update M
    M <- array(M, dim = c(p,q,nClusters))
    for (r in 1:nClusters){
      temp_zY_sum <- matrix(0, nrow=p, ncol=q)
      for(ii in 1:n){
        if(clusterAssign[ii] == r){
          temp_zY_sum <- temp_zY_sum + data[,,ii]
        }
      }
      # print(temp_zY_sum)
      temp_zeta <- n_Cluster_size[r] * (solve_Sigma2[,,r] %x% solve_W_list[,,r]) + (solve_Omega0 %x% solve_Sigma0)
      temp_xi <- c( solve_Sigma0 %*% M0 %*% solve_Omega0 + solve_W_list[,,r] %*% temp_zY_sum %*% solve_Sigma2[,,r] )
      # first generate a vector
      solve_temp_zeta <- solve(temp_zeta)
      # print(solve_temp_zeta %*% temp_xi)
      temp_M_vec <- mvtnorm::rmvnorm(1, mean = solve_temp_zeta %*% temp_xi, sigma = solve_temp_zeta)
      #
      M[,,r] <- matrix(temp_M_vec, nrow = p, ncol = q) # convert a vector to matrix by column
      # clusterAssign is a vector, convert it to Z_ij
      # rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
    
    # Update W
    
    W_bar <- matrix(0, nrow=p, ncol=p) + Sigma1[,,] 
      #
    for(ii in 1:n){
      W_bar <- W_bar + (clusterAssign[ii] == r) * ######### tbd
        (data[,,ii] - M[,,r]) %*% solve_Sigma2[,,r] %*% t(data[,,ii] - M[,,r])
      W_list[,,ii] <- MCMCpack::riwish(nu+p+q-1, W_bar)
      solve_W_list[,,ii] = solve(W_list[,,ii])
    }
      #
      # 
    # Update Sigma1
    for(r in 1:nClusters){
      Sigma1_scale <- matrix(0, nrow=q, ncol=q) + gamma*PSI1
      # 
      for(ii in 1:n){
        Sigma1_scale <- Sigma1_scale + 1 * solve_W_list[,,ii]
        
      }
      # 
      Sigma2[,,r] <- MCMCpack::riwish(2*psi+n_Cluster_size[r]*p, Sigma2_scale)
      # fix the trace of V to be q
      Sigma2[,,r] <- q * Sigma2[,,r] / sum(diag(Sigma2[,,r]))
      solve_Sigma2[,,r] = solve(Sigma2[,,r])
      # 
    }
    
    # Update Sigma2
    for(r in 1:nClusters){
      Sigma2_scale <- matrix(0, nrow=q, ncol=q) + gamma*PSI2
      # 
      for(ii in 1:n){
        Sigma2_scale <- Sigma2_scale + 1 * t(data[,,ii] - M[,,clusterAssign[ii]]) %*% solve_W %*% (data[,,ii] - M[,,clusterAssign[ii]])
        
      }
      # 
      Sigma2[,,r] <- MCMCpack::riwish(2*psi+n_Cluster_size[r]*p, Sigma2_scale)
      # fix the trace of V to be q
      Sigma2[,,r] <- q * Sigma2[,,r] / sum(diag(Sigma2[,,r]))
      solve_Sigma2[,,r] = solve(Sigma2[,,r])
      # 
    }
    
    # History[[iter]] <- list(zout = clusterAssign, phiout = phi)
    History[[iter]] <- list(zout = clusterAssign, Mout = M, Wout = W, Sigma2out = Sigma2)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

##################################################
################## E N D ############################

