
#' An alternative version to GP.simulate.curve.fast in BayesGPfit
#'
#' This function is used to generate basis functions and eigenvalues
#'
#' @param x A p by d matrix, denoting the grids of a d-dim image with p pixels.
#' This is usually produced by BayesGP::GP.generate.grids
#' @param poly_degree Integer, polynomial degrees.
#' @param a Parameter a in GP kernel \eqn{K(s,t) = \exp(-a^2*(s^2+t^2)-b*(s-t)^2 }
#' @param b Parameter b in GP kernel \eqn{K(s,t) = \exp(-a^2*(s^2+t^2)-b*(s-t)^2 }
#' @param center Center of the grid
#' @param scale Parameter for scaling the distance between grid points
#' @import BayesGPfit
#' @return a list
#' \itemize{
#' \item
#' }
  GP.simulate.curve.fast.new = function(x,poly_degree,a,b,
                                        center=NULL,scale=NULL,max_range=6){

    x = cbind(x)
    d = ncol(x)

    if(is.null(center)){
      center = apply(x,2,mean)
    }
    c_grids = t(x) - center
    if(is.null(scale)){
      max_grids =pmax(apply(c_grids,1,max),-apply(c_grids,1,min))
      scale=as.numeric(max_grids/max_range)
    }

    work_x = GP.std.grids(x,center=center,scale=scale,max_range=max_range)
    Xmat = GP.eigen.funcs.fast(grids=work_x,
                               poly_degree =poly_degree,
                               a =a ,b=b)
    lambda = GP.eigen.value(poly_degree=poly_degree,a=a,b=b,d=d)
    betacoef = rnorm(ncol(Xmat),mean=0,sd=sqrt(lambda))
    f = Xmat%*%betacoef
    return(list(f=f,x=x,work_x=work_x, eigen.func = Xmat, eigen.value = lambda))
  }


  #‘ A function used to generate a testing case, with region separation
  #' @export
  STGP_generate_theta_block = function(true.image.a,  true.image.b,
                                       sd.noise.a,   sd.noise.b, grids,Q,lambda,region_idx,L_all,
                                       n.sample){
    # set true params here:
    n_C = 2
    ca = 0; cb = 0
    zetay = runif(n_C,-2,2); zetam = runif(n_C,-2,2)/1e2
    C = matrix(rnorm(n_C * n.sample), ncol = n.sample)*5
    gamma = 0.01
    X = rnorm(n.sample)*5
    true.image.b = true.image.b*1.5

    p = P= length(true.image.a)
    L = dim(Q[[1]])[2]
    num_block = length(Q)

    theta_eta = vector(mode = "list", length = num_block)
    eta = matrix(NA, nrow=p,ncol=n.sample)
    beta_test = rep(0,p)
    alpha_test = rep(0,p)
    true_theta_alpha = rep(NA, sum(L_all))
    true_theta_beta = rep(NA, sum(L_all))
    sigma_eta=0.1
    theta_eta_final = NULL
    for(m in 1:num_block){
      print(paste("m=",m))
      L = L_all[m]
      L_end = cumsum(L_all)[m]
      theta_eta[[m]] = matrix(rnorm(n.sample*L)*sigma_eta,nrow=L, ncol = n.sample) #L by n
      idx = region_idx[[m]]
      print(paste("dim(Q[[m]])=",dim(Q[[m]]),";dim(theta_eta[[m]])=",dim(theta_eta[[m]])))
      eta[idx,] = Q[[m]]%*%theta_eta[[m]] # p_block by n
      theta_eta_final = rbind(theta_eta_final,theta_eta[[m]])
      true_theta_alpha[(L_end-L+1):L_end] = t(true.image.a[idx]%*%Q[[m]])
      true_theta_beta[(L_end-L+1):L_end] = t(true.image.b[idx]%*%Q[[m]])
      beta_test[idx] = Q[[m]]%*%true_theta_beta[(L_end-L+1):L_end]
      alpha_test[idx] = Q[[m]]%*%true_theta_alpha[(L_end-L+1):L_end]
    }
    beta_test_ST = (beta_test - sign(beta_test)*lambda)*(abs(beta_test)>lambda)
    alpha_test_ST = (alpha_test - sign(alpha_test)*lambda)*(abs(alpha_test)>lambda)
    M = ca + alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta + matrix(rnorm(n.sample * P, sd = sd.noise.a),
                                                                            ncol = n.sample)
    Y = cb + t(M )%*% beta_test_ST + gamma*X + t(C)%*%(zetay) + rnorm(n.sample,mean=0, sd = sd.noise.b)

    # R^2
    # R^2 for M at each location: higher R^2 better accuracy
    snratio = NULL
    M.nonzero = t(M)[,beta_test_ST!=0]
    snratio$R2.beta = 1-sum( (Y - cb- gamma*X - t(C)%*%zetay - M.nonzero%*%solve(crossprod(M.nonzero),t(M.nonzero)%*%Y) )^2 )/sum((Y-mean(Y))^2)
    # resY = Y-(cb + t(M )%*% beta_test_ST + gamma*X + t(C)%*%(zetay))
    snratio$sgn.beta = sd( t(M )%*% beta_test_ST )/sd( Y-(cb  + gamma*X + t(C)%*%(zetay)) )
    snratio$gamma = sd(gamma*X)/sd( Y-(cb + t(M )%*% beta_test_ST + t(C)%*%(zetay)) )
    snratio$zetay = sd(t(C)%*%zetay)/sd( Y-(cb + t(M )%*% beta_test_ST + gamma*X )  )
    snratio$Ymodel = sd( t(M )%*% beta_test_ST + gamma*X )/sd(Y)


    # M-regression(check alpha)
    M.reduced = M - (ca + alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta)
    signal = alpha_test_ST%*%t(X)
    snratio$sgn.alpha = rep(NA, P)
    for(j in 1:P){
      M_j = M.reduced[j,]
      f_j = signal[j,]
      snratio$sgn.alpha[j] = sd(f_j)/sd(M_j)
    }

    true.ll = NULL
    true.ll$beta = -sum((Y-(cb + t(M )%*% beta_test_ST + gamma*X + t(C)%*%zetay))^2)/2/sd.noise.b^2
    logll_M = 0
    for(m in 1:num_block){
      Q_t = t(Q[[m]])
      idx = region_idx[[m]]
      p_i = length(idx)
      L = L_all[m]
      L_end = cumsum(L_all)[m]
      L_idx = (L_end-L+1):L_end
      M_star_m = Q_t%*%(M[idx,] - as.matrix(rep(1,p_i))%*%t(t(C)%*%zetam)-alpha_test_ST[idx]%*%t(X)) - theta_eta_final[L_idx,];
      logll_M = logll_M -0.5*(norm(M_star_m,"f"))^2/sd.noise.a^2
    }

    norm_M = norm( M -( alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C)-eta, type="f")
    true.ll$alpha = logll_M
    return(list(Y = Y, M = M, ca = ca, cb=cb, beta=true.image.b, alpha = true.image.a, X=X,C=C,
                zetay = zetay, zetam = zetam,gamma = gamma,theta_eta=theta_eta_final,
                theta_alpha = true_theta_alpha, theta_beta = true_theta_beta,
                sigma_y = sd.noise.b, sigma_m = sd.noise.a,snratio=snratio,
                true.ll = true.ll, beta_test_ST=beta_test_ST,alpha_test_ST=alpha_test_ST,sigma_eta=sigma_eta ))
  }

  FDR = function(active_region, true_region){
    sum(active_region!=0 & true_region==0)/sum(active_region!=0)
  }
  Precision = function(active_region, true_region){
    mean(I(active_region!=0) == I(true_region!=0))
  }
  Power = function(active_region, true_region){
    sum(active_region !=0 & true_region!=0)/sum(true_region!=0)
  }

  matern_kernel = function(x,y,nu,l=1){
    d = sqrt(sum((x-y)^2))/l
    y = 2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d)^nu*besselK(sqrt(2*nu)*d,nu)
    return(y)
  }
  High_to_low = function(beta_mcmc_sample,region_idx,basis){
    n = dim(beta_mcmc_sample)[2]
    p = max(unlist(region_idx))
    num_region = length(region_idx)
    L = sum(unlist(lapply(basis$Phi_D,length)))
    est_mcmc = matrix(NA,nrow = L, ncol = n)
    dd = matrix(unlist(lapply(basis$Phi_Q,dim)),nrow=2)
    L_idx = cumsum(dd[2,])
    L_idx = c(rbind(L_idx,L_idx+1))
    L_idx = matrix(c(1,L_idx[-length(L_idx)]),nrow=2)
    for(l in 1:num_region){
      idx = region_idx[[l]]
      beta = beta_mcmc_sample[idx,]
      theta = t(basis$Phi_Q[[l]])%*%beta
      est_mcmc[L_idx[1,l]:L_idx[2,l],] = theta
    }
    return(est_mcmc)
  }

#‘ Generate a matern basis
#' @importFrom RSpectra eigs_sym
#' @export
  generate_matern_basis2 = function(grids, region_idx_list, L_vec,scale = 2,nu = 1/5,
                                    show_progress = FALSE){
    if(nu=="vec"){
      nu_vec = region_idx_list["nu_vec"]
    }
    num_block = length(region_idx_list)
    Phi_D = vector("list",num_block)
    Phi_Q = vector("list",num_block)
    Lt = NULL; pt = NULL
    for(i in 1:num_block){
      if(show_progress){
        print(paste("Computing basis for block ",i))
      }
      p_i = length(region_idx_list[[i]])
      kernel_mat = matrix(NA,nrow = p_i, ncol=p_i)
      for(l in 1:p_i){
        if(nu=="vec"){
          kernel_mat[l,] = apply(grids[region_idx_list[[i]],],1,matern_kernel,y=grids[region_idx_list[[i]],][l,],nu = nu_vec[i],l=scale)
        }else{
          kernel_mat[l,] = apply(grids[region_idx_list[[i]],],1,matern_kernel,y=grids[region_idx_list[[i]],][l,],nu = nu,l=scale)
        }
      }
      diag(kernel_mat) = 1
      K = eigs_sym(kernel_mat,L_vec[i])
      K_QR = qr(K$vectors)
      Phi_Q[[i]] = qr.Q(K_QR )
      Phi_D[[i]] = K$values
      Lt = c(Lt, length(Phi_D[[i]]))
      pt = c(pt, dim(Phi_Q[[i]])[1])
    }
    return(list(Phi_D = Phi_D,
                region_idx_block = region_idx_list,
                Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
  }
  #‘ Generate modified exponential square basis
  #' @importFrom RSpectra eigs_sym
  #' @export
  generate_sq_basis = function(grids, region_idx_list,poly_degree_vec,a = 0.01, b=10, poly_degree=20,
                               show_progress=FALSE){
    num_block = length(region_idx_list)
    Phi_D = vector("list",num_block)
    Phi_Q = vector("list",num_block)
    Lt = NULL; pt = NULL
    for(i in 1:num_block){
      if(show_progress){
        print(paste("Computing basis for block ",i))
      }
      GP = GP.simulate.curve.fast.new(x=grids[region_idx_list[[i]],], a=a ,b=b,poly_degree=poly_degree) # try to tune b, increase for better FDR
      K_esq = GP$eigen.func
      K_QR = qr(K_esq)
      Phi_Q[[i]] = qr.Q(K_QR)
      Phi_D[[i]] = GP$eigen.value
      Lt = c(Lt, length(Phi_D[[i]]))
      pt = c(pt, dim(Phi_Q[[i]])[1])
    }
    return(list(Phi_D = Phi_D,
                region_idx_block = region_idx_list,
                Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
  }

  generate_sq_basis_vec = function(grids, region_idx_list,poly_degree_vec,a_vec,b_vec, show_progress=FALSE){
    num_block = length(region_idx_list)
    Phi_D = vector("list",num_block)
    Phi_Q = vector("list",num_block)
    Lt = NULL; pt = NULL
    for(i in 1:num_block){
      if(show_progress){
        print(paste("Computing basis for block ",i))
      }
      GP = GP.simulate.curve.fast.new(x=grids[region_idx_list[[i]],], a=a_vec[i] ,b=b_vec[i],poly_degree=poly_degree_vec[i]) # try to tune b, increase for better FDR
      K_esq = GP$eigen.func
      K_QR = qr(K_esq)
      Phi_Q[[i]] = qr.Q(K_QR)
      Phi_D[[i]] = GP$eigen.value
      Lt = c(Lt, length(Phi_D[[i]]))
      pt = c(pt, dim(Phi_Q[[i]])[1])
    }
    return(list(Phi_D = Phi_D,
                region_idx_block = region_idx_list,
                Phi_Q = Phi_Q,L_all = Lt,p_length=pt))
  }

  MSE_func = function(mcmc_sample,true_beta,inclusion_map){
    S = dim(mcmc_sample)[1]
    MSE1 = NULL
    MSE0 = NULL
    for(s in 1:S){
      if(inclusion_map[s]==1){
        MSE1 = c(MSE1, mean((mcmc_sample[s,]-true_beta[s])^2) )
      }
      if(inclusion_map[s]==0){
        MSE0 = c(MSE0, mean((mcmc_sample[s,]-true_beta[s])^2) )
      }
    }
    return(list(MSE1 = MSE1, MSE0 = MSE0))
  }

  InclusionMap = function(mcmc_sample, true_beta, thresh = "auto", fdr_target = 0.1,
                          max.iter = 100){
    InclusionProb = 1 - apply(mcmc_sample, 1, function(x){mean(abs(x)==0)})
    true_beta = 1*(true_beta!=0)
    thresh_final = thresh
    fdr=NA
    if(thresh=="auto"){
      thresh = 0.5
      for(i in 1:max.iter){
        mapping = 1*(InclusionProb>thresh)
        fdr = FDR(mapping, true_beta)
        print(paste("fdr=",fdr,"thresh=",thresh))
        if(is.na(fdr)){
          print("fdr=NA, target FDR is too small")
          thresh = thresh/1.1
          mapping = 1*(InclusionProb>thresh)
          fdr = FDR(mapping, true_beta)
          print(paste("Use current fdr=",fdr,"thresh=",thresh))
          break
        }
        if(fdr<=fdr_target){
          thresh_final = thresh
          break
        }
        thresh = thresh*1.1
        if(thresh>1){
          print("New thresh>1, keep thresh at the current value and return result.")
          break
        }
      }
    }else{
      mapping = 1*(InclusionProb>thresh)
    }
    return(list(mapping = mapping, thresh=thresh_final,
                InclusionProb=InclusionProb))
  }

  CoverageProb = function(beta_sample,true_beta, lower = 0.025,upper = 0.975,error = 1e-5){
    CI = apply(beta_sample,1,function(x){quantile(x,probs = c(lower,upper))})
    return(list(cp = mean( true_beta>= CI[1,]-error & true_beta<=CI[2,]+error),
                lowerq = CI[1,], upperq = CI[2,]) )
  }

  STGP_mcmc = function(theta_mcmc_sample,region_idx,basis,lambda){
    M = dim(theta_mcmc_sample)[2]
    S = max(unlist(region_idx))
    num_region = length(region_idx)
    est_mcmc = matrix(NA,nrow = S, ncol = M)
    dd = matrix(unlist(lapply(basis$Phi_Q,dim)),nrow=2)
    L_idx = cumsum(dd[2,])
    L_idx = c(rbind(L_idx,L_idx+1))
    L_idx = matrix(c(1,L_idx[-length(L_idx)]),nrow=2)
    for(l in 1:num_region){
      idx = region_idx[[l]]
      theta = theta_mcmc_sample[L_idx[1,l]:L_idx[2,l],]
      beta = basis$Phi_Q[[l]]%*%theta
      est_mcmc[idx,] = (beta-sign(beta)*lambda)*I(abs(beta)>lambda)
    }
    return(est_mcmc)
  }

  FDR_control = function(mcmc_sample, lambda=0.5,set_FDR = 0.1,delta=0.1){

    # est_centered = apply(est_mcmc,2,function(x){(x-mean(x))/sd(x)+mean(x)})
    M = dim(mcmc_sample)[2]
    S = dim(mcmc_sample)[1]
    P_fdr = apply(mcmc_sample, 1, function(x){mean(abs(x)<=delta)})
    P_fdr[P_fdr==0] = 1/(2*M)
    P_fdr.sorted = sort(P_fdr)
    P_aim = cumsum(P_fdr.sorted)/(1:S)
    idx = which(P_aim>set_FDR)[1]-1
    signal = 1*(P_fdr<=P_aim[idx])
    return(signal)
  }

