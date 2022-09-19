

#' An alternative version to GP.simulate.curve.fast in BayesGPfit
#'
#' This function is used to generate basis functions and eigenvalues
#'
#' @param x A p by d matrix, denoting the grids of a d-dim image with p pixels.
#' This is usually produced by BayesGPfit::GP.generate.grids
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
    cb = 0
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
    M = alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta + matrix(rnorm(n.sample * P, sd = sd.noise.a),
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
    M.reduced = M - (alpha_test_ST%*%t(X) + rep(1,P)%*% t(zetam)%*%C + eta)
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
    return(list(Y = Y, M = M,  cb=cb, beta=true.image.b, alpha = true.image.a, X=X,C=C,
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

#' Bayesian Image Mediation Analysis
#'
#' BIMA function is an illustrative example for working with small data sets.
#' For large scale imaging mediation analysis, please run Y_regression_region_block_fast and
#' M_regression_region_block separately for the scalar-on-image and image-on-scalar regressions.
#' These 2 main functions need to be tuned separately to achive best performance.
#' @param Y The scalar outcome, n by 1
#' @param M The image predictor, n by p
#' @param X The scalar exposure variable, n by 1
#' @param C The q confounders, q by n
#' @param dim An integer to specify the dimension of the input image, where it is 1D 2D or 3D image.
#' This must be specified when kernel_setting = "Self-defined".
#' @param grids A matrix of size p by dim. For a 2D image,
#' \texttt{grids[j,]} represents the location (x-y coordinate) of pixel j. For 3D image, \texttt{grids[j,,]} represents
#' the location (x-y-z coordinates) of voxel j. This must be specified when kernel_setting = "Self-defined".
#' @param init_y A list of initial parameter settings for Scalar-on-image regression.
#' \itemize{
#'      \item theta_beta, A vector of length L, where L is the total number of basis functions. Default value is rep(1,L).
#'      \item a_sigma_beta Scalar, parameter a in IG(a,b) for the prior of sigma_Y and sigma_beta. Default value is 1
#'      \item b_sigma_beta Scalar, parameter b in IG(a,b) for the prior of sigma_Y and sigma_beta. Default value is 1
#'      \item sigma_Y Scalar,initial value for sigma_Y. Default value is 1.
#'      \item sigma_beta Scalar,initial value for sigma_beta. Default value is 1.
#'      \item cb Scalar, initial value for the intercept term. Default value is 0.
#'      \item zetay Vector of length q. Default value is rep(1,q).
#'      \item gamma Scalar, initial value for the coefficient of exposure variable X. Default value is 1.
#' }
#' @param init_m A list of intial parameter settings for Image-on-scalar regression.
#' \itemize{
#'     \item theta_alpha A vector of length L, where L is the total number of basis functions. Default value is rep(1,L).
#'     \item theta_eta A matrix of size L by n. Initial value for the inidividual effect eta_i. Default value is matrix(0, L,n)
#'     \item sigma_M Scalar, initial value for sigma_M. Default value is 1.
#'     \item sigma_alpha Scalar, initial value for sigma_alpha Default value is 1.
#'     \item sigma_eta Scalar, initial value for sigma_eta Default value is 0.1.
#'     \item theta_zetam A matrix of size L by m, initial value for the basis coefficients for the coefficient zetam for confounders C.
#'     Default value is  matrix(rep(0.1,L*m),L,m)).
#' }
#' @param controls_y A list of controls for running MALA algorithm for the scalar-on=image regression.
#' \itemize{
#'     \item lambda The thresholding parameter in STGP prior. Default value is 0.5
#'     \item n_mcmc Number of total MCMC iterations, default  = 1e5.
#'     \item stop_adjust The first number of iterations where step_size will be adjusted to achieve target acceptance rate,
#'     default = 0.8*1e5.
#'     \item start_joint The algorithm will update theta_beta alone for the first \texttt{start_joint} number of iterations.
#'      After that all parameters will be jointly updated. We only recommend setting start_joint>0 when the problem is very high-dimensional.
#'      default = 0.
#'     \item interval_step An integer, set the frequency of adjusting step size, default = 10
#'     \item interval_thin An integer, set the frequency of saving the MCMC for theta_beta, default = 1.
#'     \item step A numeric variable, initial step size. Default = 1e-2/n
#' }
#' @param controls_m A list of controls for running MALA algorithm for the image-on-scalar regression.
#' \itemize{
#'     \item lambda The thresholding parameter in STGP prior. Default value is 0.5
#'     \item n_mcmc Number of total MCMC iterations, default  = 2e4.
#'     \item stop_adjust The first number of iterations where step_size will be adjusted to achieve target acceptance rate,
#'     default = 0.8*2e4.
#'     \item start_joint start_joint The algorithm will update theta_alpha alone for the first \texttt{start_joint} number of iterations.
#'      After that all parameters will be jointly updated. We only recommend setting start_joint>0 when the problem is very high-dimensional.
#'      default = 0.
#'     \item interval_step An integer, set the frequency of adjusting step size, default = 10
#'     \item interval_eta An integer, set the frequency of updating the individual effect theta_eta.
#'     The accuracy of theta_eta does not effect the result too much. Larger interval_eta can potentially improve the computational speed.
#'     default = 1000.
#'     \item thinning An integer, set the frequency of saving theta_alpha, default = 1
#'     \item step A numeric variable, initial step size. Default = 1e-2/n
#' }
#' @param region_idx_list A list of length num_region. Users are allowed to pre-specify
#'  the region parcellation by putting indices for each region as each component in this list.
#'  The dafault value is list(region1 = 1:p) with only 1 region.
#' @param kernel_setting A list to specify the kernel setting by the following components
#' \itemize{
#'     \item method A character string to specify the type of kernel used. Currently implemented for
#'     Exponential_square for the modified exponential square kernel and Matern for matern kernel.
#'     Users can also set \texttt{method = "Self-defined"}, and directly use the basis functions and eigenvalues
#'     by putting them in kernel_setting$Phi_Q and kernel_setting$Phi_D respectively.
#'
#'     \item kernel_params A list object.
#'     If \texttt{method = "Exponential_square"}, the default value is \texttt{kernel_params = list(a = 0.01,b=10,poly_degree = 10)},
#'     where a and b are the parameters in modified exponential square kernel \sqn{K(s,t) = \exp(-a(s^2+t^2) - b(s-t)^2)}.
#'     poly_degree is an integer number specifying the highest degree of Hermite polynomials.
#'     If \texttt{method = "Matern"}, the default value is \texttt{kernel_params = list(L_vec = 0.5*unlist(lapply(region_idx,length)), scale = 2, nu = 1/5)}.
#'     L_vec is a vector of integers specifying the number of basis functions used in each region.
#'     See wikipedia for the definition of matern kernel.
#'     Scale represents \sqn{\rho} in this definition, and nu represents \sqn{\nu}.
#'
#'     \item Phi_Q A list object with length num_region.
#'     Each component represents one basis function, a matrix of dimension p_r by L_r for the r-th region.
#'     Needs to be specified if method = "Self-defined".
#'     \item Phi_D A list object with length num_region.
#'     Each component represents one set of eigenvalues for one region, a vector of length L_r for the r-th region.
#'     Needs to be specified if method = "Self-defined".
#' }
#' @return A list of MCMC samples. The default will output the last 20\% as the samples after burn-in.
#' \itemize{
#'    \item beta_sample
#'    \item alpha_sample
#'    \item NIE_sample MCMC sample for the natural indirect effect, product of alpha*beta
#'    \item gamma_sample MCMC sample for the natural direct effect
#' }
#' @importFrom BayesGPfit GP.generate.grids
#' @export
#' @examples
#' if(FALSE){
#'     data(one_region)
#‘     BIMA_mcmc = BIMA(one_region$Y, one_region$X, one_region$M, one_region$C)
#'     # see also ./docs/small_example.R
#' }

BIMA = function(Y,X,M,C, dim = 2, grids = NULL,
                init_y = NULL,
                init_m = NULL,
                controls_y = NULL,
                controls_m = NULL,
                region_idx_list = NULL,
                kernel_setting = list(method = "Exponential_square",
                                      kernel_params = NULL,
                                      Phi_Q = NULL, Phi_D = NULL)){
  # begin the function
  n = length(Y)
  p = dim(M)[1]
  q = dim(C)[1]


  # step1: generate kernel and grids
  if(is.null(region_idx_list)){
    region_idx_list = list(region1 = 1:p)
  }

  if(kernel_setting$method == "Exponential_square"){
    if(is.null(grids)){grids = GP.generate.grids(d=dim,num_grids=sqrt(p))}
    if(is.null(kernel_setting$kernel_params)){
      kernel_setting$kernel_params = list(a = 0.01,b=10,poly_degree = 10)
    }
    basis = generate_sq_basis(grids, region_idx_list,a = kernel_setting$kernel_params$a,
                              b = kernel_setting$kernel_params$b,
                              poly_degree = kernel_setting$kernel_params$poly_degree)
    kernel_setting$Phi_Q = basis$Phi_Q
    kernel_setting$Phi_D = basis$Phi_D
  }else if(kernel_setting$method == "Matern"){
    if(is.null(grids)){grids = GP.generate.grids(d=dim,num_grids=sqrt(p))}
    L_vec = 0.5*unlist(lapply(region_idx_list,length))
    if(is.null(kernel_setting$kernel_params)){
      kernel_setting$kernel_params = list(L_vec = L_vec, scale = 2, nu = 1/5)
    }
    basis = generate_matern_basis2(grids, region_idx_list,kernel_setting$kernel_params$L_vec,
                                   scale = kernel_setting$kernel_params$scale,
                                   nu = kernel_setting$kernel_params$nu)
    kernel_setting$Phi_Q = basis$Phi_Q
    kernel_setting$Phi_D = basis$Phi_D
  }else if(kernel_setting$method == "Self-defined"){
      if(is.null(kernel_setting$Phi_Q) | is.null(kernel_setting$Phi_D)){
        stop("kernel_setting$Phi_Q and kernel_setting$Phi_D cannot be null when using Self-defined kernels.")
      }

  }else{
      stop("kernel_setting$method: only Exponential_square and Matern are supported.
        Users are encouraged to use their own kernels by specifying Phi_Q and Phi_D in kernel_setting. ")
  }

  # step2: set initial values
  L_all = unlist(lapply(basis$Phi_D,length))
  num_region = length(basis$Phi_D)
  L = sum(L_all)

  if(is.null(init_m)){
    init_m = list(theta_alpha =  rep(1,L),
                  theta_eta = matrix(0,L,n),
                  sigma_M = 1,
                  sigma_alpha = 1,
                  sigma_eta = 0.1,
                  theta_zetam = matrix(rep(0.1,L*q),L,q))
  }
  if(is.null(init_y)){
    init_y = list(theta_beta = rep(1,L),
         a_sigma_beta = 1, b_sigma_beta = 1,
         sigma_Y = 1,
         sigma_beta = 1,
         cb = 0 ,
         zetay = rep(1,q),
         gamma = 1)
  }
  if(is.null(controls_y)){
    controls_y = list(lambda = 0.5,
                      n_mcmc = 1e5,
                      stop_adjust = 0.8*1e5,
                      start_joint = 0,
                      interval_step = 10,
                      interval_thin = 1,
                      target_accept_vec = rep(0.2,num_region),
                      step = 1e-2/n)
  }
  if(is.null(controls_m)){
    controls_m = list(lambda = 0.5,
                      n_mcmc = 2e4,
                      stop_adjust = 0.8*2e4,
                      start_joint = 0,
                      interval_step = 10,
                      interval_eta = 1000,
                      thinning = 1,
                      target_accept_vec = rep(0.2,num_region),
                      step = 1e-2/n)
  }

  # set up other parameters
  init_m$D = init_y$D = unlist(basis$Phi_D)




  region_idx_cpp = lapply(region_idx_list, function(x){x-1})
  print("Running scalar-on-image regression ....")
  sim64y = Y_regression_region_block_fast(Y = Y, M = M,
                                          X = X, C = t(C),
                                          L_all = L_all,
                                          num_region = num_region,
                                          region_idx = region_idx_cpp,
                                          n_mcmc = controls_y$n_mcmc,
                                          basis$Phi_Q,
                                          stop_burnin = controls_y$stop_adjust,
                                          start_joint = controls_y$start_joint,
                                          lambda = controls_y$lambda,
                                          target_accept_vec = controls_y$target_accept_vec,
                                          a=init_y$a_sigma_beta,b=init_y$b_sigma_beta,
                                          init = init_y,
                                          step = controls_y$step,
                                          interval_step = controls_y$interval_step,
                                          interval_thin = controls_y$interval_thin)
  print("scalar-on-image regression completed!")
  # note that in GS version, Image on scalar regression, zeta_m is implemented as a scalar instead of spatially-varying vector

  print("Running image-on-scalar regression ....")
  sim64m = M_regression_region_block(M,X, t(C),
                                     L_all,
                                     num_region = num_region ,
                                     region_idx = region_idx_cpp,
                                     n_mcmc = controls_m$n_mcmc ,
                                     basis$Phi_Q,
                                     stop_burnin = controls_m$stop_adjust,
                                     lambda = controls_m$lambda,
                                     target_accept_vec = controls_m$target_accept_vec,
                                     init = init_m,
                                     interval = controls_m$interval_step, # adjust step every 10 iter
                                     interval_eta = controls_m$interval_eta,
                                     thinning = controls_m$thinning,
                                     step = controls_m$step)
  print("image-on-scalar regression completed!")

  print("summarizing results....")
  n_mcmc = dim(sim64y$theta_beta_mcmc_thin)[2]
  theta_sample = sim64y$theta_beta_mcmc_thin[, floor(controls_y$n_mcmc*0.8/controls_y$interval_thin):floor(controls_y$n_mcmc/controls_y$interval_thin)]
  beta_sample = STGP_mcmc(theta_sample,region_idx_list,basis,lambda = controls_y$lambda)
  n_mcmc = dim(sim64m$theta_alpha_mcmc)[2]
  theta_sample = sim64m$theta_alpha_mcmc[,floor(controls_m$n_mcmc*0.8/controls_m$thinning):floor(controls_m$n_mcmc/controls_m$thinning)]
  alpha_sample = STGP_mcmc(theta_sample,region_idx_list,basis,lambda = controls_m$lambda)

  beta_sample_thin = beta_sample[,seq(1,dim(beta_sample)[2],length.out = dim(alpha_sample)[2])]
  TIE_sample = beta_sample_thin*alpha_sample

  output = list(beta_sample = beta_sample_thin, alpha_sample = alpha_sample, TIE_sample = TIE_sample,
    gamma_sample = sim64y$gamma_mcmc[ceiling(n_mcmc*0.8):n_mcmc], grids = grids, kernel_setting = kernel_setting)

  return(output)

}

