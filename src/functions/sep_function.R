#Implementation of SEP from Hutson et al. (2019)
library(gsl)
library(fGarch)
library(stats)

#helper functions
kfun<-function(beta,alpha){
  kinv<-gamma(1+((beta+1)/2))*2^(1+0.5*(1+beta))/(4*alpha*(1-alpha))
  k<-1/kinv
  return(k)
}

zfun<-function(x,theta,sigma){
  z<-(x-theta)/sigma
  return(z)
}

#SEP pdf function
sep_pdf<-function(x,theta,sigma,beta,alpha){
  k<-kfun(beta,alpha)
  z<-zfun(x,theta,sigma)
  y=(k/sigma)*exp(-0.5*(abs(z)+(2*alpha-1)*z)^(2/(1+beta)))
  return(y)
}

#SEP cdf function
sep_cdf<-function(x,theta,sigma,beta,alpha){
  z<-zfun(x,theta,sigma)
  y<-c()
  for(i in 1:length(z)){
    if(z[i]<=0){
      gam_x<-2^((2/(beta+1))-1) * ((alpha-1)*z[i])^(2/(beta+1))
      gam_a<-((beta+1)/2)
      num<-alpha*gamma_inc(gam_a,gam_x)
      den<-gamma((beta+1)/2)
      y[i]<-num/den
    }
    if(z[i]>0){
      gam_x<-2^((2/(beta+1))-1)*(1/(alpha*z[i]))^(-2/(beta+1))
      gam_a<-((beta+1)/2)
      num<-(1-alpha)*gamma_inc(gam_a,gam_x,)
      den<-gamma((beta+1)/2)
      y[i]<-1-(num/den)
    }
  }
  return(y)
}

#SEP quantile function
sep_quant<-function(u,theta,sigma,beta,alpha){
  q<-c()
  for(i in 1:length(u)){
    if(u[i]<=alpha){
      num<-(2^(0.5-(1/(beta+1))) * sqrt(qgamma((1-(u[i]/alpha)),((beta+1)/2))))^(beta+1)
      den<-1-alpha
      q[i]<-theta - sigma*(num/den)
    }
    if(u[i]>alpha){
      num<-(2^((1/(beta+1))-0.5) / sqrt(qgamma((1-((1-u[i])/(1-alpha))),((beta+1)/2))))^(-(beta+1))
      den<-alpha
      q[i]<-theta + sigma*(num/den)
    }
  }
  return(q)
}

sep_quant2<-function(u,theta,sigma,beta,alpha){
  q<-c()
  for(i in 1:length(u)){
    if(u[i]<=alpha){
      num<-(2^(0.5-(1/(beta+1))) * sqrt(gamma((beta+1)/2)/gamma_inc((beta+1)/2,u[i]/alpha)))^(beta+1)
      den<-1-alpha
      q[i]<-theta - sigma*(num/den)
    }
    if(u[i]>alpha){
      num<-(2^((1/(beta+1))-0.5) / sqrt(gamma((beta+1)/2)/gamma_inc((beta+1)/2,(1-u[i])/(1-alpha))))^(-(beta+1))
      den<-alpha
      q[i]<-theta + sigma*(num/den)
    }
  }
  return(q)
}

#SEP sampling function
sep_samp<-function(n,theta,sigma,beta,alpha){
  u<-runif(n)
  samps<-sep_quant(u,theta,sigma,beta,alpha)
  return(samps)
}


#MLE functions

#The functions within ----- are from Hutson et al., but did not work for some reason
#----------------------------------------------------------------------------
#log likelihood function per observation
#sep_ll<-function(x,theta,sigma,beta,alpha){
  #ll<-(-1)*((abs(x-theta)+(2*alpha-1)*(x-theta))^(2/(beta+1)))/(2*sigma) + log(alpha*(1-alpha)) - 
    #(1+0.5*(beta+1))*log(2) - log(gamma((beta+1)/2+1)) - log(sigma)
  #return(ll)
#}

#log likelihood function for MLE with analytic sigma estimation
#sep_mle<-function(x,pars){
  #mu=pars[1]
  #sigma=pars[2]
  #beta=pars[3]
  #alpha=pars[4]
  #ll=sum(sep_ll(x,mu,sigma,beta,alpha))
  #if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  #return(ll)
#}

#log likelihood function for MLE with analytic sigma estimation
#sep_mle_sigest<-function(x,pars){
  #mu=pars[1]
  #beta=pars[2]
  #alpha=pars[3]
  #sigma=sep_sigma_est(x,mu,beta,alpha)
  #ll=sum(sep_ll(x,mu,sigma,beta,alpha))
  #if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  #return(ll)
#}
#-----------------------------------------------------------------------------

#log likelihood function for MLE using canonical form with SEP pdf
sep_mle<-function(x,pars){
  theta=pars[1]
  sigma=pars[2]
  beta=pars[3]
  alpha=pars[4]
  ll=sum(log(sep_pdf(x,theta,sigma,beta,alpha)))
  if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  return(ll)
}

#analytic sigma estimate
sep_sigma_est<-function(x,theta,beta,alpha){
  n<-length(x)
  num<-sum((abs(x-theta)+(2*alpha-1)*(x-theta))^(2/(beta+1)))
  den<-n*(beta+1)
  sigma<-(num/den)^((beta+1)/2)
  return(sigma)
}

#log likelihood function for MLE using canonical form with SEP pdf and analytic sigma estimate
sep_mle_sigest<-function(x,pars){
  theta=pars[1]
  beta=pars[2]
  alpha=pars[3]
  sigma=sep_sigma_est(x,theta,beta,alpha)
  ll=sum(log(sep_pdf(x,theta,sigma,beta,alpha)))
  if(ll==(-Inf)|ll==Inf|is.na(ll)==T){ll<-(-1e9)}
  return(ll)
}


SEPfit<-function(x){
  #optimization constraints for theta, sigma, beta, alpha params
  par_upr<-c(5,10,10,1)
  par_est<-c(0,1,0,0.5)
  par_lwr<-c(-5,0,-1,0)

  sep_fit<-optim(par=par_est,sep_mle,x=x,
               method = 'L-BFGS-B',lower = par_lwr,upper = par_upr,
               control = list(fnscale=-1,maxit=100000))

  theta<-sep_fit$par[1]
  sigma<-sep_fit$par[2]
  beta<-sep_fit$par[3]
  alpha<-sep_fit$par[4]
  param_out<-c(theta,sigma,beta,alpha)
  names(param_out)<-c('theta','sigma','beta','alpha')
  return(param_out)
}

SEPfit_sigest<-function(x){
  #optimization constraints for theta, sigma, beta, alpha params
  par_upr<-c(5,10,1)
  par_est<-c(0,0,0.5)
  par_lwr<-c(-5,-1,0)
  
  sep_fit<-optim(par=par_est,sep_mle_sigest,x=x,
                 method = 'L-BFGS-B',lower = par_lwr,upper = par_upr,
                 control = list(fnscale=-1,maxit=100000))
  
  theta<-sep_fit$par[1]
  beta<-sep_fit$par[2]
  alpha<-sep_fit$par[3]
  sigma<-sep_sigma_est(x,theta,beta,alpha)
  param_out<-c(theta,sigma,beta,alpha)
  names(param_out)<-c('theta','sigma','beta','alpha')
  return(param_out)
}

forc_err_SEPfit<-function(err_mat,use_sigest,distn,seasonal,start_date,end_date,lv_out_yrs,parallel,use_mpi){
  if(use_mpi==FALSE){
    library(doParallel)
    parallel::detectCores()
    n.cores <- parallel::detectCores()
    my.cluster<-try(getDefaultCluster())
    if(class(my.cluster)=='try-error'|class(my.cluster)=='NULL'){
      my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')}
    print(my.cluster)
    doParallel::registerDoParallel(cl = my.cluster)
    foreach::getDoParRegistered()}
  if(use_mpi==TRUE){
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)}
  
  #date index
  ixx_fit_all <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),by='day'))
  
  wy_fun<-function(date_vec){
    wy_vec <- date_vec$year
    wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
    date_vec_wy <- date_vec
    date_vec_wy$year <- wy_vec
    return(date_vec_wy)
  }
  
  ixx_fit_all_wy <- wy_fun(ixx_fit_all)
  
  trn_idx <- !(ixx_fit_all_wy$year%in%(leave_out_years-1900))
  ixx_fit <- ixx_fit_all[trn_idx] #fit years excluding leave out years
  
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'))
  ixx_obs <- as.POSIXlt(seq(as.Date(ixx_obs[1]),as.Date(ixx_obs[length(ixx_obs)]),by='day'))
  ixx_hefs <- as.POSIXlt(seq(as.Date(ixx_hefs[1]),as.Date(ixx_hefs[length(ixx_hefs)]),by='day'))
  
  n_sites = dim(err_mat)[1]
  n_ens = dim(err_mat)[2]
  n_hefs = dim(err_mat)[3]
  n_lds = dim(err_mat)[4]
  
  err_mat_fit<-err_mat[,,ixx_fit_all%in%ixx_fit,,drop=FALSE]
  
  season_list<-seasonal
  
  if(seasonal=='monthly'){
    season_list<-vector('list',12)
    for(i in 1:12){
      season_list[[i]]<-i
    }
  }
  
  if(seasonal=='annual'){
    season_list<-list(1:12)
  }
  
  sep_params <- array(NA,c(n_sites,n_ens,length(season_list),n_lds,4))
  
  if(distn=='sep'){
  if(parallel=='FALSE'){
    for(j in 1:n_sites){
      for(e in 1:n_ens){
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          for(k in 1:n_lds){
            if(use_sigest=='FALSE'){
              sep_params[j,e,i,k,]<-SEPfit(err_mat_fit[j,e,seas_fit,k])}
            if(use_sigest=='TRUE'){
              sep_params[j,e,i,k,]<-SEPfit_sigest(err_mat_fit[j,e,seas_fit,k])}
          }
        }
      }
    }
  }

  if(parallel=='TRUE'){
    for(j in 1:n_sites){
      out<-foreach(e = 1:n_ens,.combine='c',.export=c('SEPfit','SEPfit_sigest','sep_mle','sep_mle_sigest','sep_pdf','sep_sigma_est','kfun','zfun')) %dopar% {
        sep_par<-array(NA,c(dim(sep_params)[3],dim(sep_params)[4],dim(sep_params)[5]))
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          for(k in 1:n_lds){
            if(use_sigest=='FALSE'){
              sep_par[i,k,]<-SEPfit(err_mat_fit[j,e,seas_fit,k])}
            if(use_sigest=='TRUE'){
              sep_par[i,k,]<-SEPfit_sigest(err_mat_fit[j,e,seas_fit,k])}
          }
        }
        return(list(sep_par))
      }
  
      for(e in 1:dim(err_mat)[1]){
        sep_params[j,e,,,]<-out[[e]]
      }
    }
  }
  }
  
  if(distn=='sged'){
    
  if(parallel=='FALSE'){
    for(j in 1:n_sites){
      for(e in 1:n_ens){
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          for(k in 1:n_lds){
            sep_params[j,e,i,k,]<-sgedFit(err_mat_fit[j,e,seas_fit,k])$par
          }
        }
      }
    }
  }
  
  if(parallel=='TRUE'){
    for(j in 1:n_sites){
      out<-foreach(e = 1:n_ens,.combine='c',.packages=c('fGarch')) %dopar% {
        sep_par<-array(NA,c(dim(sep_params)[3],dim(sep_params)[4],dim(sep_params)[5]))
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          for(k in 1:n_lds){
            if(use_sigest=='FALSE'){
              sep_par[i,k,]<-sgedFit(err_mat_fit[j,e,seas_fit,k])$par}
            if(use_sigest=='TRUE'){
              sep_par[i,k,]<-sgedFit(err_mat_fit[j,e,seas_fit,k])$par}
          }
        }
        return(list(sep_par))
      }
    
      for(e in 1:n_ens){
      sep_params[j,e,,,]<-out[[e]]
      }
    }
  }
  }
  
  return(sep_params)
  
  if(use_mpi==FALSE){
    stopCluster(cl = my.cluster)}
  if(use_mpi==TRUE){
    closeCluster(cl)
    Rmpi::mpi.quit()}
}

######################################END##################################