
library(BigVAR)
library(vars)
library(bigtime)

spVAR<-bigtime::sparseVAR

lagmat_fun<-function(err_mat,lag){
  lagmat<-array(0,c(dim(err_mat)[1],dim(err_mat)[2]*lag))
  for(i in 1:lag){
    lagmat[,((i-1)*dim(err_mat)[2]+1):(i*dim(err_mat)[2])]<-rbind(matrix(0,nrow=i,ncol=dim(err_mat)[2]),err_mat[1:(dim(err_mat)[1]-i),])
  }
  return(lagmat)
}

var_resids_fun<-function(err_mat,var_coef,lag,mean_vec,use_mean){
  lagmat<-lagmat_fun(err_mat,lag)
  if(use_mean=='FALSE'){
    pred<-lagmat%*%t(var_coef)
    resids<-err_mat - pred
  }
  if(use_mean=='TRUE'){
    lagmat<-cbind(lagmat,rep(1,dim(err_mat)[1]))
    var_coef<-cbind(var_coef,mean_vec)
    pred<-lagmat%*%t(var_coef)
    resids<-err_mat-pred
  }
  return(resids)
}

#function to fit VAR models
var_fit<-function(err_mat,lag,var_mod,use_mean,seasonal,start_date,end_date,lv_out_yrs,parallel,use_mpi){
  ##parallelization code
  if(use_mpi==FALSE){
    library(doParallel)
    parallel::detectCores()
    n.cores <- parallel::detectCores()
    my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')
    print(my.cluster)
    setDefaultCluster(cl = my.cluster)
  
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
  
  #output arrays
  at_mat <- array(NA,dim(err_mat))
  
  if(use_mean=='FALSE'){
    var_coefs <- array(NA,c(n_sites,n_ens,length(season_list),n_lds,n_lds*lag))}
  if(use_mean=='TRUE'){
    var_coefs <- array(NA,c(n_sites,n_ens,length(season_list),n_lds,n_lds*lag+1))}
  
  #1) 'bigtime' sparseVAR model
  if(var_mod=='svar'){
    
  if(parallel=='FALSE'){
    for(j in 1:n_sites){
      for(e in 1:n_ens){
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          var_fit<-sparseVAR(err_mat_fit[j,e,seas_fit,],p=lag,VARpen = 'L1',selection = 'cv',check_std = F)
          var_resids<-var_resids_fun(err_mat[j,e,seas_all,],var_fit$Phihat,lag=lag,var_fit$phi0hat,use_mean)
          at_mat[j,e,seas_all,]<-var_resids
          if(use_mean=='FALSE'){
            var_coefs[j,e,i,,]<-var_fit$Phihat}
          if(use_mean=='TRUE'){
            var_coefs[j,e,i,,]<-cbind(var_fit$Phihat,var_fit$phi0hat)}
          var_coefs[var_coefs>1|var_coefs<(-1)]<-0
        }
      }
    }
  }
  
  if(parallel=='TRUE'){
    for(j in 1:n_sites){
      out<-foreach(e = 1:n_ens,.combine='c',.export = c('spVAR','var_resids_fun','lagmat_fun')) %dopar% {
        var_coef<-array(NA,c(dim(var_coefs)[3],dim(var_coefs)[5],dim(var_coefs)[5]))
        at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          var_fit<-spVAR(err_mat_fit[j,e,seas_fit,],p=lag,VARpen = 'L1',selection = 'cv',check_std = F)
          var_resids<-var_resids_fun(err_mat[j,e,seas_all,],var_fit$Phihat,lag=lag,var_fit$phi0hat,use_mean)
          at[seas_all,]<-var_resids
          if(use_mean=='FALSE'){
            var_coef[i,,]<-var_fit$Phihat}
          if(use_mean=='TRUE'){
            var_coef[i,,]<-cbind(var_fit$Phihat,var_fit$phi0hat)}
        }
        return(list(at,var_coef))
      }
    
      for(e in 1:n_ens){
        at_mat[j,e,,]<-out[[(e*2-1)]]
        var_coefs[j,e,,,]<-out[[e*2]]
      }
    }
    var_coefs[var_coefs>1|var_coefs<(-1)]<-0
  }
  
  return(list(at_mat,var_coefs))
  }

  #2) Fit 'BigVAR' lasso penalized VAR model
  if(var_mod=='bvar'){

  if(parallel=='FALSE'){
    for(j in 1:n_sites){
      for(e in 1:n_ens){
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          var_mod = constructModel(err_mat_fit[j,e,seas_fit,], p = lag, struct = "Basic", gran = c(10, 10),IC = F,
              verbose = F, VARX = list(),separate_lambdas = F, rolling_oos=F,model.controls=list(intercept = use_mean, MN=F))
          var_fit = cv.BigVAR(var_mod)
        
          if(use_mean=='FALSE'){
            var_coefs[j,e,i,,]<-var_fit@betaPred[,2:dim(var_fit@betaPred)[2]]
            at_mat[j,e,seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coefs[j,e,i,,],lag,c(0),use_mean)}
          if(use_mean=='TRUE'){
            var_coefs[j,e,i,,]<-cbind(var_fit@betaPred[,2:dim(var_fit@betaPred)[2]],var_fit@betaPred[,1])
            at_mat[j,e,seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_fit@betaPred[,2:dim(var_fit@betaPred)[2]],lag,var_fit@betaPred[,1],use_mean)}
          #var_coefs[var_coefs>1|var_coefs<(-1)]<-0
        }
      }
    }
  }
    
  
  if(parallel=='TRUE'){
    for(j in 1:n_sites){
      out<-foreach(e = 1:n_ens,.combine='c',.packages = c('BigVAR'),.export=c('var_resids_fun','lagmat_fun')) %dopar% {
        var_coef<-array(NA,c(dim(var_coefs)[3],dim(var_coefs)[4],dim(var_coefs)[5]))
        at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
        for(i in 1:length(season_list)){
          seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
          seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          var_mod = constructModel(err_mat_fit[j,e,seas_fit,], p = lag, struct = "Basic", gran = c(10, 10),IC = F,
                            verbose = F, VARX = list(),separate_lambdas = F, rolling_oos=F,model.controls=list(intercept = use_mean, MN=F))
          var_fit = cv.BigVAR(var_mod)
      
          if(use_mean=='FALSE'){
            var_coef[i,,]<-var_fit@betaPred[,2:dim(var_fit@betaPred)[2]]
            at[seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coef[i,,],lag,c(0),use_mean)}
          if(use_mean=='TRUE'){
            var_coef[i,,]<-cbind(var_fit@betaPred[,2:dim(var_fit@betaPred)[2]],var_fit@betaPred[,1])
            at[seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_fit@betaPred[,2:dim(var_fit@betaPred)[2]],lag,var_fit@betaPred[,1],use_mean)}
        }
      return(list(at,var_coef))
      }
    
      for(e in 1:n_ens){
        at_mat[j,e,,]<-out[[(e*2-1)]]
        var_coefs[j,e,,,]<-out[[e*2]]
      }
    #var_coefs[var_coefs>1|var_coefs<(-1)]<-0
    }
  }
  return(list(at_mat,var_coefs))
  }

  #3) Fit rmgarch 'varx' robust least squares VAR model
  if(var_mod=='varx'){
    
  if(parallel=='FALSE'){
    for(j in 1:n_sites){
      for(e in 1:n_ens){
        var_coef<-array(NA,c(dim(var_coefs)[3],dim(var_coefs)[4],dim(var_coefs)[5]))
        at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
        for(i in 1:length(season_list)){
          seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          #try robust least squares VAR fit
          var_fit<-try(varxfit(err_mat[j,e,seas,],lag,constant=use_mean,robust = T),T)
          at[seas,]<-try(rbind(err_mat[j,e,seas,][1:lag,],var_fit$xresiduals),T)
          var_coef[i,,]<-try(var_fit$Bcoef,T)
          #if it doesnt work, try OLS VAR
          if(class(var_fit)=='try-error'){
            var_fit<-try(varxfit(err_mat[j,e,seas,],lag,constant=use_mean,robust = F),T)
            at[seas,]<-try(rbind(err_mat[j,e,seas,][1:lag,],var_fit$xresiduals),T)
            var_coef[i,,]<-try(var_fit$Bcoef,T)
            print(paste('ens',e,'mth',i,'trying OLS VAR'))
          }
          #as last resort, set VAR residuals to normalized errors and VAR coefficients all to zero
          if(class(var_fit)=='try-error' | anyNA(as.numeric(at[seas,]))==T | anyNA(as.numeric(var_coef[i,,]))==T){
            at[seas,]<-err_mat[j,e,seas,]
            var_coef[i,,]<-array(0,dim(var_coef[i,,]))
            print(paste('ens',e,'mth',i,'all zero VAR'))
          }
        }
        #correct missing or erroneous data as needed
        res_out<-at
        res_out[res_out==Inf|res_out==-Inf]<-0
        at_mat[j,e,,]<-apply(res_out,c(1:2),as.numeric)
        coef_out<-var_coef
        coef_out<-apply(coef_out,c(1:3),as.numeric)
        #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
        coef_out[coef_out>=1]<-0
        coef_out[coef_out<=(-1)]<-0
        var_coefs[j,e,,,]<-coef_out
      }
    }
  }
  
  if(parallel=='TRUE'){
    for(j in 1:n_sites){
      out<-foreach(e = 1:n_ens,.combine='c',.packages = c('rmgarch')) %dopar% {
        var_coef<-array(NA,c(dim(var_coefs)[3],dim(var_coefs)[4],dim(var_coefs)[5]))
        at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
        for(i in 1:length(season_list)){
          seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
          #try robust least squares VAR fit
          var_fit<-try(varxfit(err_mat[j,e,seas,],lag,constant=use_mean,robust = T),T)
          at[seas,]<-try(rbind(err_mat[j,e,seas,][1:lag,],var_fit$xresiduals),T)
          var_coef[i,,]<-try(var_fit$Bcoef,T)
          #if it doesnt work, try OLS VAR
          if(class(var_fit)=='try-error'){
          var_fit<-try(varxfit(err_mat[j,e,seas,],lag,constant=use_mean,robust = F),T)
          at[seas,]<-try(rbind(err_mat[j,e,seas,][1:lag,],var_fit$xresiduals),T)
          var_coef[i,,]<-try(var_fit$Bcoef,T)
          print(paste('ens',e,'mth',i,'trying OLS VAR'))
        }
        #as last resort, set VAR residuals to normalized errors and VAR coefficients all to zero
        if(class(var_fit)=='try-error' | anyNA(as.numeric(at[seas,]))==T | anyNA(as.numeric(var_coef[i,,]))==T){
          at[seas,]<-err_mat[j,e,seas,]
          var_coef[i,,]<-array(0,dim(var_coef[i,,]))
          print(paste('ens',e,'mth',i,'all zero VAR'))
        }
      }
      
      return(list(at,var_coef))
      }
    
      for(e in 1:n_ens){
        res_out<-out[[(e*2-1)]]
        res_out[res_out==Inf|res_out==-Inf]<-0
        at_mat[j,e,,]<-apply(res_out,c(1:2),as.numeric)
        coef_out<-out[[e*2]]
        coef_out<-apply(coef_out,c(1:3),as.numeric)
        #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
        coef_out[coef_out>=1]<-0
        coef_out[coef_out<=(-1)]<-0
        var_coefs[j,e,,,]<-coef_out
      }
    }
  }
  
  return(list(at_mat,var_coefs))
  }
  
  #4) Fit vars 'VAR' least squares VAR model
  if(var_mod=='var'){
    
    if(parallel=='FALSE'){
      for(j in 1:n_sites){
        for(e in 1:n_ens){
          for(i in 1:length(season_list)){
            seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
            seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
            #try robust least squares VAR fit
            if(use_mean==TRUE){
              typ='const'
              var_fit<-VAR(err_mat_fit[j,e,seas_fit,],p=lag,type=typ)
              coefs<-coef(var_fit)
              for(k in 1:length(coefs)){
                var_coefs[j,e,i,k,]<-coefs[[k]][,'Estimate']}
              at_mat[j,e,seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coefs[j,e,i,,-c(dim(var_coefs)[5])],lag,var_coefs[j,e,i,,dim(var_coefs)[5]],use_mean)}
            if(use_mean==FALSE){
              typ='none'
              var_fit<-VAR(err_mat_fit[j,e,seas_fit,],p=lag,type=typ)
              coefs<-coef(var_fit)
              for(k in 1:length(coefs)){
                var_coefs[j,e,i,k,]<-coefs[[k]][,'Estimate']}
              at_mat[j,e,seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coefs[j,e,i,,],lag,c(0),use_mean)}
          }
        #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
        var_coefs[var_coefs>=1]<-0
        var_coefs[var_coefs<=(-1)]<-0
        }
      }
    }
    
    if(parallel=='TRUE'){
      for(j in 1:n_sites){
        out<-foreach(e = 1:n_ens,.combine='c',.packages = c('vars'),.export=c('var_resids_fun','lagmat_fun')) %dopar% {
          var_coef<-array(NA,c(dim(var_coefs)[3],dim(var_coefs)[4],dim(var_coefs)[5]))
          at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
          for(i in 1:length(season_list)){
            seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
            seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
            if(use_mean==TRUE){
              typ='const'
              var_fit<-VAR(err_mat_fit[j,e,seas_fit,],p=lag,type=typ)
              coefs<-coef(var_fit)
              for(k in 1:length(coefs)){
                var_coef[i,k,]<-coefs[[k]][,'Estimate']}
              at[seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coef[i,,-c(dim(var_coef)[3])],lag,var_coef[i,,dim(var_coef)[3]],use_mean)}
            if(use_mean==FALSE){
              typ='none'
              var_fit<-VAR(err_mat_fit[j,e,seas_fit,],p=lag,type=typ)
              coefs<-coef(var_fit)
              for(k in 1:length(coefs)){
                var_coef[i,k,]<-coefs[[k]][,'Estimate']}
              at[seas_all,]<-var_resids_fun(err_mat[j,e,seas_all,],var_coef[i,,],lag,c(0),use_mean)}
          }
        return(list(at,var_coef))
        }
      
        for(e in 1:1:dim(err_mat)[1]){
          at_mat[j,e,,]<-out[[(e*2-1)]]
          coef_out<-out[[e*2]]
          #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
          coef_out[coef_out>=1]<-0
          coef_out[coef_out<=(-1)]<-0
          var_coefs[j,e,,,]<-coef_out
        }
      }
    }
    
    return(list(at_mat,var_coefs))
  }
  
  #5) Fit vanilla AR(p) model
  if(var_mod=='ar'){
    #output arrays
    at_mat <- array(NA,dim(err_mat))
    if(use_mean=='FALSE'){
      ar_coefs <- array(NA,c(n_sites,n_ens,length(season_list),n_lds,lag))}
    if(use_mean=='TRUE'){
      ar_coefs <- array(NA,c(n_sites,n_ens,length(season_list),n_lds,lag+1))}
    
    ar_fun<-function(x){
      out<-arima(x,order=c(lag,0,0),include.mean = use_mean,method='CSS')
      return(out)}
    
    if(parallel=='FALSE'){
      for(j in 1:n_sites){
        for(e in 1:n_ens){
          for(i in 1:length(season_list)){
            seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
            seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
            #AR(P) fit
            ar_coefs[j,e,i,,]<-t(apply(err_mat_fit[j,e,seas_fit,],2,function(x){out<-ar_fun(x);return(out$coef)}))
            for(k in 1:n_lds){
            at_mat[j,e,seas_all,k]<-arima(err_mat[j,e,seas_all,k],order=c(lag,0,0),fixed=ar_coefs[j,e,i,k,],include.mean = use_mean,method='CSS')$residuals}
          }
        }
      }
    }
    
    if(parallel=='TRUE'){
      for(j in 1:n_sites){
        out<-foreach(e = 1:n_ens,.combine='c',.packages = c('stats')) %dopar% {
          ar_coef<-array(NA,c(dim(ar_coefs)[3],dim(ar_coefs)[4],dim(ar_coefs)[5]))
          at<-array(NA,c(dim(at_mat)[3],dim(at_mat)[4]))
          for(i in 1:length(season_list)){
            seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
            seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
            #try robust least squares VAR fit
            ar_coef[i,,]<-t(apply(err_mat_fit[j,e,seas_fit,],2,function(x){out<-ar_fun(x);return(out$coef)}))
            for(k in 1:dim(err_mat)[3]){
              at[seas_all,k]<-arima(err_mat[j,e,seas_all,k],order=c(lag,0,0),fixed=ar_coef[i,k,],include.mean = use_mean,method='CSS')$residuals}
          }
          return(list(at,ar_coef))
        }
      
        for(e in 1:n_ens){
          at_mat[j,e,,]<-out[[(e*2-1)]]
          ar_coefs[j,e,,,]<-out[[e*2]]
        }
      }
    }
    return(list(at_mat,ar_coefs))
  }
  
  if(use_mpi==FALSE){
    stopCluster(cl = my.cluster)}
  if(use_mpi==TRUE){
    closeCluster(cl)
    Rmpi::mpi.quit()}
}

###########################################END################################