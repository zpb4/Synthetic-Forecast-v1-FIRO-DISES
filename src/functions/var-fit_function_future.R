
#library(fGarch) #required for SGED fit
#library(BigVAR)
#library(rmgarch) #required for varxfit
#library(vars)
library(bigtime)
library(doParallel)
library(future)
library(Rmpi)

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
var_fit<-function(err_mat,lag,var_mod,use_mean,seasonal,start_date,end_date,parallel){
  #parallelization code
  parallel::detectCores()

  n.cores <- parallel::detectCores()-2
  my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')

  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()

  #date index
  ixx_fit <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),by='day'),tz = "UTC")

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
    var_coefs <- array(NA,c(dim(err_mat)[1],length(season_list),dim(err_mat)[3],dim(err_mat)[3]*lag))}
  if(use_mean=='TRUE'){
    var_coefs <- array(NA,c(dim(err_mat)[1],length(season_list),dim(err_mat)[3],dim(err_mat)[3]*lag+1))}
  
  #1) 'bigtime' sparseVAR model
  if(var_mod=='svar'){
    
  if(parallel=='FALSE'){
    plan(cluster)
    out<-future(for(e in 1:dim(err_mat)[1]){
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        var_fit<-sparseVAR(err_mat[e,seas,],p=lag,VARpen = 'L1',selection = 'cv',check_std = F)
        var_resids<-var_resids_fun(err_mat[e,seas,],var_fit$Phihat,lag=lag,var_fit$phi0hat,use_mean)
        at_mat[e,seas,]<-var_resids
        if(use_mean=='FALSE'){
          var_coefs[e,i,,]<-var_fit$Phihat}
        if(use_mean=='TRUE'){
          var_coefs[e,i,,]<-cbind(var_fit$Phihat,var_fit$phi0hat)}
        var_coefs[var_coefs>1|var_coefs<(-1)]<-0
      }
    })
    at_mat<-value(out)
  }
  
  if(parallel=='TRUE'){
    out<-foreach(e = 1:dim(err_mat)[1],.combine='c',.packages = c('bigtime'),.export = c('var_resids_fun','lagmat_fun')) %dopar% {
      var_coef<-array(NA,c(dim(var_coefs)[2],dim(var_coefs)[3],dim(var_coefs)[4]))
      at<-array(NA,c(dim(at_mat)[2],dim(at_mat)[3]))
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        var_fit<-sparseVAR(err_mat[e,seas,],p=lag,VARpen = 'L1',selection = 'cv',check_std = F)
        var_resids<-var_resids_fun(err_mat[e,seas,],var_fit$Phihat,lag=lag,var_fit$phi0hat,use_mean)
        at[seas,]<-var_resids
        if(use_mean=='FALSE'){
          var_coef[i,,]<-var_fit$Phihat}
        if(use_mean=='TRUE'){
          var_coef[i,,]<-cbind(var_fit$Phihat,var_fit$phi0hat)}
      }
      return(list(at,var_coef))
    }
    
    for(e in 1:dim(err_mat)[1]){
      at_mat[e,,]<-out[[(e*2-1)]]
      var_coefs[e,,,]<-out[[e*2]]
    }
    var_coefs[var_coefs>1|var_coefs<(-1)]<-0
  }
  
  return(list(at_mat,var_coefs))
  }

  #2) Fit 'BigVAR' lasso penalized VAR model
  if(var_mod=='bvar'){

  if(parallel=='FALSE'){
    for(e in 1:dim(err_mat)[1]){
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        m1 = constructModel(err_mat[e,seas,], p = lag, struct = "Basic", gran = c(10, 10),IC = F,
        verbose = F, VARX = list(),separate_lambdas = F, rolling_oos=F,model.controls=list(intercept = use_mean, MN=F))
        m1_res = cv.BigVAR(m1)
        
        at_mat[e,seas,]<-rbind(err_mat[e,seas,][1:lag,],m1_res@resids)
        if(use_mean=='FALSE'){
          var_coefs[e,i,,]<-var_fit$Bcoef[,2:dim(var_fit$Bcoef)[2]]}
        if(use_mean=='TRUE'){
          var_coefs[e,i,,]<-cbind(var_fit$Bcoef[,2:dim(var_fit$Bcoef)[2]],var_fit$Bcoef[,1])}
        var_coefs[var_coefs>1|var_coefs<(-1)]<-0
      }
    }
  }
  
  if(parallel=='TRUE'){
    out<-foreach(e = 1:dim(err_mat)[1],.combine='c',.packages = c('bigtime'),.export = c('var_resids_fun','lagmat_fun')) %dopar% {
      var_coef<-array(NA,c(dim(var_coefs)[2],dim(var_coefs)[3],dim(var_coefs)[4]))
      at<-array(NA,c(dim(at_mat)[2],dim(at_mat)[3]))
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        m1 = constructModel(err_mat[e,seas,], p = lag, struct = "Basic", gran = c(10, 10),IC = F,
                            verbose = F, VARX = list(),separate_lambdas = F, rolling_oos=F,model.controls=list(intercept = use_mean, MN=F))
        m1_res = cv.BigVAR(m1)
        
        at[seas,]<-rbind(err_mat[e,seas,][1:lag,],m1_res@resids)
        if(use_mean=='FALSE'){
          var_coefs[e,i,,]<-var_fit$Bcoef[,2:dim(var_fit$Bcoef)[2]]}
        if(use_mean=='TRUE'){
          var_coefs[e,i,,]<-cbind(var_fit$Bcoef[,2:dim(var_fit$Bcoef)[2]],var_fit$Bcoef[,1])}
      }
      return(list(at,var_coef))
    }
    
    for(e in 1:dim(err_mat)[1]){
      at_mat[e,,]<-out[[(e*2-1)]]
      var_coefs[e,,,]<-out[[e*2]]
    }
    var_coefs[var_coefs>1|var_coefs<(-1)]<-0
  }
  
  return(list(at_mat,var_coefs))
  }

  #2) Fit rmgarch 'varx' robust least squares VAR model
  if(var_mod=='varx'){
    
  if(parallel=='FALSE'){
    for(e in 1:dim(err_mat)[1]){
      var_coef<-array(NA,c(dim(var_coefs)[2],dim(var_coefs)[3],dim(var_coefs)[4]))
      at<-array(NA,c(dim(at_mat)[2],dim(at_mat)[3]))
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        #try robust least squares VAR fit
        var_fit<-try(varxfit(err_mat[e,seas,],lag,constant=use_mean,robust = T),T)
        at[seas,]<-try(rbind(err_mat[e,seas,][1:lag,],var_fit$xresiduals),T)
        var_coef[i,,]<-try(var_fit$Bcoef,T)
        #if it doesnt work, try OLS VAR
        if(class(var_fit)=='try-error'){
          var_fit<-try(varxfit(err_mat[e,seas,],lag,constant=use_mean,robust = F),T)
          at[seas,]<-try(rbind(err_mat[e,seas,][1:lag,],var_fit$xresiduals),T)
          var_coef[i,,]<-try(var_fit$Bcoef,T)
          print(paste('ens',e,'mth',i,'trying OLS VAR'))
        }
        #as last resort, set VAR residuals to normalized errors and VAR coefficients all to zero
        if(class(var_fit)=='try-error' | anyNA(as.numeric(at[seas,]))==T | anyNA(as.numeric(var_coef[i,,]))==T){
          at[seas,]<-err_mat[e,seas,]
          var_coef[i,,]<-array(0,dim(var_coef[i,,]))
          print(paste('ens',e,'mth',i,'all zero VAR'))
        }
      }
      #correct missing or erroneous data as needed
      res_out<-at
      res_out[res_out==Inf|res_out==-Inf]<-0
      at_mat[e,,]<-apply(res_out,c(1:2),as.numeric)
      coef_out<-var_coef
      coef_out<-apply(coef_out,c(1:3),as.numeric)
      #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
      coef_out[coef_out>=1]<-0
      coef_out[coef_out<=(-1)]<-0
      var_coefs[e,,,]<-coef_out
    }
  }
  
  if(parallel=='TRUE'){
    out<-foreach(e = 1:1:dim(err_mat)[1],.combine='c',.packages = c('rmgarch')) %dopar% {
      var_coef<-array(NA,c(dim(var_coefs)[2],dim(var_coefs)[3],dim(var_coefs)[4]))
      at<-array(NA,c(dim(at_mat)[2],dim(at_mat)[3]))
      for(i in 1:length(season_list)){
        seas<-which(ixx_fit$mon%in%(season_list[[i]]-1))
        #try robust least squares VAR fit
        var_fit<-try(varxfit(err_mat[e,seas,],lag,constant=use_mean,robust = T),T)
        at[seas,]<-try(rbind(err_mat[e,seas,][1:lag,],var_fit$xresiduals),T)
        var_coef[i,,]<-try(var_fit$Bcoef,T)
        #if it doesnt work, try OLS VAR
        if(class(var_fit)=='try-error'){
          var_fit<-try(varxfit(err_mat[e,seas,],lag,constant=use_mean,robust = F),T)
          at[seas,]<-try(rbind(err_mat[e,seas,][1:lag,],var_fit$xresiduals),T)
          var_coef[i,,]<-try(var_fit$Bcoef,T)
          print(paste('ens',e,'mth',i,'trying OLS VAR'))
        }
        #as last resort, set VAR residuals to normalized errors and VAR coefficients all to zero
        if(class(var_fit)=='try-error' | anyNA(as.numeric(at[seas,]))==T | anyNA(as.numeric(var_coef[i,,]))==T){
          at[seas,]<-err_mat[e,seas,]
          var_coef[i,,]<-array(0,dim(var_coef[i,,]))
          print(paste('ens',e,'mth',i,'all zero VAR'))
        }
      }
      
      return(list(at,var_coef))
    }
    
    for(e in 1:ens_num){
      res_out<-out[[(e*2-1)]]
      res_out[res_out==Inf|res_out==-Inf]<-0
      at_mat[e,,]<-apply(res_out,c(1:2),as.numeric)
      coef_out<-out[[e*2]]
      coef_out<-apply(coef_out,c(1:3),as.numeric)
      #OLS VAR estimates in low flow months often unstable coefficients, replace them with zeroes
      coef_out[coef_out>=1]<-0
      coef_out[coef_out<=(-1)]<-0
      var_coefs[e,,,]<-coef_out
    }
  }
  
  return(list(at_mat,var_coefs))
  }

}

###########################################END################################