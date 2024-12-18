library(abind)
source('./src/functions/hybrid-loess-fit_fun.R')

#function to rearrange forward looking forecast to obs-synched forecast
fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}

logreg_pred<-function(x,coeffs){
  p <- exp(coeffs[1] +  x * coeffs[2]) / (1 + exp(coeffs[1] +  x * coeffs[2]))
  prd<-c()
  for(i in 1:length(x)){
    prd[i]<-rbinom(1,1,p)
  }
  return(prd)
}

#Function to estimate conditional expectation of forecasts given observations
#forecasts must be arranged in 'fwd look' format
#derived based on ensemble mean or median of forecasts
fit_logreg<-function(obs,forecast,seasonal,fit_start,fit_end,lv_out_yrs){
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'))
  
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
  
  #process obs and forecasts when only 1 site
  if(is.null(dim(obs)[2])==TRUE){
    obs_mat<-matrix(rep(obs,dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    obs_mat_fit<-obs_mat[ixx_obs%in%ixx_fit,]
    forecast_fit<-forecast[1,,-c(idx_lv_out),]
  }
  
  #concatenate across lead time dimension if more than one site
  if(dim(obs)[2]>1){
    obs_mat<-matrix(rep(obs[,1],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    obs_mat_fit<-obs_mat[ixx_obs%in%ixx_fit,]
    forecast_out<-forecast[1,,,]
    forecast_fit<-forecast_out[,ixx_hefs%in%ixx_fit,]
    for(i in 2:dim(obs)[2]){
      obs_sset<-matrix(rep(obs[,i],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
      obs_mat<-cbind(obs_mat,obs_sset)
      obs_fit<-obs_sset[ixx_obs%in%ixx_fit,]
      obs_mat_fit<-cbind(obs_mat_fit,obs_fit)}
    
    for(i in 2:dim(obs)[2]){
      forc<-forecast[i,,,]
      forc_fit<-forc[,ixx_hefs%in%ixx_fit,]
      forecast_fit<-abind(forecast_fit,forc_fit,along=3)
      forecast_out<-abind(forecast_out,forc,along=3)}
  }
  
  #determine seasonal subsets
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
  
  lreg_fit<-vector('list',length(season_list))
  lreg_fit_sub<-vector('list',dim(forecast_fit)[3])
  
  #fit hybrid LOESS model to the data (obs v forecasts)
  for(i in 1:length(season_list)){
    seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
    seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
    lreg_fit[[i]]<-lreg_fit_sub
    #fit by month and lead time across both sites
    for(k in 1:(dim(forecast_fit)[3])){
      forc_vec<-as.vector(t(forecast_fit[,seas_fit,k]))
      forc_vec[forc_vec>0]<-1
      obs_vec<-rep(obs_mat_fit[seas_fit,k],dim(forecast_fit)[1])
      lrfit=try(glm(forc_vec~obs_vec,family='binomial'))
      lr_out=try(logreg_pred(obs_mat[seas_all,k],lrfit$coefficients))
      lrfit2=
      lreg_fit[[i]][[k]]<-glm(forc_vec~obs_vec,family='binomial')$coefficients
      
      plot(1:775,obs_vec[776:(775*2)],typ='l')
      lines(1:775,forc_vec[776:(775*2)],col='green')
      
      #model <- logistic_reg(mixture = double(1), penalty = double(1)) %>%
        #set_engine("glmnet") %>%
        #set_mode("classification") %>%
        #fit(y ~ ., data = train)
      
      #forc_vec<-as.vector(t(forecast_out[,seas_all,k]))
      #forc_vec[forc_vec>0]<-1
      #obs_vec<-rep(obs_mat[seas_all,k],dim(forecast_out)[1])
      #lreg_fit[[i]][[k]]<-glm(forc_vec~obs_vec,family='binomial')$coefficients
    }
  }
  return(lreg_fit)
}

fit_zero_flow_mod<-function(obs,forecast,seasonal,fit_start,fit_end,lv_out_yrs){
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'))
  
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
  
  if(is.null(dim(obs)[2])==TRUE){obs <- matrix(obs,ncol=1)}
  n_sites = dim(obs)[2]
  n_ens = dim(forecast)[2]
  n_hefs = dim(forecast)[3]
  n_lds = dim(forecast)[4]
  
  
  #process obs and forecasts when only 1 site
  ##if(is.null(dim(obs)[2])==TRUE){
    ##obs_mat<-matrix(rep(obs,dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    ##obs_mat_fit<-obs_mat[ixx_obs%in%ixx_fit,]
    ##forecast_fit<-forecast[1,,-c(idx_lv_out),]
  ##}
  
  #concatenate across lead time dimension if more than one site
  ##if(dim(obs)[2]>1){
    ##obs_mat<-matrix(rep(obs[,1],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    ##obs_mat_fit<-obs_mat[-c(idx_lv_out),]
    ##forecast_out<-forecast[1,,,]
    ##forecast_fit<-forecast_out[,-c(idx_lv_out),]
    ##for(i in 2:dim(obs)[2]){
      ##obs_sset<-matrix(rep(obs[,i],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
      ##obs_mat<-cbind(obs_mat,obs_sset)
      ##obs_fit<-obs_sset[-c(idx_lv_out),]
      ##obs_mat_fit<-cbind(obs_mat_fit,obs_fit)}
    
    ##for(i in 2:dim(obs)[2]){
      ##forc<-forecast[i,,,]
      ##forc_fit<-forc[,-c(idx_lv_out),]
      ##forecast_fit<-abind(forecast_fit,forc_fit,along=3)
      ##forecast_out<-abind(forecast_out,forc,along=3)}
  ##}
  
  forecast_fit <- forecast[,,ixx_fit_all%in%ixx_fit,,drop=FALSE]
  obs_fit <- obs[ixx_fit_all%in%ixx_fit,,drop=FALSE]
  
  frac_fun<-function(x){out = sum(x)/length(x);return(out)}
  
  zero_mod_fit<-vector('list',n_sites)
  
  #fit hybrid LOESS model to the data (obs v forecasts)

  #fit by month and lead time across both sites
  for(j in 1:n_sites){
    zero_mod_fit[[j]] <- vector('list',n_lds)
    for(k in 1:n_lds){
      forc_vec<-forecast_fit[j,,,k]
      forc_vec[forc_vec>0]<-1
    
      forc_nzero_frac = apply(forc_vec,2,frac_fun)
    
      obs_vec<-obs_fit[,j]
    
      q = quantile(obs_vec,probs = seq(0, 1, 0.1))
      qvec = unique(q)
      qvec[length(qvec)]=Inf
      ct_idx = cut(obs_vec,breaks=qvec,include.lowest = T,labels=1:(length(qvec)-1))
      ct_idx = as.numeric(ct_idx)
    
      frac_vec<-c()
      for(i in 1:(length(qvec)-1)){
        ix = which(ct_idx==i)
        frac_vec[i]=mean(forc_nzero_frac[ix])
      }
      zero_mod_sub = list(qvec,frac_vec)
      zero_mod_fit[[j]][[k]] = zero_mod_sub
    }
  }
  
  return(zero_mod_fit)
}

fit_cexp<-function(obs,forecast,ctr,seasonal,fit_start,fit_end,lv_out_yrs,noise_reg=T){
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'))
  
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
  
  if(is.null(dim(obs)[2])==TRUE){obs <- matrix(obs,ncol=1)}
  n_sites = dim(obs)[2]
  n_ens = dim(forecast)[2]
  n_hefs = dim(forecast)[3]
  n_lds = dim(forecast)[4]
  
  #process obs and forecasts when only 1 site
  ##if(is.null(dim(obs)[2])==TRUE){
    ##obs_mat<-matrix(rep(obs,dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    ##obs_mat_fit<-obs_mat[-c(idx_lv_out),]
    ##forecast_fit<-forecast[1,,-c(idx_lv_out),]
  ##}

  #concatenate across lead time dimension if more than one site
  ##if(dim(obs)[2]>1){
    ##obs_mat<-matrix(rep(obs[,1],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
    ##obs_mat_fit<-obs_mat[-c(idx_lv_out),]
    ##forecast_out<-forecast[1,,,]
    ##forecast_fit<-forecast_out[,-c(idx_lv_out),]
    ##for(i in 2:dim(obs)[2]){
      ##obs_sset<-matrix(rep(obs[,i],dim(forecast)[4]),ncol=dim(forecast)[4],byrow=F)
      ##obs_mat<-cbind(obs_mat,obs_sset)
      ##obs_fit<-obs_sset[-c(idx_lv_out),]
      ##obs_mat_fit<-cbind(obs_mat_fit,obs_fit)}
  
    ##for(i in 2:dim(obs)[2]){
      ##forc<-forecast[i,,,]
      ##forc_fit<-forc[,-c(idx_lv_out),]
      ##forecast_fit<-abind(forecast_fit,forc_fit,along=3)
      ##forecast_out<-abind(forecast_out,forc,along=3)}
  ##}
  
  forecast_fit <- forecast[,,ixx_fit_all%in%ixx_fit,,drop=FALSE]
  forecast_out <- forecast[,,,,drop=FALSE]
  obs_fit <- obs[ixx_fit_all%in%ixx_fit,,drop=FALSE]
  obs_all <- obs
  
  #determine seasonal subsets
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

  #calculate ensemble ctr
  if(ctr=='median'){
    forc_ctr_mat<-apply(forecast_fit,c(1,3,4),median)
  }

  if(ctr=='mean'){
    forc_ctr_mat<-apply(forecast_fit,c(1,3,4),mean)
  }

  forc_cexp<-array(NA,c(n_sites,n_hefs,n_lds))

  #list to store hybrid LOESS model parameters
  cexp_fit<-vector('list',length(season_list))

  #fit hybrid LOESS model to the data (obs v forecasts)
  for(i in 1:length(season_list)){
    seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
    seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
    cexp_fit[[i]]<-vector('list',n_sites)
    #fit by month and lead time across both sites
    for(j in 1:n_sites){
      cexp_fit[[i]][[j]] <- vector('list',n_lds)
      for(k in 1:n_lds){
        #try symmetric (non-parametric) LOESS fit
        hy_fit<-try(hyb_loess_fit(obs_fit[seas_fit,j],forc_ctr_mat[j,seas_fit,k],0.75,2,'symmetric'),T)
        #if that doesn't work, try a gaussian fit
        if(class(hy_fit)=='try-error'){
          #print(paste('mth',k,'K=',i,'trying gaussian'))
          hy_fit<-try(hyb_loess_fit(obs_fit[seas_fit,j],forc_ctr_mat[j,seas_fit,k],0.75,2,'gaussian'),T)
        }
        #if that doesn't work reduce to a local linear model (degree 1) fit
        if(class(hy_fit)=='try-error'){
          #print(paste('mth',k,'K=',i,'trying deg1'))
          hy_fit<-try(hyb_loess_fit(obs_fit[seas_fit,j],forc_ctr_mat[j,seas_fit,k],0.75,1,'gaussian'),T)
        }
      cexp_fit[[i]][[j]][[k]]<-hy_fit
      cexp<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],obs_all[seas_all,j])
      forc_cexp[j,seas_all,k]<-cexp
      }
    }
  }

  #2b) calculate debiased errors as difference between conditional expectation and forecasts across lead times
  forc_err<-array(NA,dim(forecast_out))
  
  for(j in 1:n_sites){
    for(i in 1:n_ens){
      obs_mat<-matrix(rep(obs_all[,j],n_lds),ncol=n_lds,byrow=F)
      forc_cex<-forc_cexp[j,,]
      na_idx<-which(is.na(forc_cex)==T)
      forc_cex[na_idx]<-obs_mat[na_idx]
      neg_idx<-which(forc_cex<0)
      forc_cex[neg_idx]<-obs_mat[neg_idx]
      forc_err[j,i,,]<-forc_cex - forecast_out[j,i,,]
    }
  }

  if(noise_reg==T){
    for(j in 1:n_sites){
      obs_srt<-sort(obs[ixx_fit_all%in%ixx_fit,j])
      obs_lwr<-obs_srt[1:round(0.5*length(obs_srt))]
      ob_sdev<-sd(obs_lwr)
      forc_err_sset<-forc_err[j,,,]
      zero_idx<-which(forc_err_sset==0)
      nzero_idx<-which(forc_err_sset!=0)
      forc_err_sset[zero_idx]<-rnorm(length(zero_idx),sd=ob_sdev)
      forc_err_sset[nzero_idx]<- forc_err_sset[nzero_idx] + forc_err_sset[nzero_idx]*rnorm(length(nzero_idx),sd=0.1)
      forc_err[j,,,]<-forc_err_sset
    }
  }

return(list(cexp_fit,forc_cexp,forc_err))
}

#-------------------------------------------------------------------------------
#3) Fit linear heteroscedasticity model to debiased residuals

linear_hsked_fit<-function(err_mat,cexp_mat,obs,scale_ref,seasonal,min_quantile,fit_start,fit_end,lv_out_yrs){
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'))
  
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
  
  if(is.null(dim(obs)[2])==TRUE){obs <- matrix(obs,ncol=1)}
  obs_fit <- obs[ixx_fit_all%in%ixx_fit,,drop=FALSE]
  obs_all <- obs
  
  cexp_mat_fit<-cexp_mat[,ixx_fit_all%in%ixx_fit,,drop=FALSE] 
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
  
  lbs<-c(0,0) #lower bounds for intcpt, slope
  ubs<-c(max(abs(err_mat)),10) #upper bounds for intcpt, slope
  sts<-c(0,1) #starting parameters for intcpt, slope

  #list to store linear model parameters across months and lead times
  norm_fit<-vector('list',length(season_list))
  sd_arr<-array(NA,c(length(season_list),n_sites,n_lds))
  err_norm_mat<-array(NA,dim(err_mat))

  for(i in 1:length(season_list)){
    seas_all<-which(ixx_fit_all$mon%in%(season_list[[i]]-1))
    seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
    if(scale_ref=='cexp'){
      scale_mat_fit<-cexp_mat_fit[,seas_fit,,drop=FALSE]
      scale_mat<-cexp_mat[,seas_all,,drop=FALSE]
    }
    if(scale_ref=='obs'){
      scale_mat_fit<-obs_fit[seas_fit,,drop=FALSE]
      scale_mat<-obs_all[seas_all,,drop=FALSE]
    }
    std_errs<-abs(err_mat_fit[,,seas_fit,,drop=FALSE]) #absolute value of errors as proxy for stdev
    mean_std_errs<-apply(std_errs,c(1,3,4),mean) #take mean of above across ensembles to simplify prediction
    norm_fit[[i]]<-vector('list',n_sites)
    for(j in 1:n_sites){
      norm_fit[[i]][[j]]<-vector('list',n_lds)
      for(k in 1:n_lds){
        #------------------------
        #1. Fit linear model to abs value of errors
        if(scale_ref=='cexp'){
          scale_vec_fit<-scale_mat_fit[j,,k]
          scale_vec<-scale_mat[j,,k]}
        if(scale_ref=='obs'){
          scale_vec_fit<-scale_mat_fit[,k]
          scale_vec<-scale_mat[,k]}
        res<-mean_std_errs[j,,k]
        res_sort<-sort(res)
        lbs[1]<-mean(res_sort[1:round(min_quantile*length(seas_fit))]) #lower bound of intcpt inferred from lower 80th % of data, prevents extremely small normalization values
        ubs[1]<-max(res) #ensure intercept stays within range of the data
        sd_arr[i,j,k]<-mean(res) #global stdev as backup for linear model instability
    
        #try to fit linear model
        lfit<-try(optim(par=sts,lin_fit,x=scale_vec_fit,y=res,
                method = 'L-BFGS-B',lower = lbs,upper = ubs,
                control = list(fnscale=1,maxit=100000)),T)
        #if issues, use global stdev as intcpt with no heteroscedastic slope
        if(class(lfit)=='try-error'){linpars<-c(mean(res),0)}
        else{linpars<-lfit$par}
    
        norm_fit[[i]][[j]][[k]]<-linpars
      
        #-------------------------------
        #2. Calculate normalized errors
        norm_vec<-lin_mod(norm_fit[[i]][[j]][[k]][1],norm_fit[[i]][[j]][[k]][2],scale_vec)
      
        inp_resids<-err_mat[j,,seas_all,k]/matrix(rep(norm_vec,n_ens),nrow=n_ens,byrow=T)
        #erroneous normalized errors (e.g. NA, Inf due to division by 0) replaced with random gaussian (0,1) errors
        inp_resids[is.na(inp_resids)==T]<-rnorm(length(which(is.na(inp_resids)==T)))
        inp_resids[inp_resids==0]<-rnorm(length(which(inp_resids==0)))
        inp_resids[inp_resids==Inf | inp_resids==-Inf]<-rnorm(length(which(inp_resids==Inf | inp_resids==-Inf)))
        inp_resids[abs(inp_resids)>10*max(abs(err_mat[j,,seas_all,k]))]<-rnorm(length(which(abs(inp_resids)>10*max(abs(err_mat[j,,seas_all,k])))))
      
        err_norm_mat[j,,seas_all,k]<-inp_resids
      }
    }
  }
  return(list(norm_fit,sd_arr,err_norm_mat))
}

###########################################END################################