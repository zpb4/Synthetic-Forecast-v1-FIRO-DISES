
source('./src/functions/hybrid-loess-fit_fun.R')

#function to rearrange forward looking forecast to obs-synched forecast
fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}

#Function to estimate conditional expectation of forecasts given observations
#forecasts must be arranged in 'fwd look' format
#derived based on ensemble mean or median of forecasts
#

fit_cexp<-function(obs,forecast,ctr,seasonal,start_date,end_date,lv_out_yrs){
  
  idx<-seq(as.Date(start_date),as.Date(end_date),by='day')
  ixx_all <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),by='day'))
  idx_lv_out<-c()
  for(i in 1:length(lv_out_yrs)){
    idx_lv_out<-c(idx_lv_out,which(idx==paste(lv_out_yrs[i]-1,'-10-01',sep='')):which(idx==paste(lv_out_yrs[i],'-09-30',sep='')))
  }
  ixx_fit<-ixx_all[-c(idx_lv_out)]
  
  obs_mat<-matrix(rep(obs,dim(forecast)[3]),ncol=dim(forecast)[3],byrow=F)
  obs_mat_fit<-obs_mat[-c(idx_lv_out),]
  
  forecast_fit<-forecast[,-c(idx_lv_out),]
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
  forc_ctr_mat<-apply(forecast_fit,c(2,3),median)
}

if(ctr=='mean'){
  forc_ctr_mat<-apply(forecast_fit,c(2,3),mean)
}

forc_cexp<-array(NA,dim(obs_mat))

#list to store hybrid LOESS model parameters
cexp_fit<-vector('list',length(season_list))
cexp_fit_sub<-vector('list',dim(forecast)[3])

#fit hybrid LOESS model to the data (obs v forecasts)
for(i in 1:length(season_list)){
  seas_all<-which(ixx_all$mon%in%(season_list[[i]]-1))
  seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
  cexp_fit[[i]]<-cexp_fit_sub
    #fit by month and lead time across both sites
    for(k in 1:(dim(forecast)[3])){
      #try symmetric (non-parametric) LOESS fit
      hy_fit<-try(hyb_loess_fit(obs_mat_fit[seas_fit,k],forc_ctr_mat[seas_fit,k],0.75,2,'symmetric'),T)
      #if that doesn't work, try a gaussian fit
      if(class(hy_fit)=='try-error'){
        #print(paste('mth',k,'K=',i,'trying gaussian'))
        hy_fit<-try(hyb_loess_fit(obs_mat_fit[seas_fit,k],forc_ctr_mat[seas_fit,k],0.75,2,'gaussian'),T)
      }
      #if that doesn't work reduce to a local linear model (degree 1) fit
      if(class(hy_fit)=='try-error'){
        #print(paste('mth',k,'K=',i,'trying deg1'))
        hy_fit<-try(hyb_loess_fit(obs_mat_fit[seas_fit,k],forc_ctr_mat[seas_fit,k],0.75,1,'gaussian'),T)
      }
      cexp_fit[[i]][[k]]<-hy_fit
      cexp<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],obs_mat[seas_all,k])
      forc_cexp[seas_all,k]<-cexp
    }
}

#set any NAs or negative condition expectation values to the observed value
na_idx<-which(is.na(forc_cexp)==T)
neg_idx<-which(forc_cexp<0)

forc_cexp[na_idx]<-obs_mat[na_idx]
forc_cexp[neg_idx]<-obs_mat[neg_idx]

#2b) calculate debiased errors as difference between conditional expectation and forecasts across lead times
forc_err<-array(NA,dim(forecast))

for(i in 1:dim(forecast)[1]){
  forc_err[i,,]<-forc_cexp - forecast[i,,]
}

return(list(cexp_fit,forc_cexp,forc_err))
}

#-------------------------------------------------------------------------------
#3) Fit linear heteroscedasticity model to debiased residuals

linear_hsked_fit<-function(err_mat,cexp_mat,obs,scale_ref,seasonal,min_quantile,start_date,end_date,lv_out_yrs){

  idx<-seq(as.Date(start_date),as.Date(end_date),by='day')
  ixx_all <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),by='day'))
  idx_lv_out<-c()
  for(i in 1:length(lv_out_yrs)){
    idx_lv_out<-c(idx_lv_out,which(idx==paste(lv_out_yrs[i]-1,'-10-01',sep='')):which(idx==paste(lv_out_yrs[i],'-09-30',sep='')))
  }
  ixx_fit<-ixx_all[-c(idx_lv_out)]
  
  obs_mat<-matrix(rep(obs,dim(cexp_mat)[2]),ncol=dim(cexp_mat)[2],byrow=F)
  obs_mat_fit<-obs_mat[-c(idx_lv_out),] 
  
  cexp_mat_fit<-cexp_mat[-c(idx_lv_out),] 
  err_mat_fit<-err_mat[,-c(idx_lv_out),]
  
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
  norm_fit_sub<-vector('list',dim(cexp_mat)[2])
  sd_arr<-array(NA,c(length(season_list),dim(cexp_mat)[2]))
  err_norm_mat<-array(NA,dim(err_mat))

  for(i in 1:length(season_list)){
    seas_all<-which(ixx_all$mon%in%(season_list[[i]]-1))
    seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
    if(scale_ref=='cexp'){
      scale_mat_fit<-cexp_mat_fit[seas_fit,]
      scale_mat<-cexp_mat[seas_all,]
    }
    if(scale_ref=='obs'){
      scale_mat_fit<-obs_mat_fit[seas_fit,]
      scale_mat<-obs_mat[seas_all,]
    }
    std_errs<-abs(err_mat_fit[,seas_fit,]) #absolute value of errors as proxy for stdev
    mean_std_errs<-apply(std_errs,c(2,3),mean) #take mean of above across ensembles to simplify prediction
    norm_fit[[i]]<-norm_fit_sub
    
    for(k in 1:dim(cexp_mat)[2]){
      #------------------------
      #1. Fit linear model to abs value of errors
      scale_vec_fit<-scale_mat_fit[,k]
      scale_vec<-scale_mat[,k]
      res<-mean_std_errs[,k]
      res_sort<-sort(res)
      lbs[1]<-mean(res_sort[1:round(min_quantile*length(seas_fit))]) #lower bound of intcpt inferred from lower 80th % of data, prevents extremely small normalization values
      ubs[1]<-max(res) #ensure intercept stays within range of the data
      sd_arr[i,k]<-mean(res) #global stdev as backup for linear model instability
    
      #try to fit linear model
      lfit<-try(optim(par=sts,lin_fit,x=scale_vec_fit,y=res,
                method = 'L-BFGS-B',lower = lbs,upper = ubs,
                control = list(fnscale=1,maxit=100000)),T)
      #if issues, use global stdev as intcpt with no heteroscedastic slope
      if(class(lfit)=='try-error'){linpars<-c(mean(res),0)}
      else{linpars<-lfit$par}
    
      norm_fit[[i]][[k]]<-linpars
      
      #-------------------------------
      #2. Calculate normalized errors
      norm_vec<-lin_mod(norm_fit[[i]][[k]][1],norm_fit[[i]][[k]][2],scale_vec)
      
      inp_resids<-err_mat[,seas_all,k]/matrix(rep(norm_vec,dim(err_mat)[1]),nrow=dim(err_mat)[1],byrow=T)
      #erroneous normalized errors (e.g. NA, Inf due to division by 0) replaced with random gaussian (0,1) errors
      inp_resids[is.na(inp_resids)==T]<-rnorm(length(which(is.na(inp_resids)==T)))
      inp_resids[inp_resids==0]<-rnorm(length(which(inp_resids==0)))
      inp_resids[inp_resids==Inf | inp_resids==-Inf]<-rnorm(length(which(inp_resids==Inf | inp_resids==-Inf)))
      inp_resids[abs(inp_resids)>10*max(abs(err_mat[,seas_all,k]))]<-rnorm(length(which(abs(inp_resids)>10*max(abs(err_mat[,seas_all,k])))))
      
      err_norm_mat[,seas_all,k]<-inp_resids
    }
  }
  return(list(norm_fit,sd_arr,err_norm_mat))
}

###########################################END################################